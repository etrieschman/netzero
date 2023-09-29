import pandas as pd

# SUMMARIZE BY COUNTS
def summarize_id_counts_byyear(df, ids):
    ldf = (df.groupby('year')[ids]
        .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
    for id in ids:
        prev = set()
        counts = {}
        counts['n'], counts['n_new'], counts['n_drop'] = [], [], []
        for i, row in ldf.iterrows():
            curr = row[id]
            # Calculate overlap and new IDs
            counts['n'] += [len(curr)]
            counts['n_new'] += [len(curr.difference(prev))]
            counts['n_drop'] += [len(prev.difference(curr))]
            # Update prev_ids for the next iteration
            prev = curr

        # Add new columns to the DataFrame
        for k, v in counts.items():
            ldf[f'{id}_{k}'] = v
        ldf = ldf.drop(columns=[id])
    return ldf


# %%
# SUMMARIZE MERGE HELPER FUNCTION
def summarize_merge(df, left_on, right_on, r_hasyear):
    summ = pd.DataFrame([])
    for y in df.year.drop_duplicates().sort_values().dropna().values:
        s = {}
        s['year'] = int(y)
        s['left'] = [len(df.loc[(df.year == y), left_on].drop_duplicates())]
        s['left_nright'] = [len(df.loc[(df.year == y) & df[right_on].isna().all(axis=1), 
                                    left_on].drop_duplicates())]
        s['overlap'] = [len(df.loc[(df.year == y) & df[left_on].notna().all(axis=1) & df[right_on].notna().all(axis=1),
                                    left_on + right_on].drop_duplicates())]
        if r_hasyear:
            y_condition = (df.year == y)
        else:
            y_condition = True
        s['right_nleft'] = [len(df.loc[y_condition & df[left_on].isna().all(axis=1) & df[right_on].notna().all(axis=1), right_on].drop_duplicates())]
        s['right'] = [len(df.loc[y_condition & df[right_on].notna().all(axis=1), right_on].drop_duplicates())]
        summ = pd.concat([summ, pd.DataFrame(s)], ignore_index=True)
    return summ


# helper summary function for in-epa merge
def summarize_plant_inepa(df, group):
    summ_plant = (df.groupby(['utility_id', 'plant_code', 'state_plant', 'entity_type', 'year', 'plant_in_epa'])
                  .agg({'nameplate_capacity_mw':'sum', 'key_gen':'nunique', 'not_renewable':'mean', 'gen_age_yrs':'mean'})
                  .reset_index())
    group_cols = ['year', 'plant_in_epa']
    if group is not None:
        group_cols += [group] 
    summ = (summ_plant.groupby(group_cols)
            .agg({'plant_code':['nunique'], 'key_gen':['sum', 'mean', 'std'], 'nameplate_capacity_mw':['sum', 'mean', 'std'],
                   'not_renewable':['mean', 'std'], 'gen_age_yrs':['mean', 'std']})
            .reset_index())
    summ.rename(columns={group:'value', 'not_renewable':'pct_gens_not_renewable', 'gen_age_yrs':'mean_gen_age_yrs', 
                         'plant_code':'plants', 'key_gen':'generators'}, inplace=True)
    if group is None:
        summ['value'] = 'total'
    return summ


# helper summary function
def summarize_gen_inepa(df, group):
    group_cols = ['year', 'gen_in_epa']
    if group is not None:
        group_cols += [group]
    summ = (df.groupby(group_cols)
                .agg({'plant_code':'nunique', 'key_gen':['nunique'], 'nameplate_capacity_mw':['sum', 'mean', 'std'], 'gen_age_yrs':['mean', 'std']})
                .reset_index())
    summ.rename(columns={group:'value', 'key_gen':'generators', 'plant_code':'plants'}, inplace=True)
    if group is None:
        summ['value'] = 'total'
    return summ
