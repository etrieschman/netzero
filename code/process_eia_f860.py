# %%
import pandas as pd
import numpy as np
from tqdm import tqdm

from utils import PATH_EIA, PATH_PROCESSED
from utils import readin_eia

YR_START, YR_END = 2018, 2021

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 25)

# %%
# UTILITY DATA
yr_start = 2018
vars_keep = ['utility_id', 'utility_name', 'city', 'state', 'zip', 'entity_type']
udf = readin_eia(YR_START, YR_END, f'{PATH_EIA}f860', '1___Utility_Y', vars_keep, 
                 readin_params={'header':1})
udf['utility_id'] = pd.to_numeric(udf.utility_id).astype('Int64')
udf.to_csv(PATH_PROCESSED + 'eia_f860_utility.csv', index=False)

# %%
# PLANT DATA
vars_keep = ['utility_id', 'plant_code', 'plant_name', 'latitude',
       'longitude', 'primary_purpose_naics_code']
pdf = readin_eia(YR_START, YR_END, f'{PATH_EIA}f860', '2___Plant_Y', vars_keep, 
                 readin_params={'header':1})
pdf['utility_id'] = pd.to_numeric(pdf.utility_id).astype('Int64')
pdf['plant_code'] = pd.to_numeric(pdf.plant_code).astype('Int64')
pdf.to_csv(PATH_PROCESSED + 'eia_f860_plant.csv', index=False)

# %%
# GENERATOR DATA
sheets = ['Operable', 'Retired and Canceled']
gen_tags = {'1_Generator':{'sheets':sheets + ['Proposed']}, 
            '2_Wind':{'sheets':sheets}, 
            '3_Solar':{'sheets':sheets}, 
            '4_Energy_Storage':{'sheets':sheets}, 
            '5_Multifuel':{'sheets':sheets + ['Proposed']}}
vars_date = ['operating_month', 'operating_year',
             'current_month', 'current_year',
             'planned_retirement_month', 'planned_retirement_year',
             'retirement_month', 'retirement_year']
vars_keep = ['utility_id', 'plant_code', 'generator_id', 'status', 'ownership', 'sector', 
             'technology', 'prime_mover', 'nameplate_capacity_mw'] + vars_date
gdf = pd.DataFrame({})
for g in gen_tags.keys():
    print('Generator type:', g)
    for s in gen_tags[g]['sheets']:
        df = readin_eia(YR_START, YR_END, f'{PATH_EIA}f860', f'3_{g}_Y', vars_keep, 
                        readin_params={'header':1, 'sheet_name':s})
        df['sheet_gen'] = s.lower().replace(' ', '_')
        df['utility_id'] = pd.to_numeric(df.utility_id, errors='coerce').astype('Int64')
        df['plant_code'] = pd.to_numeric(df.plant_code, errors='coerce').astype('Int64')
        df['generator_id'] = df.generator_id.astype(str, errors='ignore')
        df['nameplate_capacity_mw'] = pd.to_numeric(df.nameplate_capacity_mw, errors='coerce').astype('Float64')
        for var in df.columns.intersection(vars_date):
            df[var] = pd.to_numeric(df[var], errors='coerce').astype('Int64')
        df = df.loc[df.utility_id.notna()]
        gdf = pd.concat([df, gdf], axis=0, ignore_index=True)
gdf.to_csv(PATH_PROCESSED + 'eia_f860_generator.csv', index=False)

# %%
# OWNER DATA
vars_keep = ['utility_id', 'plant_code', 'generator_id', 
             'ownership_id', 'status', 'owner_name',
             'owner_city', 'owner_state', 'owner_zip', 'percent_owned']
odf = readin_eia(YR_START, YR_END, f'{PATH_EIA}f860', '4___Owner_y', vars_keep,
                 readin_params={'header':1})
odf['percent_owned'] = pd.to_numeric(odf.percent_owned.astype(str).str.strip())
odf.to_csv(PATH_PROCESSED + 'eia_f860_ownership.csv', index=False)





# %%
# SUMMARIZE BY COUNTS
def summarize_id_counts(df, ids):
    for id in ids:
        prev = set()
        counts = {}
        counts['n'], counts['n_new'], counts['n_drop'] = [], [], []
        for i, row in df.iterrows():
            curr = row[id]
            # Calculate overlap and new IDs
            counts['n'] += [len(curr)]
            counts['n_new'] += [len(curr.difference(prev))]
            counts['n_drop'] += [len(prev.difference(curr))]
            # Update prev_ids for the next iteration
            prev = curr

        # Add new columns to the DataFrame
        for k, v in counts.items():
            df[f'{id}_{k}'] = v
        df = df.drop(columns=[id])
    return df

# %%
# GENERATE UTILITY SUMMARY
# utility data
ludf = (udf.groupby('year')[['utility_id']]
           .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
ludf = ludf.rename(columns={'utility_id':'uid'})
print('Utility dataset:')
summarize_id_counts(ludf.copy(), ['uid'])

# %%
# GENERATE PLANT SUMMARY
# plant data
pdf['pid'] = pdf.utility_id.astype(str) + '.' + pdf.plant_code.astype(str)
lpdf = (pdf.groupby('year')[['utility_id', 'pid']]
        .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
lpdf = lpdf.rename(columns={'utility_id':'uid'})
print('Plant dataset:')
summarize_id_counts(lpdf.copy(), ['uid', 'pid'])

# %%
# GENERATE GENERATOR SUMMARY
# unit data
gdf['pid'] = gdf.utility_id.astype(str) + '.' + gdf.plant_code.astype(str)
gdf['gid'] = gdf.pid + '.' + gdf.generator_id
lgdf = (gdf.groupby(['year'])[['utility_id', 'pid', 'gid']]
        .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
lgdf = lgdf.rename(columns={'utility_id':'uid'})
print('Generator dataset:')
summarize_id_counts(lgdf.copy(), ['uid', 'pid', 'gid'])

# %%
