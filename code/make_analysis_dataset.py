# %%
import numpy as np
import pandas as pd

from utils import PATH_DATA, PATH_PROCESSED

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 15)

# %%
# IMPORT
eia = pd.read_parquet(PATH_PROCESSED + 'eia_final.parquet')
epa = pd.read_csv(PATH_PROCESSED + 'epa_emissions.csv')
xwalk = pd.read_csv(PATH_DATA + 'epa_eia_crosswalk.csv')
eia['key_gen'] = eia.plant_code.astype(str) + '_' + eia.generator_id
eia['key_own'] = eia.key_gen + '_' + eia.ownership_id.astype(str)

# %%
# CLEAN CROSSWALK
xwalk.columns = xwalk.columns.str.lower()
nrows_xwalk_raw = len(xwalk)
xwalk = xwalk.loc[~xwalk.match_type_gen.isin(
    ['Manual CAMD Excluded', 'CAMD Unmatched'])]
print('CAMD units dropped, either excluded or unmatched:', nrows_xwalk_raw - len(xwalk))
xwalk = xwalk.astype({'eia_plant_id':'Int64', 'camd_plant_id':'Int64'})
# DECISION drop boiler ID and epa generator ID and dedup on remaining crosswalk
keepcols_xwalk = ['camd_plant_id', 'camd_unit_id',
       'eia_plant_id', 'eia_generator_id']
nrows_xwalk_boiler = len(xwalk)
xwalk = xwalk[keepcols_xwalk].drop_duplicates()
print('Crosswalk rows dropped for unique EIA boilers and EIA generators:', nrows_xwalk_boiler - len(xwalk))

# %%
# SUMMARIZE UNIT-TO-GENERATOR MAPPINGS
xwalk_summ = xwalk.copy()
xwalk_summ['n_eia_plants_pepa'] = xwalk.groupby(['camd_plant_id'])['eia_plant_id'].transform('nunique')
xwalk_summ['n_epa_plants_peia'] = xwalk.groupby(['eia_plant_id'])['camd_plant_id'].transform('nunique')
display(xwalk_summ[['camd_plant_id', 'eia_plant_id', 'n_eia_plants_pepa']].drop_duplicates()['n_eia_plants_pepa'].value_counts())
display(xwalk_summ[['camd_plant_id', 'eia_plant_id', 'n_epa_plants_peia']].drop_duplicates()['n_epa_plants_peia'].value_counts())
# NOTE: Anecdotally, we get multiple EIA plants to an EPA plant when an old plant goes offline and is replaced by a new plant or new owner
# display(xwalk.loc[xwalk.groupby(['camd_plant_id'])['eia_plant_id'].transform('nunique') > 1, ['camd_plant_id', 'eia_plant_id']].drop_duplicates())
# eia.loc[eia.plant_code.isin([302, 59002])]
# epa.loc[epa.facility_id == 302]
xwalk_summ['n_eia_gen_pepa'] = xwalk.groupby(['camd_plant_id', 'camd_unit_id'])['eia_generator_id'].transform('nunique')
xwalk_summ['n_epa_unit_peia'] = xwalk.groupby(['eia_plant_id', 'eia_generator_id'])['camd_unit_id'].transform('nunique')
xwalk_summ[['n_eia_gen_pepa', 'n_epa_unit_peia']].value_counts().reset_index()

# %%
# DECISION: BASED ON THESE FINDINGS, PROPOSE MAPPING EPA TO EIA AT THE FACILITY LEVEL
xwalk = xwalk['']


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


# %%
# UNDERSTAND MERGE EIA ON CROSSWALK
eia_on = ['plant_code', 'generator_id']
xwalk_eia_on = ['eia_plant_id', 'eia_generator_id'] 
m_eia = pd.merge(left=eia, right=xwalk, how='outer', 
             left_on=eia_on, 
             right_on=xwalk_eia_on)
summarize_merge(m_eia, eia_on, xwalk_eia_on, r_hasyear=False)

# %%
# UNDERSTAND MERGE EPA ON CROSSWALK
# DECISION: COLLAPSE TO EPA FACILITY
epa_col = epa.groupby(['facility_id', 'year'])[['co2_mass_short_tons', 'gross_load_mwh']].sum().reset_index()
epa_on = ['facility_id', 'unit_id']
xwalk_epa_on = ['camd_plant_id', 'camd_unit_id']
m_epa = pd.merge(left=epa_col, right=xwalk, how='outer',
              left_on=epa_on, right_on=xwalk_epa_on)
m_epa['year'] = m_epa.year.astype('Int64')
print('Merge summary')
display(summarize_merge(m_epa, epa_on, xwalk_epa_on, r_hasyear=False))


# %% 
# AGGREGATE EPA TO EIA GENERATOR
vars_group = ['facility_id', 'unit_id', 'year', 
              'camd_plant_id', 'eia_plant_id', 'eia_generator_id']
vars_col = ['sum_of_the_operating_time', 'gross_load_mwh', 'steam_load_1000_lb',
    'so2_mass_short_tons', 'co2_mass_short_tons', 'nox_mass_short_tons',
    'heat_input_mmbtu']
m_epa_col = m_epa.groupby(vars_group)[vars_col].sum().reset_index()
m = pd.merge(left=eia, right=m_epa_col, how='outer',
             left_on=eia_on+['year'], right_on=xwalk_eia_on+['year'])




summarize_merge(m, eia_on, xwalk_eia_on, r_hasyear=True)

# %%
m_epa['count_unique'] = m_epa.groupby(['facility_id'])['eia_plant_id'].transform('nunique')
m_epa.loc[m_epa.count_unique > 1].sort_values(['facility_id', 'eia_plant_id'])
# %%
m_epa.groupby(['facility_id', 'unit_id'])[['camd_generator_id', 'eia_generator_id']].nunique()

# %%
