# %%
# INITIALIZE
import pandas as pd
import numpy as np
import os
import re

from utils import PATH_DATA, PATH_PROCESSED

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 25)

# NOTE: continuity across time looks good. For each year, let's look at overlap


# %%
# READIN PROCESSED DATA
# make dictionary to automate summaries
dfs = {
    'epa_emi':{'merge_from':['facility_id', 'unit_id'], 'merge_to':['camd_plant_id', 'camd_unit_id']},
    'epa_fac':{'merge_from':['facility_id', 'unit_id'], 'merge_to':['camd_plant_id', 'camd_unit_id']},
    'eia_own':{'merge_from':['plant_code', 'generator_id'], 'merge_to':['eia_plant_id', 'eia_generator_id']},
    'eia_pla':{'merge_from':['plant_code'], 'merge_to':['eia_plant_id']},
    'eia_gen':{'merge_from':['plant_code', 'generator_id'], 'merge_to':['eia_plant_id', 'eia_generator_id']},
    'eia_uop':{'merge_from':None},
    'eia_uti':{'merge_from':None}}
files = [f for f in os.listdir(PATH_PROCESSED)]
for f in files:
    print(f)
    fk = re.sub(r'f.{3}_|\.csv', '', f)[:7]
    dfs[fk]['data'] = pd.read_csv(PATH_PROCESSED + f)
    if dfs[fk]['merge_from'] is None:
        continue
    # dfs[fk]['data'][dfs[fk]['merge_from']] = dfs[fk]['data'][dfs[fk]['merge_from']].astype(int)
dfs.keys()

# %%
# READIN CROSSWALK
xwalk = pd.read_csv(PATH_DATA + 'epa_eia_crosswalk.csv')
xwalk.columns = xwalk.columns.str.lower()
# drop CAMD excluded rows
xwalk = xwalk.loc[xwalk.match_type_gen != 'Manual CAMD Excluded']
xwalk['eia_plant_id'] = xwalk['eia_plant_id'].astype('Int64')

# %%
# SUMMARIZE OVERLAP WITH CROSSWALK
vars_id = ['camd_plant_id', 'camd_unit_id', 'eia_plant_id', 'eia_generator_id']
summ = pd.DataFrame({})
for k in dfs.keys():
    if dfs[k]['merge_from'] is None:
        continue
    df = dfs[k]['data']
    m = (pd.merge(left=df, right=xwalk, how='outer', left_on=dfs[k]['merge_from'], right_on=dfs[k]['merge_to']))    
    s = {}
    s['frame'] = [k]
    for y in m.year.drop_duplicates().sort_values().dropna().values:
        s['year'] = y
        s['left'] = [len(df.loc[df.year == y, dfs[k]['merge_from']].drop_duplicates())]
        s['left_nright'] = [len(m.loc[(m.year == y) & m[dfs[k]['merge_from']].notna().all(axis=1) & m[dfs[k]['merge_to']].isna().all(axis=1), dfs[k]['merge_from']].drop_duplicates())]
        s['overlap'] = [m.loc[m.year == y, dfs[k]['merge_from'] + dfs[k]['merge_to'] + vars_id].drop_duplicates().notna().all(axis=1).sum()]
        s['right_nleft'] = [len(m.loc[m[dfs[k]['merge_from']].isna().all(axis=1), vars_id].drop_duplicates())]
        s['right'] = [len(xwalk[vars_id].drop_duplicates())]
        summ = pd.concat([summ, pd.DataFrame(s)], ignore_index=True)

summ


# %%
