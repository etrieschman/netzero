# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, re

from utils import PATH_PROCESSED

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 50)

# %%
# READIN DATA
generator = pd.read_csv(PATH_PROCESSED + 'eia_f860_generator.csv',
                        dtype={'nameplate_capacity_mw':object,
                               'ownership': str})
plant = pd.read_csv(PATH_PROCESSED + 'eia_f860_plant.csv')
utility = pd.read_csv(PATH_PROCESSED + 'eia_f860_utility.csv')
owner = pd.read_csv(PATH_PROCESSED + 'eia_f860_ownership.csv')

# %%
# Q: DOES OWNERSHIP ADD UP? No
print(owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'])['percent_owned'].sum().value_counts())
# DECISIONS:
# 1. Set all missing ownership to 100%
# 2. Manual fix upon inspection of ownership ID 7239 (should be 1., not 0.1)
# 3. Could not figure out what is happening with the Conemaugh Hydro plant (52173). 
#    Looks like possible consolidation of ownership, or possibly omitted reporting. Leaving as is
# 4. Decided not to change Long Beach Generation since it's a retired plant
owner.loc[owner.percent_owned.isna(), 'percent_owned'] = 1.
owner.loc[owner.ownership_id == 7239, 'percent_owned'] = 1.
owner['pct_ownership_total'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('sum').round(5)
owner['owner_count'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('count')
display(owner.loc[owner.pct_ownership_total != 1.])
owner = owner.drop(columns=['pct_ownership_total', 'owner_count'])

# %%
# MERGE
print('generator ID-years:\t', len(generator))
print('ownership ID-years:\t', len(owner))
m_go = pd.merge(left=generator, right=owner, how='outer',
         on=['utility_id', 'plant_code', 'generator_id', 'year', 'status'])
print('CHECK: merged ID-years:\t', len(m_go))
print('utility ID-years:\t', len(utility))
m_gou = pd.merge(left=m_go, right=utility, on=['utility_id', 'year'])
print('CHECK: merged ID-years:\t', len(m_gou))
assert len(m_go) == len(m_gou)
print('plant ID-years:\t\t', len(plant))
m = pd.merge(left=m_gou, right=plant, on=['utility_id', 'plant_code', 'year'])
print('CHECK: merged ID-years\t', len(m))
assert len(m_gou) == len(m)

# %% 
# UPDATE COLUMNS: OWNERSHIP TYPE
m_cln = m.copy()
update_mask = m_cln.ownership_id.isna().values
update_to_cols = ['owner_name', 'owner_city', 'owner_state', 'owner_zip', 'ownership_id']
update_from_cols = ['utility_name', 'city', 'state', 'zip']
for col_f, col_t in zip(update_to_cols, update_from_cols + ['utility_id']):
    m_cln.loc[update_mask, col_f] = m_cln.loc[update_mask, col_t]
m_cln.loc[update_mask, 'ownership'] = 'S' # Single ownership by respondent
m_cln.loc[~update_mask & (m_cln.percent_owned == 1.), 'ownership'] = 'W' # Wholly owned by an entity other than respondent
m_cln.loc[~update_mask & (m_cln.percent_owned != 1.), 'ownership'] = 'J' # Jointly owned with another entity
m_cln.loc[update_mask, 'percent_owned'] = 1.
# m_cln = m_cln.drop(columns=update_from_cols)

# %% 
# UPDATE COLUMNS: OPERATION STATUS
mask_op = (m_cln.status == 'OP')
mask_pl = (m_cln.status.isin(('TS', 'P', 'L', 'T', 'U', 'V')))
mask_sb = (m_cln.status.isin(('SB', 'OA', 'OS')))
mask_ds = (m_cln.status.isin(('RE', 'CN', 'IP', 'OT')))
m_cln.loc[mask_op, 'status_simp'] = 'operating'
m_cln.loc[mask_pl, 'status_simp'] = 'planned'
m_cln.loc[mask_sb, 'status_simp'] = 'standby'
m_cln.loc[mask_ds, 'status_simp'] = 'discontinued'
# OPERATION DATES
m_cln.loc[m_cln.operating_month > 12, 'operating_month'] = 1
m_cln.loc[mask_op | mask_sb | mask_ds , 'dt_operation_start'] = (
    pd.to_datetime(m_cln.loc[mask_op | mask_sb | mask_ds, 'operating_year'].astype('Int64').astype(str) + 
                   m_cln.loc[mask_op | mask_sb | mask_ds, 'operating_month'].astype('Int64').astype(str), format='%Y%m', errors='coerce'))
m_cln.loc[mask_op | mask_sb, 'dt_operation_end'] = (
    pd.to_datetime(m_cln.loc[mask_op | mask_sb, 'planned_retirement_year'].astype('Int64').astype(str) + 
                   m_cln.loc[mask_op | mask_sb, 'planned_retirement_month'].astype('Int64').astype(str), format='%Y%m', errors='coerce'))
m_cln.loc[mask_pl, 'dt_operation_start'] = (
    pd.to_datetime(m_cln.loc[mask_pl, 'current_year'].astype('Int64').astype(str) + 
                   m_cln.loc[mask_pl, 'current_month'].astype('Int64').astype(str), format='%Y%m'))
m_cln.loc[mask_ds, 'dt_operation_end'] = (
    pd.to_datetime(m_cln.loc[mask_ds, 'retirement_year'].astype('Int64').astype(str) + 
                   m_cln.loc[mask_ds, 'retirement_month'].astype('Int64').astype(str), format='%Y%m', errors='coerce'))
# SUMMARIZE
summ_status = (m_cln.groupby('status', dropna=False)
 [['utility_id', 'current_month', 'current_year', 'operating_month', 'operating_year', 
   'retirement_month', 'retirement_year', 'planned_retirement_month', 'planned_retirement_year']]
 .agg(lambda x: (x.notna() & (x != ' ')).sum()))
summ_status_simp = (m_cln.groupby('status_simp', dropna=False)
 [['utility_id', 'dt_operation_start', 'dt_operation_end']]
 .agg(lambda x: x.notna().sum()))
m_cln = m_cln.drop(columns=['operating_month', 'operating_year', 
                    'planned_retirement_month', 'planned_retirement_year', 
                    'current_month', 'current_year', 'retirement_month', 'retirement_year'])
display(summ_status)
display(summ_status_simp)


#%%
# REORDER COLUMNS

m_final = m_cln[[
    'year', 'ownership_id', 'utility_id', 'plant_code', 'generator_id',
    'ownership', 'owner_name', 'owner_city', 'owner_state', 'owner_zip', 'percent_owned',
    'entity_type', 'primary_purpose_naics_code', 'sector', 'latitude', 'longitude',
    'technology', 'prime_mover', 'nameplate_capacity_mw', 'status', 'status_simp',
    'dt_operation_start', 'dt_operation_end']]
m_final

generator.loc[(generator.plant_code == 492)].sort_values(['generator_id', 'year'], ascending=False)
# %%
m_final
# %%
