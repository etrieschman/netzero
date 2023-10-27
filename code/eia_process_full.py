# %%
import numpy as np
import pandas as pd

from utils import PATH_PROCESSED

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 50)

# %%
# READIN DATA
generator = pd.read_parquet(PATH_PROCESSED + 'eia_f860_generator.parquet')
plant = pd.read_parquet(PATH_PROCESSED + 'eia_f860_plant.parquet')
utility = pd.read_parquet(PATH_PROCESSED + 'eia_f860_utility.parquet')
owner = pd.read_parquet(PATH_PROCESSED + 'eia_f860_ownership.parquet')

# %%
# # Q: DOES OWNERSHIP ADD UP? No
# # DECISIONS:
# # 0. Divide by 100 in settings where decimal is off
# # 1. If owner_count is 1 and percent ownership is "close" to 1, nudge to 1
# # 2. Set all missing ownership to 100% if only 1 owner; 
# #    Set missing ownership to 0% if missing, but sum of ownership is 1
# # 3. TODO: 2010 is a messy year that needs to be reconciled
# owner['pct_ownership_total'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('sum').round(5)
# owner['owner_count'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('count')
# owner.loc[owner.pct_ownership_total == 100., 'percent_owned'] /= 100.
# owner.loc[(owner.owner_count == 1) & 
#     (np.abs(owner.pct_ownership_total - 1) <= 0.02), 'percent_owned'] = 1.
# owner.loc[(owner.percent_owned.isna()) & (owner.owner_count == 0), 'percent_owned'] = 1.
# owner.loc[(owner.percent_owned.isna()) & (owner.owner_count > 0), 'percent_owned'] = 0.
# # print(owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'])['percent_owned'].sum().value_counts().sort_index())
# owner['pct_ownership_total'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('sum').round(5)
# owner['owner_count'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('count')
# print("Count of generators where ownership doesn't add up after cleaning, by year:")
# display(owner.loc[np.abs(owner.pct_ownership_total - 1.) >= 0.02].groupby('year')['generator_id'].agg(['count', 'nunique']))
# owner = owner.drop(columns=['pct_ownership_total', 'owner_count'])

# %%
# MERGE GENERATORS <--> OWNERS
print('generator ID-years:\t', len(generator))
print('ownership ID-years:\t', len(owner))
# NOTE: drop utility ID from ownership data since it's messed up in 2010 
#       and not needed since generators are unique on plant+gen ids
m_go = pd.merge(left=generator, right=owner.drop(columns=['utility_id', 'status']), how='outer',
         on=['plant_code', 'generator_id', 'year'])
print('CHECK: merged ID-years:\t', len(m_go))

# %%
# MERGE UTILITY --> (GEN/OWNER)
utility['merge_key'] = utility.utility_id.astype(str) + "_" + utility.year.astype(str)
m_go['merge_key'] = m_go.utility_id.astype(str) + '_' + m_go.year.astype(str)
print('utility ID in generator, but not in utility:')
# m_go.loc[~m_go.merge_key.isin(utility.merge_key.drop_duplicates())].sort_values(['utility_id', 'plant_code', 'generator_id', 'year'])
display(m_go.loc[~m_go.merge_key.isin(utility.merge_key.drop_duplicates())].groupby('year')['utility_id'].count())
# NOTE: It is the usual suspect, 2010. Ignore the fact that ~150 utilities apear in generator dataset, but not in utility
print('utility ID not in merge dataset')
display(utility.loc[~utility.utility_id.isin(m_go.utility_id.drop_duplicates())].groupby('year')['utility_id'].count())
# FINDING: Looks like a chunk of utilities exist in '06-'09 that don't have associated generators
# FINDING: Looks like many of the utility names are wind/solar
# DECISION: Left merge onto gen/owner dataset
print('utility ID-years:\t', len(utility))
m_gou = pd.merge(left=m_go, right=utility, how='left', on=['utility_id', 'year'])
print('CHECK: merged ID-years:\t', len(m_gou))
assert(len(m_gou) == len(m_go))


# %%
# MERGE PLANT --> (UTILITY/GEN/OWNER)
plant['merge_key'] = plant.plant_code.astype(str) + "_" + plant.year.astype(str)
m_gou['merge_key'] = m_gou.plant_code.astype(str) + '_' + m_gou.year.astype(str)
print('plant ID in merge dataset, but not in plant:')
display(m_gou.loc[~m_gou.merge_key.isin(plant.merge_key.drop_duplicates())].groupby('year')['utility_id'].count())
# FINDING: these are three proposed solar plants that only appear in 2016. Fine to disregard
print('plant ID not in merge dataset')
display(plant.loc[~plant.merge_key.isin(m_gou.merge_key.drop_duplicates())].groupby('year')['merge_key'].count())
# NOTE: Handful of plant IDs without associated generators or owners or utilities
#       Looks to be moslty an issue in earlier years, especially 2010, the usual suspect
# DECISION: Inner join on plant id and year!
plant.drop(columns='merge_key', inplace=True)
m_gou.drop(columns='merge_key', inplace=True)
print('plant ID-years:\t\t', len(plant))
m_goup = pd.merge(left=m_gou, right=plant, how='inner', on=['utility_id', 'plant_code', 'year'])
print('CHECK: merged ID-years\t', len(m_goup))

# %%
m = m_goup

# %% 
# # UPDATE COLUMNS: OWNERSHIP TYPE
# m_cln = m.copy()
# update_mask = m_cln.ownership_id.isna().values
# update_to_cols = ['owner_name', 'city_owner', 'state_owner', 'zip_owner', 'ownership_id']
# update_from_cols = ['utility_name', 'city_util', 'state_util', 'zip_util']
# for col_f, col_t in zip(update_to_cols, update_from_cols + ['utility_id']):
#     m_cln.loc[update_mask, col_f] = m_cln.loc[update_mask, col_t]
# m_cln.loc[update_mask, 'ownership'] = 'S' # Single ownership by respondent
# m_cln.loc[~update_mask & (m_cln.percent_owned == 1.), 'ownership'] = 'W' # Wholly owned by an entity other than respondent
# m_cln.loc[~update_mask & (m_cln.percent_owned != 1.), 'ownership'] = 'J' # Jointly owned with another entity
# m_cln.loc[update_mask, 'percent_owned'] = 1.

# %% 
# UPDATE COLUMNS: OPERATION STATUS
# FINDING: MISSING OR UNKOWN STATUS ONLY IN 2006
print('unknown statuses:')
display(m_cln.loc[m_cln.status == 'nan'].groupby(['sheet', 'gen_category'])['utility_id'].count())
# DECISION: These all look like 2006 proposed generators. updating status to `P`
m_cln.loc[m_cln.status == 'nan', 'status'] = 'IP'
mask_op = (m_cln.status == 'OP')
mask_pl = (m_cln.status.isin(('TS', 'P', 'L', 'T', 'U', 'V')))
mask_sb = (m_cln.status.isin(('SB', 'OA', 'OS', 'BU'))) # BU is from 2006 only, stands for "backup"
mask_ds = (m_cln.status.isin(('RE', 'CN', 'IP', 'OT')))
m_cln.loc[mask_op, 'status_simp'] = 'operating'
m_cln.loc[mask_pl, 'status_simp'] = 'planned'
m_cln.loc[mask_sb, 'status_simp'] = 'standby'
m_cln.loc[mask_ds, 'status_simp'] = 'discontinued'

# %%
# UPDATE COLUMNS: OPERATION DATES
vars_month = ['operating_month', 'current_month', 'retirement_month', 'planned_retirement_month']
vars_year = ['operating_year', 'current_year', 'retirement_year', 'planned_retirement_year']
for vm in vars_month:
    m_cln.loc[(m_cln[vm] > 12) | (m_cln[vm] < 1), vm] = 1
for vm, vy in zip(vars_month, vars_year):
    m_cln.loc[m_cln[vy].notna() & m_cln[vm].isna(), vm] = 1
    m_cln.loc[m_cln[vy] == 0, vy] = pd.NA
    m_cln.loc[m_cln[vy].isna(), vm] = pd.NA

m_cln.loc[mask_op | mask_sb | mask_ds , 'dt_operation_start'] = (
    pd.to_datetime(m_cln.loc[mask_op | mask_sb | mask_ds, 'operating_year'].astype(str) + 
                   m_cln.loc[mask_op | mask_sb | mask_ds, 'operating_month'].astype(str), format='%Y%m', errors='coerce'))
m_cln.loc[mask_op | mask_sb, 'dt_operation_end'] = (
    pd.to_datetime(m_cln.loc[mask_op | mask_sb, 'planned_retirement_year'].astype('Int64').astype(str) + 
                   m_cln.loc[mask_op | mask_sb, 'planned_retirement_month'].astype('Int64').astype(str), format='%Y%m', errors='coerce'))
m_cln.loc[mask_pl, 'dt_operation_start'] = (
    pd.to_datetime(m_cln.loc[mask_pl, 'current_year'].astype('Int64').astype(str) + 
                   m_cln.loc[mask_pl, 'current_month'].astype('Int64').astype(str), format='%Y%m', errors='coerce'))
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
display(summ_status)
display(summ_status_simp)

# %%
# UPDATE COLUMNS: PRIME MOVER AND ENERGY SOURCE
m_cln['prime_mover'] = m_cln.prime_mover.str.strip().str.lower()
print('Generators missing prime mover value:', m_cln.prime_mover.isna().sum())
print('Generators missing energy source value:', m_cln.energy_source_1.isna().sum())
display(m_cln.loc[m_cln.energy_source_1.isna(), ])
# DECISION: NO NEED TO UPDATE THIS. ONLY 7 ROWS MAPPING TO GENERATORS THAT WERE INDEFINITELY DISCONTINUED
# Create more general fuel categories
lbl_coal = ['ANT','BIT','LIG','SGC','SUB','WC','RC']
lbl_ng = ['BFG', 'NG', 'OG']
lbl_pet = ['DFO','JF','KER','PC','PG','RFO','SGP','WO']
lbl_ren = ['SUN','WND','GEO','WAT', 'MWH']
lbl_othfuel = ['AB','MSW','OBS','WDS','OBL','SLW','BLQ','WDL','LFG','OBG', 'TDF', 'OTH']
m_cln.loc[m_cln.energy_source_1.isin(lbl_coal), 'energy_source_simp'] = 'coal'
m_cln.loc[m_cln.energy_source_1.isin(lbl_ng), 'energy_source_simp'] = 'natural_gas'
m_cln.loc[m_cln.energy_source_1.isin(lbl_pet), 'energy_source_simp'] = 'petroleum'
m_cln.loc[m_cln.energy_source_1.isin(lbl_ren), 'energy_source_simp'] = 'renewables'
m_cln.loc[m_cln.energy_source_1.isin(lbl_othfuel), 'energy_source_simp'] = 'other_carbon_fuel'
m_cln.loc[m_cln.energy_source_1 == 'NUC', 'energy_source_simp'] = 'nuclear'


#%%
# REORDER COLUMNS
m_final = m_cln[[
    'year', 'ownership_id', 'utility_id', 'plant_code', 'generator_id',
    'ownership', 'owner_name', 'city_owner', 'state_owner', 'zip_owner', 
    'utility_name', 'state_util', 'zip_util',
    'plant_name', 'state_plant', 'zip_plant', 'latitude', 'longitude', 'percent_owned',
    'entity_type', 'naics_primary', 'sector', 
    'prime_mover', 'energy_source_1', 'energy_source_simp', 'cofire_energy_source_1', 'nameplate_capacity_mw', 
    'status', 'status_simp', 'dt_operation_start', 'dt_operation_end',
    # 'generation_kwh', 'total_fuel_consumption_mmbtu', 'fuel_consumption_for_electric_generation_mmbtu',
    # 'tons_of_co2_emissions'
    ]]

# %%
# CONFIRM NO OTHER DUPLICATES
# drop ownership columns
m_final_do = m_final[[col for col in m_final.columns if 'own' not in col]].drop_duplicates()
dup_labels = (m_final_do.groupby(['plant_code', 'generator_id', 'year']).transform('count') > 1)
dup_labels = dup_labels.sum(axis=1) > 0
assert(0 == dup_labels.sum())

# %%
# WRITE TO FILE
m_final.to_parquet(PATH_PROCESSED + 'eia_final.parquet')
