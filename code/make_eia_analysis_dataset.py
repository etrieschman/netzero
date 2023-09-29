# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import itertools
import os, re

from utils import PATH_PROCESSED, PATH_RESULTS, summarize_id_counts_byyear

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 120)

# %%
# READIN DATA
generator = pd.read_parquet(PATH_PROCESSED + 'eia_f860_generator.parquet')
plant = pd.read_parquet(PATH_PROCESSED + 'eia_f860_plant.parquet')
utility = pd.read_parquet(PATH_PROCESSED + 'eia_f860_utility.parquet')
owner = pd.read_parquet(PATH_PROCESSED + 'eia_f860_ownership.parquet')

# %%
# Q: DOES OWNERSHIP ADD UP? No
# DECISIONS:
# 0. Divide by 100 in settings where decimal is off
# 1. If owner_count is 1 and percent ownership is "close" to 1, nudge to 1
# 2. Set all missing ownership to 100% if only 1 owner; 
#    Set missing ownership to 0% if missing, but sum of ownership is 1
# 3. TODO: 2010 is a messy year that needs to be reconciled
owner['pct_ownership_total'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('sum').round(5)
owner['owner_count'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('count')
owner.loc[owner.pct_ownership_total == 100., 'percent_owned'] /= 100.
owner.loc[(owner.owner_count == 1) & 
    (np.abs(owner.pct_ownership_total - 1) <= 0.02), 'percent_owned'] = 1.
owner.loc[(owner.percent_owned.isna()) & (owner.owner_count == 0), 'percent_owned'] = 1.
owner.loc[(owner.percent_owned.isna()) & (owner.owner_count > 0), 'percent_owned'] = 0.
# print(owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'])['percent_owned'].sum().value_counts().sort_index())
owner['pct_ownership_total'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('sum').round(5)
owner['owner_count'] = owner.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)['percent_owned'].transform('count')
print("Count of generators where ownership doesn't add up after cleaning, by year:")
display(owner.loc[np.abs(owner.pct_ownership_total - 1.) >= 0.02].groupby('year')['generator_id'].agg(['count', 'nunique']))
owner = owner.drop(columns=['pct_ownership_total', 'owner_count'])

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
# UPDATE COLUMNS: OWNERSHIP TYPE
m_cln = m.copy()
update_mask = m_cln.ownership_id.isna().values
update_to_cols = ['owner_name', 'city_owner', 'state_owner', 'zip_owner', 'ownership_id']
update_from_cols = ['utility_name', 'city_util', 'state_util', 'zip_util']
for col_f, col_t in zip(update_to_cols, update_from_cols + ['utility_id']):
    m_cln.loc[update_mask, col_f] = m_cln.loc[update_mask, col_t]
m_cln.loc[update_mask, 'ownership'] = 'S' # Single ownership by respondent
m_cln.loc[~update_mask & (m_cln.percent_owned == 1.), 'ownership'] = 'W' # Wholly owned by an entity other than respondent
m_cln.loc[~update_mask & (m_cln.percent_owned != 1.), 'ownership'] = 'J' # Jointly owned with another entity
m_cln.loc[update_mask, 'percent_owned'] = 1.

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
display(m_cln.loc[m_cln.energy_source_1.isna(), ].groupby('status')['utility_id'].count())
# DECISION: NO NEED TO UPDATE THIS. ONLY 6 ROWS MAPPING TO GENERATORS THAT WERE INDEFINITELY DISCONTINUED
# Create more general categories
lbl_coal = ['ANT','BIT','LIG','SGC','SUB','WC','RC']
lbl_ng = ['BFG', 'NG', 'OG']
lbl_pet = ['DFO','JF','KER','PC','PG','RFO','SGP','WO']
lbl_ren = ['SUN','WND','GEO','WAT']
lbl_othfuel = ['AB','MSW','OBS','WDS','OBL','SLW','BLQ','WDL','LFG','OBG']
m_cln.loc[m_cln.energy_source_1.isin(lbl_coal), 'energy_source_simp'] = 'coal'
m_cln.loc[m_cln.energy_source_1.isin(lbl_ng), 'energy_source_simp'] = 'natgas'
m_cln.loc[m_cln.energy_source_1.isin(lbl_pet), 'energy_source_simp'] = 'petrl'
m_cln.loc[m_cln.energy_source_1.isin(lbl_ren), 'energy_source_simp'] = 'renew'
m_cln.loc[m_cln.energy_source_1.isin(lbl_othfuel), 'energy_source_simp'] = 'othfuel'
m_cln.loc[m_cln.energy_source_1 == 'NUC', 'energy_source_simp'] = 'nuclear'
m_cln[['energy_source_1', 'energy_source_simp']].value_counts(dropna=False)

#%%
# REORDER COLUMNS
m_final = m_cln[[
    'year', 'ownership_id', 'utility_id', 'plant_code', 'generator_id',
    'ownership', 'owner_name', 'city_owner', 'state_owner', 'zip_owner', 
    'utility_name', 'state_util', 'zip_util',
    'plant_name', 'state_plant', 'zip_plant', 'latitude', 'longitude', 'percent_owned',
    'entity_type', 'naics_primary', 'sector', 
    'prime_mover', 'energy_source_1', 'cofire_energy_source_1', 'nameplate_capacity_mw', 
    'status', 'status_simp', 'dt_operation_start', 'dt_operation_end',
    # 'generation_kwh', 'total_fuel_consumption_mmbtu', 'fuel_consumption_for_electric_generation_mmbtu',
    # 'tons_of_co2_emissions'
    ]]

m_final.to_parquet(PATH_PROCESSED + 'eia_final.parquet')

# %%
# SUMMARIZE IDS BY YEAR
m_final['pid'] = m_final.utility_id.astype(str) + '_' + m_final.plant_code.astype(str)
m_final['oid'] = m_final.ownership_id
m_final['gid'] = m_final.pid + '_' + m_final.generator_id
summarize_id_counts_byyear(m_final, ['oid', 'pid', 'gid'])

# %%
#  SUMMARIZE GENERATORS BY SECTOR
m_final.loc[:, 'chp'] = np.where(m_final.sector.isin([3, 5, 7]), True, False)
m_final.loc[m_final.sector == 1, 'sector_ii'] = 'electric_utility'
m_final.loc[m_final.sector.isin([2,3]), 'sector_ii'] = 'ipp'
m_final.loc[m_final.sector.isin([4,5]), 'sector_ii'] = 'commercial'
m_final.loc[m_final.sector.isin([6,7]), 'sector_ii'] = 'industrial'
(m_final.loc[m_final.year.isin([2013, 2021]) & (m_final.status_simp == 'operating'), ['year', 'sector_ii', 'chp', 'gid', 'nameplate_capacity_mw']]
 .drop_duplicates()
 .groupby(['year', 'sector_ii', 'chp'])[['gid', 'nameplate_capacity_mw']]
 .agg({'gid':'count', 'nameplate_capacity_mw':'sum'}))

# %%
# VISUALIZE FACILITY DATA
import geopandas as gpd
from shapely.geometry import Point
m_final['latitude'] = pd.to_numeric(m_final.latitude, errors='coerce')
m_final['longitude'] = pd.to_numeric(m_final.longitude, errors='coerce')
m_final_op = (m_final
        .loc[(m_final.year == 2021) & (m_final.status_simp == 'operating') ,
             ['plant_code', 'prime_mover', 'generator_id', 'latitude', 'longitude', 'nameplate_capacity_mw']]
             .drop_duplicates()
        .groupby(['latitude', 'longitude', 'prime_mover'])
        .agg({'nameplate_capacity_mw':'count'})
        .reset_index())
geometry = [Point(xy) for xy in zip(m_final_op['longitude'], m_final_op['latitude'])]
m_final_op = gpd.GeoDataFrame(m_final_op, geometry=geometry)
m_final_op.plot(column='prime_mover', cmap='Paired', alpha=0.25, markersize=3, legend=False)
plt.show()

# %%
# MAKE DATASET FOR SANKEY DIAGRAM
# make year-gen dataset
m_util = m_final.drop(columns=['ownership_id', 'ownership', 'owner_name', 'city_owner', 'state_owner', 'zip_owner', 'percent_owned']).drop_duplicates()
m_util['id'] = m_util.plant_code.astype(str) + '_' + m_util.generator_id
base = pd.DataFrame(itertools.product(m_util.id.drop_duplicates().values, m_util.year.drop_duplicates().values), columns=['id', 'year'])
mu_all = pd.merge(left=base, right=m_util[['id', 'year', 'status', 'status_simp']], on=['id', 'year'], how='left').drop_duplicates().sort_values(['id', 'year'])
mu_all['status_prev'] = mu_all.groupby('id')['status_simp'].shift(1)
mu_all['status'] = mu_all.status_simp
mu_all_summ = mu_all.groupby(['year', 'status_prev', 'status'], dropna=False)['id'].count().reset_index()
mu_all_summ.loc[mu_all_summ.status_prev.isna(), 'status_prev'] = 'nan'
mu_all_summ.loc[mu_all_summ.status.isna(), 'status'] = 'nan'
mu_all_summ['status_prev'] = (mu_all_summ.year - 1).astype(str) + '-' + mu_all_summ.status_prev
mu_all_summ['status'] = mu_all_summ.year.astype(str) + '-' + mu_all_summ.status
# make status index mapping
unique_statusyears = np.sort(pd.concat([mu_all_summ.status_prev, mu_all_summ.status]).unique())
value_to_index = {value:index for index, value in enumerate(unique_statusyears)}
mu_all_summ['status_prev_idx'] = mu_all_summ.status_prev.map(value_to_index)
mu_all_summ['status_idx'] = mu_all_summ.status.map(value_to_index)

# %%
# PLOT SANKEY FIGURE
node_colors = np.empty(len(unique_statusyears), dtype=object)
for i, k in enumerate(unique_statusyears):
    if 'nan' in k: node_colors[i] = 'grey'
    elif 'discontinued' in k: node_colors[i] = 'red'
    elif 'standby' in k: node_colors[i] = 'yellow'
    elif 'planned' in k: node_colors[i] = 'blue'
    elif 'operating' in k: node_colors[i] = 'green'
link_colors = np.empty(len(mu_all_summ), dtype=object)
for i, k in enumerate(mu_all_summ.status_prev.values):
    if 'nan' in k: link_colors[i] = 'grey'
    elif 'discontinued' in k: link_colors[i] = 'red'
    elif 'standby' in k: link_colors[i] = 'yellow'
    elif 'planned' in k: link_colors[i] = 'blue'
    elif 'operating' in k: link_colors[i] = 'green'

fig = go.Figure(data=[go.Sankey(
    node = dict(
    #   pad = 15,
    #   thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = unique_statusyears,
      color = node_colors
    ),
    link = dict(
      source = mu_all_summ.status_prev_idx.values, # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = mu_all_summ.status_idx.values,
      value = mu_all_summ.id.values,
      color = link_colors
  ))])

fig.update_layout(title_text='Generator status over time', font_size=8)
fig.show()
fig.write_html(PATH_RESULTS + 'plt_sankey_status.html')

# %%
# MAKE DATASET FOR SANKEY DIAGRAM
# make year-gen dataset
m_own = m_final.copy()
m_own['id'] = m_own.plant_code.astype(str) + '_' + m_own.generator_id
m_own = m_own.groupby(['id', 'year', 'ownership'])['ownership_id'].apply(set).reset_index()
base = pd.DataFrame(itertools.product(m_own.id.drop_duplicates().values, m_own.year.drop_duplicates().values), columns=['id', 'year'])

s = pd.merge(left=base, right=m_own, on=['id', 'year'], how='left').sort_values(['id', 'year'])
s['ownership_id_prev'] = s.groupby('id')['ownership_id'].shift(1)
s['ownership_prev'] = s.groupby('id')['ownership'].shift(1)
s.loc[s.ownership_prev.isna(), 'ownership_prev'] = 'nan'
s.loc[s.ownership.isna(), 'ownership'] = 'nan'
s['ownership_prev'] = (s.year - 1).astype(str) + '-' + s.ownership_prev
s['ownership'] = s.year.astype(str) + '-' + s.ownership

# %%
# TABLE SUMMARY
s.loc[s.ownership_id_prev.isna(), 'status'] = 'new generator'
s.loc[s.ownership_id_prev.notna() & (s.ownership_id == s.ownership_id_prev), 'status'] = 'same owner'
s.loc[s.ownership_id_prev.notna() & (s.ownership_id != s.ownership_id_prev), 'status'] = 'change in ownership'
s.groupby(['year', 'status'])['id'].count().reset_index().pivot(index='year', columns='status', values='id').plot()
plt.title('Count of generators, by how ownership changed in current year ')

# %% 
# SANKEY DATASET
ss = s.groupby(['year', 'ownership_prev', 'ownership'], dropna=False)['id'].count().reset_index()
# make ownership index mapping
unique = np.sort(pd.concat([ss.ownership_prev, ss.ownership]).unique().astype(str))
value_to_index = {value:index for index, value in enumerate(unique)}
ss['ownership_prev_idx'] = ss.ownership_prev.map(value_to_index)
ss['ownership_idx'] = ss.ownership.map(value_to_index)

# %%
# PLOT SANKEY FIGURE
node_colors = np.empty(len(unique), dtype=object)
for i, k in enumerate(unique):
    if 'nan' in k: node_colors[i] = 'grey'
    elif 'S' in k: node_colors[i] = 'blue'
    elif 'J' in k: node_colors[i] = 'purple'
    elif 'W' in k: node_colors[i] = 'red'
link_colors = np.empty(len(ss), dtype=object)
for i, k in enumerate(ss.ownership_prev.values):
    if 'nan' in k: link_colors[i] = 'grey'
    elif 'S' in k: link_colors[i] = 'blue'
    elif 'J' in k: link_colors[i] = 'purple'
    elif 'W' in k: link_colors[i] = 'red'
fig = go.Figure(data=[go.Sankey(
    node = dict(
    #   pad = 15,
    #   thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = unique,
      color = node_colors
    ),
    link = dict(
      source = ss.ownership_prev_idx.values, # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = ss.ownership_idx.values,
      value = ss.id.values,
      color = link_colors
  ))])

fig.update_layout(title_text='Generator ownership status over time', font_size=10)
fig.show()
fig.write_html(PATH_RESULTS + 'plt_sankey_ownership.html')
# %%
