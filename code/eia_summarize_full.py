# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import itertools

from utils import PATH_PROCESSED, PATH_RESULTS
from netzero.code.transform.utils_summ import summarize_id_counts_byyear

# %%
m_final = pd.read_parquet(PATH_PROCESSED + 'eia_final.parquet')

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
