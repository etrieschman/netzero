# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

# global variables
PATH_DATA = '../data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
PATH_RESULTS = '../results/'
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)

# read in data
gdf = pd.read_parquet(PATH_PROCESSED + 'eia860_generator.parquet')
gendf = pd.read_parquet(PATH_PROCESSED + 'df_generation.parquet')
emdf = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'eia860_plant.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'eia860_utility.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'eia860_ownership.parquet')

# %%
# GET UTILITIES, PLANTS, AND GENERATORS
duke_utilities = (udf
                  .loc[udf.utility_name.str.lower().str.contains('duke|degs'), ['utility_id', 'utility_name']]
                  .drop_duplicates().reset_index(drop=True))
duke_plants = (pdf
               .loc[pdf.utility_id.isin(duke_utilities.utility_id.values), ['utility_id', 'plant_code']]
               .drop_duplicates().reset_index(drop=True))
duke_generators = (gdf
                   .loc[gdf.plant_code.isin(duke_plants.plant_code.values), ['utility_id', 'plant_code', 'generator_id']]
                   .drop_duplicates().reset_index(drop=True))
duke_generators['gid'] = duke_generators.plant_code.astype(str) + '_' + duke_generators.generator_id


# %%
# MAKE ANALYSIS DATASET
dg = gdf.copy()
dg['gid'] = dg.plant_code.astype(str) + '_' + dg.generator_id
dg = dg.loc[dg.gid.isin(duke_generators.gid.values)]
dgg = pd.merge(left=dg, right=gendf.drop(columns=['energy_source_1', 'status']), 
               how='left', on=['year', 'utility_id', 'plant_code', 'generator_id'])
assert (0 == dgg[['year', 'utility_id', 'plant_code', 'generator_id']].duplicated().sum())
dgge = pd.merge(left=dgg, right=emdf.drop(columns='gid'), how='left', on=['year', 'plant_code', 'generator_id'])
assert (0 == dgge[['year', 'utility_id', 'plant_code', 'generator_id']].duplicated().sum())
dgge = pd.merge(left=dgge, right=pdf[['year', 'plant_code', 'state_plant', 'latitude', 'longitude']], how='left', 
                on=['year', 'plant_code'])
assert (0 == dgge[['year', 'utility_id', 'plant_code', 'generator_id']].duplicated().sum())


# %%
# CREATE FLAGS
start_year = 2015
dgge['is_duke'] = dgge.utility_id.isin(duke_utilities.utility_id.values)
dgge['is_duke_startyear'] = (dgge.year == start_year) & dgge.utility_id.isin(duke_utilities.utility_id.values)
dgge['is_duke_startyear'] = dgge.groupby(['plant_code', 'generator_id'])['is_duke_startyear'].transform('max')
dgge['capacity_factor'] = dgge.net_gen_tot_an / (dgge.nameplate_capacity_mw*8760)
# subset to generators that are either operating or have emissoins/generation data
dgge['is_operating_status'] = dgge.status.isin(['OP', 'OA', 'OS', 'SB'])
dgge['is_operating'] = (
    (dgge.co2_mass_short_tons_gen_923 > 0)
    | (dgge.co2_mass_short_tons_gen > 0)
    | (dgge.net_gen_tot_an > 0)
)
dgge_op = dgge.loc[dgge.is_operating | dgge.is_operating_status]


# %%
# TOPLINE OPERATIONS TABLE
groupers = [['state_plant', 'year'], ['energy_source_1', 'year'], ['year']]
summ_df = pd.DataFrame()
for group in groupers:
    df = (dgge_op.groupby(group)
        .agg({'utility_id':'nunique', 'plant_code':'nunique', 'gid':'nunique',
            'nameplate_capacity_mw':'sum', 'co2_mass_short_tons_gen':'sum',
            'co2_mass_short_tons_gen_923':'sum',
            'net_gen_tot_an':'sum',
            'capacity_factor':['mean', 'std']})
            .reset_index())
    summ_df = pd.concat([df, summ_df], ignore_index=True)

vars_id = [('year',''), ('energy_source_1',''), ('state_plant','')]
vars_other = [col for col in summ_df.columns if col not in vars_id]
summ_df[vars_id] = summ_df[vars_id].fillna('TOTAL')
summ_df[vars_id + vars_other].to_csv(PATH_RESULTS + 'duke/summary.csv', index=False)

# %%
# GENERATION MIX BY STATE
import seaborn as sns
gm = (dgge_op
      .groupby(['year', 'state_plant', 'energy_source_1'])
      .agg({'net_gen_tot_an':'sum'}) / 1e6).reset_index()
gm['net_gen_tot_an_state'] = gm.groupby(['year', 'state_plant'])['net_gen_tot_an'].transform('sum')
gm['net_gen_pct'] = gm.net_gen_tot_an / gm.net_gen_tot_an_state

# Pivot the data to get energy type breakdown for each state by year
gmt = gm.pivot_table(
    values='net_gen_pct', index=['state_plant', 'year'], 
    columns='energy_source_1', aggfunc='sum')
# Fill missing values with 0 (if any)
gmt = gmt.fillna(0).reset_index()
# Now we plot a stacked bar chart for each state
states = gm['state_plant'].unique()
states = ['NC', 'SC', 'FL', 'IN', 'OH', 'KY']

# Set up the matplotlib figure and axes
figscale=5
ncols=2
fig, axes = plt.subplots(
    nrows=len(states)//ncols, ncols=ncols, figsize=(5*ncols, 4*len(states)//ncols), sharex=True, sharey=True)
# Loop through each state and create a stacked bar chart
for ax, state in zip(axes.flatten(), states):
    state_data = gmt.loc[gmt.state_plant == state]
    state_data.plot(kind='bar', x='year', legend=False, stacked=True, cmap='tab20', ax=ax, title=f'Net Generation in {state}')
    ax.set_xlabel('Year')
    ax.set_ylabel('Net Generation (pct)')
axes[2,1].legend(title='Energy Type')

# Adjust the layout to prevent overlap
plt.tight_layout()
plt.show()



# %%
# MAP OF OPERATIONS
year = 2018

# Read in the US states shapefile using GeoPandas
gdf_states = gpd.read_file(PATH_DATA + 'resources/cb_state_boundaries/cb_2018_us_state_20m.shp')
gdf_states = gdf_states.loc[~gdf_states.STUSPS.isin(['HI', 'AK', 'PR'])]

# MAP GENERATION
# Convert the DataFrame to a GeoDataFrame
dgge_map = (dgge_op
            .loc[(dgge_op.year == year)]
            .groupby(['latitude', 'longitude'])
            [['net_gen_tot_an']].sum()).reset_index()
gdf = gpd.GeoDataFrame(dgge_map, geometry=gpd.points_from_xy(dgge_map.longitude, dgge_map.latitude))
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_states.geometry.plot(ax=ax, linewidth=0.8, color='lightgrey', alpha=0.5, edgecolor='black')
gdf.plot(ax=ax, 
         column='net_gen_tot_an', cmap='YlOrRd', markersize=np.sqrt(gdf.net_gen_tot_an)/10,
         legend=True, legend_kwds={'label': "Electricity Generation (MWh)",
                                   'shrink':0.75},
         alpha=0.75)
plt.xlim((-107, -75))
plt.title(f'All Duke Energy plants in {year}, by total annual generation')
plt.show()

# MAP FUEL TYPES
dgge_map = (dgge_op
            .loc[(dgge_op.year == year) & (dgge_op.net_gen_tot_an > 0),
                 ['latitude', 'longitude', 'energy_source_1']]
             .drop_duplicates().reset_index(drop=True))
gdf = gpd.GeoDataFrame(dgge_map, geometry=gpd.points_from_xy(dgge_map.longitude, dgge_map.latitude))
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_states.geometry.plot(ax=ax, linewidth=0.8, color='lightgrey', edgecolor='black')
gdf.plot(ax=ax, 
         column='energy_source_1', s=250, cmap='tab10', 
         legend=True,
         alpha=0.75)
plt.xlim((-107, -75))
plt.title(f'Any Duke Energy plants producing electricity in {year}, by energy source')
plt.show()

# %%
# RETIREMENTS AND NEW PLANTS
mask_gen_or_em = ((dgge.net_gen_tot_an.notna() & (dgge.net_gen_tot_an != 0))
                   | (dgge.co2_mass_short_tons_gen > 0))
dgge['retired_ops'] = (dgge.year >= dgge.dt_operation_end.dt.year)
dgge['retired_status'] = (dgge.status == 'RE')
dgge.loc[(dgge.retired_ops | dgge.retired_status) & mask_gen_or_em]


# %%
# DIVESTMENTS
