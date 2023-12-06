# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
from tqdm import tqdm

# global variables
PATH_DATA = '../data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
PATH_RESULTS = '../results/analysis/duke/'
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)

# read in data
gdf = pd.read_parquet(PATH_PROCESSED + 'df_generators.parquet')
gendf = pd.read_parquet(PATH_PROCESSED + 'df_generation.parquet')
emdf = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'df_plants.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'df_utilities.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'df_owners.parquet')

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
dgg = pd.merge(left=dg, right=gendf, 
               how='left', on=['year', 'utility_id', 'plant_code', 'generator_id'])
assert (0 == dgg[['year', 'utility_id', 'plant_code', 'generator_id']].duplicated().sum())
dgge = pd.merge(left=dgg, right=emdf.drop(columns=['gid', 'nameplate_capacity_mw']), how='left', on=['year', 'plant_code', 'generator_id'])
assert (0 == dgge[['year', 'utility_id', 'plant_code', 'generator_id']].duplicated().sum())
dgge = pd.merge(left=dgge, right=pdf[['year', 'plant_code', 'state_plant', 'latitude', 'longitude']], how='left', 
                on=['year', 'plant_code'])
assert (0 == dgge[['year', 'utility_id', 'plant_code', 'generator_id']].duplicated().sum())


# %%
# CREATE FLAGS
start_year = 2013
dgge['is_duke'] = dgge.utility_id.isin(duke_utilities.utility_id.values)
dgge['is_duke_startyear'] = (dgge.year == start_year) & dgge.utility_id.isin(duke_utilities.utility_id.values)
dgge['is_duke_startyear'] = dgge.groupby(['plant_code', 'generator_id'])['is_duke_startyear'].transform('max')
dgge['capacity_factor'] = dgge.net_gen_tot_an / (dgge.nameplate_capacity_mw*8760)
# subset to generators that are either operating or have emissoins/generation data
dgge['is_operating_status'] = dgge.status.isin(['OP', 'OA', 'OS', 'SB'])
dgge['is_operating'] = (
    (dgge.co2_mass_tons_gen_923 > 0)
    | (dgge.co2_mass_tons_gen > 0)
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
            'nameplate_capacity_mw':'sum', 'co2_mass_tons_gen':'sum',
            'co2_mass_tons_gen_923':'sum',
            'net_gen_tot_an':'sum',
            'capacity_factor':['mean', 'std']})
            .reset_index())
    summ_df = pd.concat([df, summ_df], ignore_index=True)

vars_id = [('year',''), ('energy_source_1',''), ('state_plant','')]
vars_other = [col for col in summ_df.columns if col not in vars_id]
summ_df[vars_id] = summ_df[vars_id].fillna('TOTAL')
summ_df[vars_id + vars_other].to_csv(PATH_RESULTS + 'summary.csv', index=False)

# %%
# GENERATION MIX BY STATE
import seaborn as sns
var = 'net_gen_tot_an'
show_pct = True
options = {
    'nameplate_capacity_mw':{'units':'MW', 'denom':1},
    'net_gen_tot_an':{'units':'TWh', 'denom':1e6}
}
gm = (dgge_op
      .groupby(['year', 'state_plant', 'energy_source_1_subcat'])
      .agg({var:'sum'}) / options[var]['denom']).reset_index()
gm[f'{var}_state'] = gm.groupby(['year', 'state_plant'])[var].transform('sum')
gm[f'{var}_pct'] = gm[var] / gm[f'{var}_state']

# Pivot the data to get energy type breakdown for each state by year
var_summ = var if not show_pct else f'{var}_pct'
var_summ_unit = 'pct' if show_pct else options[var]['units']
gmt = gm.pivot_table(
    values=var_summ, index=['state_plant', 'year'], 
    columns='energy_source_1_subcat', aggfunc='sum')
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
    state_data.plot(kind='bar', x='year', legend=False, stacked=True, cmap='tab20', ax=ax, title=f'{var} in {state}')
    ax.set_xlabel('Year')
    ax.set_ylabel(f'{var_summ_unit}')
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
# STATUS OVER TIME
# plot
# Initialize a figure
def plot_egu_status(df, status_mapping, ax, ylab):
    pivot_df = df.pivot_table(index='gid', columns='year', values='status_code', fill_value=None).astype('Int64')
    pivot_df.sort_values(pivot_df.columns.to_list(), ascending=True, inplace=True)
    # Assign a color to each status
    colors = plt.cm.viridis(np.linspace(0, 1, len(status_mapping)))

    # Iterate through each EGU and plot
    # Iterate through each EGU and plot segments
    for i, (egu_id, data) in enumerate(pivot_df.iterrows()):
        # Previous year and status to start the first segment
        prev_year = data.index[0]
        prev_status = data.iloc[0]
        
        # Iterate through each year for the EGU
        for year, status in data.items():
            # Skip if data is missing
            if pd.isna(prev_status):
                prev_year = year
                prev_status = status
                continue
            # Plot the segment from the previous year to the current year
            ax.plot([prev_year, year], [i, i], lw=1, c=colors[prev_status])

            # Update previous year and status
            prev_year = year
            prev_status = status
        if pd.notna(prev_status):
            ax.plot([prev_year, prev_year+1], [i, i], lw=1, c=colors[prev_status])
    # Customizing the plot
    # Create a colorbar
    norm = mcolors.Normalize(vmin=min(status_mapping.values()), vmax=max(status_mapping.values()))
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
    sm.set_array([])
    # Add the colorbar to the figure
    cbar = plt.colorbar(sm, ax=ax, ticks=range(len(status_mapping)), label='Status')
    cbar.set_ticklabels(list(status_mapping.keys()))  
    ax.set_ylabel(f'{ylab} EGUs')
    return ax

# Create a mapping for the statuses

status_order = ['OP', 'SB', 'OA', 'P', 'L', 'T', 'U', 'V', 'IP', 'CN', 'RE']
status_mapping = {status:i for i, status in enumerate(status_order)}
dgge['status_code'] = dgge.status.map(status_mapping)

FIGSIZE=3
states = dgge.state_plant.unique()
states = ['NC', 'SC', 'FL', 'IN', 'OH', 'KY']
for state in tqdm(states):
    dgge_state = dgge.loc[dgge.state_plant == state]
    fuels = dgge_state.energy_source_1_subcat.unique()
    fig, ax = plt.subplots(nrows=len(fuels), sharex=True, figsize=(3*FIGSIZE, len(fuels)*FIGSIZE))
    for i, fuel in enumerate(fuels):
        plot_egu_status(dgge_state.loc[(dgge_state.energy_source_1_subcat == fuel)], 
                        status_mapping, ax[i], fuel)
    ax[0].set_title(f'EGU status in {state}, by fuel type')
    plt.savefig(PATH_RESULTS + f'fig_egu_status_in_{state}.png', dpi=300, bbox_inches='tight')


# %%
