# %%
# SET UP
import pandas as pd
import networkx as nx

# global variables
PATH_DATA = '../../data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 150)

# %%
# READIN DATA
gdf = pd.read_parquet(PATH_PROCESSED + 'eia860_generator.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'eia860_plant.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'eia860_utility.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'eia860_ownership.parquet')
edf = pd.read_parquet(PATH_PROCESSED + 'epa_emissions.parquet')
xw = pd.read_csv(PATH_DATA + 'epa_eia_crosswalk.csv')

# %% 
# CLEAN UP CROSSWALK
def clean_xw(xw):
    # drop and rename columns
    xw.columns = xw.columns.str.lower()
    xw.drop(columns=[col for col in xw.columns if col.startswith('mod_')], inplace=True)
    xw.drop(columns=['sequence_number', 'camd_state', 'camd_facility_name',
                    'camd_latitude', 'camd_longitude', 'camd_nameplate_capacity',
                    'eia_state', 'eia_plant_name', 'eia_latitude', 'eia_longitude',
                    'eia_nameplate_capacity', 'plant_id_change_flag',
                    'match_type_boiler'], inplace=True)
    # check: make sure we're not dropping a unique identifying column
    assert len(xw) == len(xw.drop_duplicates())

    # drop CAMD excluded and unmatched
    print('Dropping records labeled "excluded or unmatched"...')
    print('Rows dropped:\t', len(xw.loc[xw.match_type_gen.isin(['Manual CAMD Excluded', 'CAMD Unmatched'])]))
    xw = xw.loc[~xw.match_type_gen.isin(['Manual CAMD Excluded', 'CAMD Unmatched'])]
    # drop boiler id and collapse on remaining unique ids
    print('Dropping EIA boiler IDs and recollapsing...')
    print('Rows dropped:\t', xw.drop(columns='eia_boiler_id').duplicated().sum())
    xw = xw.drop(columns='eia_boiler_id').drop_duplicates()
    # drop nameplate capacity information and recollapse
    print('Dropping CAMD generator IDs and recollapsing...')
    print('Rows dropped:\t', xw.drop(columns=['camd_generator_id']).duplicated().sum())
    xw = xw.drop(columns='camd_generator_id').drop_duplicates()

    # update datatypes
    xw['eia_plant_id'] = pd.to_numeric(xw.eia_plant_id).astype(int)

    # make ids
    xw['camd_uid'] = 'camd_' + xw.camd_plant_id.astype(str) + '_' + xw.camd_unit_id
    xw['eia_gid'] = 'eia_' + xw.eia_plant_id.astype(str) + '_' + xw.eia_generator_id

    # check: rows unique on the identifiers we want
    assert len(xw) == len(xw[['camd_uid', 'eia_gid']].drop_duplicates())
    print('Returning clean and deduped crosswalk...')
    return xw

xwc = clean_xw(xw.copy())

# %%
# add merge summary
xwcm = pd.merge(left=xwc, right=pd.Series(xwc.camd_uid.value_counts(), name='merge_camd'),
               how='outer', left_on='camd_uid', right_index=True).reset_index(drop=True)


xwcm = pd.merge(left=xwcm, right=pd.Series(xwc.eia_gid.value_counts(), name='merge_eia'),
               how='outer', left_on='eia_gid', right_index=True).reset_index(drop=True)

xwcm[['merge_camd', 'merge_eia']].value_counts()


# %%
# Create subplant associations
# make edgegraph and assert that it's bipartite
graph = nx.from_pandas_edgelist(
    xwcm, source='camd_uid', target='eia_gid', edge_attr=True)
assert nx.algorithms.bipartite.is_bipartite(graph)

# label subgraphs
for i, node_set in enumerate(nx.connected_components(graph)):
    subgraph = graph.subgraph(node_set)
    nx.set_edge_attributes(subgraph, name='subplant_id', values=f'sp_{i}')

nx.to_pandas_edgelist(graph).rename({'source':'camd_uid', 'target':'eia_gid'})

# %%
# Allocation method 1
# Assign unit emissions proportional to EGU nameplate capacity
# B. Merge Crosswalk onto emissions data
#  - How much emissions are not in crosswalk?
# C. Merge gdf onto emisisons data
# - what percent of generation is accounted for by year?
# D. Calculate percent to allocate to generators
# - allocate


# %%
# Allocation method 2
# Assign subplant emissions proportional to EGU nameplate capacity
# A. Merge subplant crosswalk onto emissions data
#  - how much is unaccounted for?
# B. Aggregate to subplant emissions
# C. Merge on generator data
# D. Allocate emissions to generation


# %%
# CALCULATE EMISSIONS FROM EIA DATA
# fuel-specific emisisons factors
# https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/emission_factors_for_co2_ch4_n2o.csv