# %%
import pandas as pd
import matplotlib.pyplot as plt

# global variables
PATH_DATA = '../data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 150)

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
dgge = pd.merge(left=dgge, right=pdf[['year', 'plant_code', 'state_plant']], how='left', 
                on=['year', 'plant_code'])
assert (0 == dgge[['year', 'utility_id', 'plant_code', 'generator_id']].duplicated().sum())

# %%
# CREATE FLAGS
start_year = 2015
dgge['is_duke'] = dgge.utility_id.isin(duke_utilities.utility_id.values)
dgge['is_duke_startyear'] = (dgge.year == start_year) & dgge.utility_id.isin(duke_utilities.utility_id.values)
dgge['is_duke_startyear'] = dgge.groupby(['plant_code', 'generator_id'])['is_duke_startyear'].transform('max')
dgge['capacity_factor'] = dgge.net_gen_tot_an / (dgge.nameplate_capacity_mw*8760)



# %%
groupers = [['state_plant', 'year'], ['energy_source_1', 'year'], ['year']]
summ_df = pd.DataFrame()
for group in groupers:
    df = (dgge.groupby(group)
        .agg({'utility_id':'nunique', 'plant_code':'nunique', 'gid':'nunique',
            'nameplate_capacity_mw':'sum', 'co2_mass_short_tons_gen':'sum',
            'co2_mass_short_tons_gen_923':'sum',
            'net_gen_tot_an':'sum',
            'capacity_factor':['mean', 'std']})
            .reset_index())
    summ_df = pd.concat([df, summ_df], ignore_index=True)

# %%
vars_id = [('year',''), ('energy_source_1',''), ('state_plant','')]
vars_other = [col for col in summ_df.columns if col not in vars_id]
summ_df[vars_id] = summ_df[vars_id].fillna('TOTAL')
summ_df[vars_id + vars_other]
# %%
