# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

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
# PREP DATASET FOR SUMMARY
# 0. drop ownership
print('N records with ownership:', len(eia))
eia_do = eia[[col for col in eia.columns if 'own' not in col]].drop(columns='cofire_energy_source_1').copy()
eia_do = eia_do.drop_duplicates()
print('N records w/o ownership:', len(eia_do))
eia_do.loc[(eia_do.groupby(['key_gen', 'year']).transform('count') > 1).any(axis=1)]
# TODO: Fix disparate energy sources and disparate operation start/end dates

# %%
# 1. TODO: Create an epa merge flag in eia (placeholders used here)
eia_do['gen_in_epa'] = np.random.choice([True, False], len(eia_do))
eia_do['plant_in_epa'] = eia_do.groupby(['plant_code', 'year'])['gen_in_epa'].transform(lambda x: x.any())
eia_do['fuel_simp'] = np.random.choice([1,2,3,4], len(eia_do))
eia_do['not_renewable'] = eia_do.fuel_simp != 4
eia_do['plant_has_any_ff'] = eia_do.groupby(['plant_code', 'year'])['not_renewable'].transform('sum') > 0
eia_do['gen_age_yrs'] = (dt.datetime.today() - eia_do.dt_operation_start).dt.days / 365.25

# %%
# 2. TODO: drop relevant statuses (e.g., discontinued)
display(eia_do.groupby(['status', 'gen_in_epa'])['key_gen'].nunique())

# %%
# 3. Summarize at the plant level
# helper summary function
def get_inepa_plantsumm(df, group):
    summ_plant = (df.groupby(['utility_id', 'plant_code', 'state_plant', 'entity_type', 'year', 'plant_in_epa'])
                  .agg({'nameplate_capacity_mw':'sum', 'key_gen':'nunique', 'not_renewable':'mean', 'gen_age_yrs':'mean'})
                  .reset_index())
    summ = (summ_plant.groupby([group, 'plant_in_epa'])
            .agg({'plant_code':['nunique'], 'key_gen':['sum', 'mean', 'std'], 'nameplate_capacity_mw':['sum', 'mean', 'std'],
                   'not_renewable':['mean', 'std'], 'gen_age_yrs':['mean', 'std']})
            .reset_index())
    summ.rename(columns={group:'value', 'not_renewable':'pct_gens_not_renewable', 'gen_age_yrs':'mean_gen_age_yrs', 
                         'plant_code':'plants', 'key_gen':'generators'}, inplace=True)
    return summ

# drop plants that have no fossil fuels
# NOTE: THIS WOULD INCLUDE NON-FOSSIL FUEL CAPACITY AT PLANTS WITH ANY FOSSIL FUELS
print(eia_do.groupby('plant_has_any_ff')['plant_code'].nunique())
eia_do_dp = eia_do.loc[eia_do.plant_has_any_ff]
# summarize
plant_summ = pd.DataFrame([])
groups = ['year', 'entity_type', 'state_plant']
for g in groups:
    plant_summ = pd.concat([get_inepa_plantsumm(eia_do_dp,g), plant_summ])
plant_summ


# %%
# 4. Summarize at the generator level
# helper summary function
def get_inepa_gensumm(df, group):
    summ = (df.groupby([group, 'gen_in_epa'])
                .agg({'plant_code':'nunique', 'key_gen':['nunique'], 'nameplate_capacity_mw':['sum', 'mean', 'std'], 'gen_age_yrs':['mean', 'std']})
                .reset_index())
    summ.rename(columns={group:'value', 'key_gen':'generators', 'plant_code':'plants'}, inplace=True)
    return summ

# drop nonrenewables
print('Dropping nonrenewable generators:\n', eia_do.groupby('not_renewable')['key_gen'].nunique())
eia_do_dr = eia_do.loc[eia_do.not_renewable]
# summarize
gen_summ = pd.DataFrame([])
groups = ['year', 'entity_type', 'state_plant']
for g in groups:
    gen_summ = pd.concat([get_inepa_gensumm(eia_do_dr,g), gen_summ])
gen_summ

# %%
