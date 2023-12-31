# %%
import sqlite3
import pandas as pd
import numpy as np
import os
from pathlib import Path
from datetime import datetime as dt
pd.set_option('display.max_columns', None)

YR_START = 2018

# get paths
PATH_DATA = '../../data/pudl/'
path = Path(PATH_DATA + 'md5sums.txt').read_text()
path = path[path.find('pudl'):path.find('.tzg')-4]
PATH_PUDL = PATH_DATA + path + '/pudl_data/'
PATH_PUDL_SQL = PATH_PUDL + 'sqlite/pudl.sqlite'

# %%
def get_sqltable(table_name:str, where:str=None, path=PATH_PUDL_SQL) -> pd.DataFrame:
    conn = sqlite3.connect(path)
    cursor = conn.cursor()
    if where is None:
        cursor.execute(f'select * from {table_name}')
    else:
        cursor.execute(f'select * from {table_name} where {where}')
    df = pd.DataFrame(cursor.fetchall(), 
                      columns=[column[0] for column in cursor.description])
    if 'report_date' in df.columns:
        df['year'] = pd.to_datetime(df.report_date).dt.year
    conn.close()
    return df

# %%
# GET TABLE NAMES
table_names = get_sqltable('sqlite_master', "type='table'").name
table_names

# %%
# SAVE USEFUL KEYS
key_energysource = get_sqltable('energy_sources_eia')
key_fueltype = get_sqltable('fuel_types_aer_eia')
key_status = get_sqltable('operational_status_eia')
key_primemover = get_sqltable('prime_movers_eia')
key_opstatus = get_sqltable('operational_status_eia')



# %% 
# GET ALL GENERATORS AT PLANTS WITH NATURAL GAS GENERATORS
# get all generator data
gen_yrs = get_sqltable('generators_eia860', f"report_date >= '{YR_START}-01-01'")
gen_ent = get_sqltable('generators_entity_eia')
gens = pd.merge(left=gen_yrs, right=gen_ent, on=['plant_id_eia', 'generator_id'])
assert len(gen_yrs) == len(gens)
# get gas generators
# gens = gens.loc[(gens.operational_status == 'existing')] # existing
gens['is_gas'] = ((gens.fuel_type_code_pudl == 'gas')) # & 
                #   gens.prime_mover_code.isin(['CA', 'CC', 'CS', 'CT', 'GT', 'IC', 'ST']))
# gens['is_gas'] = gens.fuel_type_code_pudl == 'gas'
# gens['is_gas'] = gens.energy_source_code_1 == 'NG'
gens_gas = gens.loc[gens.is_gas]
gens_at_gasplant = gens.loc[gens.plant_id_eia.isin(gens_gas.plant_id_eia)] # at gas plant
gas_plant_ids = gens.plant_id_eia.drop_duplicates().values

# %%
# GET ALL PLANT DATA FOR PLANTS WITH GAS GENERATORS
plant_yrs = get_sqltable('plants_eia860', f"report_date >= '{YR_START}-01-01'")
plant_ent = get_sqltable('plants_entity_eia')
plants = pd.merge(left=plant_yrs, right=plant_ent, on=['plant_id_eia'])
assert len(plant_yrs) == len(plants)
gas_plant = plants.loc[plants.plant_id_eia.isin(gas_plant_ids)]

# %%
# GET ALL GENERATION FROM THESE PLANTS AND THESE GENERATORS
# GENERATION-LEVEL DATA: NOTE: `pwr_gen` is unique on plantID, genID, and year
pwr_gen = get_sqltable('generation_eia923', f"report_date >= '{YR_START}-01-01'")
pwr_gen['month'] = pd.to_datetime(pwr_gen.report_date).dt.month
# PLANT-LEVEL DATA: NOTE: `pwr_plant` is unique on plantID, year, energy source, and prime mover
pwr_plant = get_sqltable('generation_fuel_eia923', f"report_date >= '{YR_START}-01-01'")
pwr_plant['month'] = pd.to_datetime(pwr_gen.report_date).dt.month


# %%
# MERGE ANNUAL GENERATION STATS ONTO GENERATOR INFO
pwr_gen_ann = pwr_gen.groupby(['plant_id_eia', 'generator_id', 'year']).agg({'net_generation_mwh':'sum'})
print('generators at gasplant:\t', len(gens_at_gasplant))
gens_at_gasplant_m = pd.merge(left=gens_at_gasplant, right=pwr_gen_ann, how='left', on=['plant_id_eia', 'generator_id', 'year'])
print('after merge:\t', len(gens_at_gasplant_m))

# %%
# MERGE PLANT AND GENERATOR DATA
m = pd.merge(left=gas_plant, right=gens_at_gasplant_m, on=['plant_id_eia', 'utility_id_eia', 'year', 'report_date'])
print('generator rows:\t', len(gens_at_gasplant_m))
print('merge rows:\t', len(m))
# assert len(m) == len(gens_at_gasplant)

# %%
# GET GENERATOR-LEVEL SUMMARY (to compare to EIA)
# make additional columns
bucket_cap = {0:'0-24', 24:'25-49', 49:'50-99', 99:'100-149', 149:'150-249', 
              249:'250-499', 499:'500-749', 749:'750-999', 999:'1000-1500', 1500:'1500+'}
m['capacity_mw_bucket'] = pd.cut(m.capacity_mw, bins=list(bucket_cap.keys()), 
                                 labels=list(bucket_cap.values())[:-1], right=True, include_lowest=True)
m['age'] = (pd.to_datetime('today') - pd.to_datetime(m.operating_date)).dt.days // 365
m['capacity_factor'] = (m.net_generation_mwh / (m.capacity_mw*365*24)).clip(0, 1)
m.loc[m.capacity_factor.isna(), 'capacity_factor'] = 0.0
m['total'] = 'total'
# subset dataset
m_sub = m.loc[(m.year == 2021) & (m.is_gas)]
# sumarize dataset
summ_dict = {'generator_id':['count'], 
                  'age':['mean'],
                  'summer_capacity_mw':['mean', 'sum'],
                  'capacity_factor':['mean']}
m_summ = m_sub.groupby(['capacity_mw_bucket']).agg(summ_dict)
m_summ = pd.concat([m_summ, m_sub.groupby('total').agg(summ_dict)], axis=0)
for col in [('generator_id', 'count'), ('summer_capacity_mw', 'sum')]:
    m_summ[(col[0], 'pct_total')] = m_summ[col] / m_summ.loc['total', col].sum()

m_summ.round(3)

# %%
# GET GENERATOR VALUE BY CAPACITY FACTOR AND SIZE
bucket_cap = np.arange(0, 1600, 100)
bucket_cap_labs = [
    str(round(bot)) + '-' + str(round(top-1)) for bot, top in 
    zip(bucket_cap[:-1], bucket_cap[1:])]
bucket_capfac = np.arange(0, 1.1, 0.1)
bucket_capfac_labs = [
    str(round(bot, 3)) + '-' + str(round(top-0.01, 2)) for bot, top in 
    zip(bucket_capfac[:-1], bucket_capfac[1:])]
m_sub['capacity_mw_bucket'] = pd.cut(
    m_sub.capacity_mw, bins=bucket_cap, labels=bucket_cap_labs, 
    right=True, include_lowest=True)
m_sub['capfac_bucket'] = pd.cut(
    m_sub.capacity_factor, bins=bucket_capfac, labels=bucket_capfac_labs, 
    right=True, include_lowest=True)

# summarize raw
m_summ = m_sub.groupby(['capacity_mw_bucket', 'capfac_bucket']).agg({'generator_id':'count', 'summer_capacity_mw':'sum', 'net_generation_mwh':'sum'}).reset_index()
var = 'generator_id'
m_summ = m_summ.pivot(index='capacity_mw_bucket', columns='capfac_bucket', values=var)
m_summ.sort_index(ascending=False, inplace=True)
m_summ = m_summ[m_summ.columns.sort_values(ascending=False).values]
# summarize cummulative
m_summ_cumm = m_summ.cumsum()
m_summ_cumm = m_summ_cumm.cumsum(axis=1)
m_summ_cumm_pct = m_summ_cumm / m_summ_cumm.iloc[-1, -1]
display(m_summ.round(3))
display(m_summ_cumm_pct.round(3))

# %%
# PLANT-LEVEL SUMMARY
# By year
#    - Capacity (% of capacity that is natural gas)
#    - Generation (% of generation from natural gas)
#    - Capacity factor
#    - Emissions
#    - Num generators, average age generators
