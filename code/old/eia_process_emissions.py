# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

from utils import PATH_EIA, PATH_PROCESSED, START_YEAR, END_YEAR
from utils_data_eia import readin_eia

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 70)

# %%
# READIN EMISSIONS DATA
vars_keep = ['plant_code', 'prime_mover', 'fuel_code', 'state', 'sector',
    'generation_kwh', 'total_fuel_consumption_mmbtu', 
    'fuel_consumption_for_electric_generation_mmbtu', 'tons_of_co2_emissions',
    'nerc_region', 'balancing_authority_code',
    'eia_balancing_authority_region', 'year']
readin_dict={year:{} for year in range(START_YEAR, END_YEAR+1)}
for year in readin_dict.keys():
    readin_dict[year]['vars_keep'] = vars_keep
    readin_dict[year]['path_file'] = f'emissions{year}.xlsx'
    readin_dict[year]['excel_params'] = {'header':1}
    readin_dict[year]['rename_vars'] = {'state':'state_plant', 'sector_code':'sector'}

edf = readin_eia(path_folder=f'{PATH_EIA}emissions/', readin_dict=readin_dict)
# %%
edf['plant_code'] = pd.to_numeric(edf.plant_code, errors='coerce').astype('Int64')
edf = edf.loc[edf.plant_code.notna() & (edf.plant_code != 99999)]
edf.to_parquet(PATH_PROCESSED + 'eia_emissions.parquet', index=False)
# %%
edf.loc[edf.plant_code == 6002]
# %%
