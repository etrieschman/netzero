# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

from utils_transform import (
    PATH_RAW, PATH_INTERIM, PATH_PROCESSED, START_YEAR, END_YEAR)
from utils_transform import readin_eia_years
from utils_summ import summarize_id_counts_byyear

# %%
# EIA form 923 OPERATIONS DATA
# Data documentation: https://www.eia.gov/electricity/data/eia923/
readin_dict = {}
readin_dict[2021] = {
    'files': [f'{2021}/EIA923_Schedules_2_3_4_5_M_12_{2021}_Final_Revision.xlsx'],
    'excel_params':{'header':5, 'na_values': '.',
                    'sheet_name':['Page 1 Generation and Fuel Data', 'Page 4 Generator Data']},
    'rename_vars': None
}
# cut corner: 2013+ is all the same
for yr in range(2011, 2021):
    readin_dict[yr] = readin_dict[2021].copy()
    readin_dict[yr]['files'] = [f'{yr}/EIA923_Schedules_2_3_4_5_M_12_{yr}_Final_Revision.xlsx']

readin_dict[2013]['files'] = [f'EIA923_Schedules_2_3_4_5_{2013}_Final_Revision.xlsx']
readin_dict[2011]['files'] = [f'EIA923_Schedules_2_3_4_5_{2011}_Final_Revision.xlsx']

readin_dict[2010] = {
    'files': [f'EIA923 SCHEDULES 2_3_4_5 Final {2010}.xls'],
    'excel_params': {'header':7, 'na_values':'.',
                     'sheet_name':['Page 1 Generation and Fuel Data', 'Page 4 Generator Data']},
    'rename_vars': None
}
# I AM HERE
readin_dict[2009] = {
    'files': [f'SCHEDULE 3A 5A 8A 8B 8C 8D 8E 8F REVISED {2009} 04112011.xls'],
    'excel_params': {'header':7, 'na_values':'.',
                     'sheet_name':['Page 1 Generation and Fuel Data', 'Page 4 Generator Data']},
    'rename_vars': None
}

# %%
if __name__ == '__main__':
    # read-in parameters
    print('Reading in data...')
    vars_keep = []
    df_raw = readin_eia_years(f'{PATH_RAW}eia/f923/', readin_dict, 2014)
    # drop state-level fuel increments
    df = df_raw.loc[(df_raw.plant_id != 99999) & (df_raw.plant_name != 'State-Fuel Level Increment')].copy()
    display(df.groupby(['year', 'sheet']).agg({'file':'count'}))

    # confirm unique observation ID
    print('Confirming unique observation ids...')
    plant_groupcols = ['year', 'plant_id', 'nuclear_unit_id', 'reported_prime_mover', 
                       'reported_fuel_type_code', 'combined_heat_and_power_plant']
    print('Plant-level duplicates by groupcols:\n', plant_groupcols)
    display(df
        .loc[df.sheet == 'page_1_generation_and_fuel_data']
        .groupby(plant_groupcols, dropna=False)
        .agg({'file':'count'}).value_counts()
    )
    gen_groupcols = ['year', 'plant_id', 'combined_heat_and_power_plant', 'generator_id']
    print('Generator-level duplicates by groupcols:\n', gen_groupcols)
    display(df
        .loc[df.sheet != 'page_1_generation_and_fuel_data']
        .groupby(gen_groupcols, dropna=False)
        .agg({'file':'count'}).value_counts()
    )

    # transpose wide to long
    print('Transposing wide to long...')
    vars_id = ['year', 'file', 'sheet', 'operator_id', 'plant_id', 
               'generator_id', 'nuclear_unit_id',
               'combined_heat_and_power_plant',
               'reported_prime_mover', 'reported_fuel_type_code']
    vars_val_pre = ('net_generation', 'elec_mmbtu', 'netgen')
    vars_val = [col for col in df.columns if col.startswith(vars_val_pre)]
    dfl = df.melt(id_vars=vars_id, value_vars=vars_val, var_name='variable_full')
    dfl['date'] = pd.to_datetime(dfl.year.astype(str) + 
                                dfl['variable_full'].str.split('_').str[-1],
                                format='%Y%B', errors='coerce')
    dfl['variable'] = dfl['variable_full'].str.split('_').str[:-1].str.join('_')
    dfl = dfl.loc[dfl.date.notna()].drop(columns='variable_full')
    dfl.loc[dfl.variable == 'netgen', 'variable'] = 'net_generation'

    # summarize unique ids over time
    print('Summarizing unique IDs over time...')
    dfl.loc[dfl.sheet == 'page_1_generation_and_fuel_data', 'puid'] = (
        dfl.plant_id.astype(str) + dfl.reported_prime_mover + 
        dfl.reported_fuel_type_code + dfl.combined_heat_and_power_plant)
    dfl.loc[dfl.sheet == 'page_4_generator_data', 'gid'] = (
        dfl.plant_id.astype(str) + dfl.generator_id.astype(str))
    summarize_id_counts_byyear(dfl.rename(columns={'plant_id':'pid'}), ['pid', 'puid', 'gid'])
# %%
