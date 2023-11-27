# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import re

from utils_transform import readin_eia_years
from utils_summ import summarize_id_counts_byyear

def rename_cols(s):
    sn = s
    sn = re.sub(r'(mmbtus[_]?)(?!ep)', 'mmbtu', sn)
    sn = re.sub(r'(mmbtu[_]?)(?=[^\s$])', 'mmbtu_', sn)
    sn = re.sub(r'(mmbtu[_]?per[_]?unit[_]?)', 'mmbtu_per_unit_', sn)
    sn = re.sub(r'(net[_]?generation[_]?)', 'net_gen_', sn)
    sn = re.sub(r'(net[_]?gen[_]?)', 'net_gen_', sn)
    sn = re.sub('net_gen_year_to_date|net_gen_megawatthours', 'total_net_gen_mwh', sn)
    sn = re.sub('elec_fuel_consumption_mmbtu', 'total_elec_fuel_consumption_mmbtu', sn)
    sn = re.sub('electric_fuel_consumption_quantity', 'total_elec_fuel_consumption_quantity', sn)
    sn = re.sub('eia_sector_number', 'sector_number', sn)
    sn = re.sub('prime_mover_type', 'reported_prime_mover', sn)

    months = ['january', 'february', 'march', 'april', 'may', 'june', 
            'july', 'august', 'september', 'october', 'november', 'december']
    for m in months:
        sn = re.sub(f'({m})|({m[:3]})', m[:3], sn)
    return sn

# %%
# EIA form 923 OPERATIONS DATA
# Data documentation: https://www.eia.gov/electricity/data/eia923/
readin_dict = {}
readin_dict[2021] = {
    'files': [f'{2021}/EIA923_Schedules_2_3_4_5_M_12_{2021}_Final_Revision.xlsx'],
    'excel_params':{'header':5, 'na_values': '.',
                    'sheet_name':['Page 1 Generation and Fuel Data', 'Page 4 Generator Data']},
    'rename_vars': rename_cols
}
# cut corner: 2011+ is mostly the same
for yr in range(2011, 2021):
    readin_dict[yr] = readin_dict[2021].copy()
    readin_dict[yr]['files'] = [f'{yr}/EIA923_Schedules_2_3_4_5_M_12_{yr}_Final_Revision.xlsx']

readin_dict[2013]['files'] = [f'{2013}/EIA923_Schedules_2_3_4_5_{2013}_Final_Revision.xlsx']
readin_dict[2011]['files'] = [f'{2011}/EIA923_Schedules_2_3_4_5_{2011}_Final_Revision.xlsx']

readin_dict[2010] = {
    'files': [f'{2010}/EIA923 SCHEDULES 2_3_4_5 Final {2010}.xls'],
    'excel_params': {'header':7, 'na_values':'.',
                     'sheet_name':['Page 1 Generation and Fuel Data', 'Page 4 Generator Data']},
    'rename_vars': rename_cols
}
readin_dict[2009] = {
    'files': [f'{2009}/EIA923 SCHEDULES 2_3_4_5 M Final {2009} REVISED 05252011.xls'],
    'excel_params': {'header':7, 'na_values':'.',
                     'sheet_name':['Page 1 Generation and Fuel Data', 'Page 4 Generator Data']},
    'rename_vars': rename_cols
}
readin_dict[2008] = {
    'files': [f'{2008}/eia923December{2008}.xls'],
    'excel_params': {'header':7, 'na_values':'.',
                     'sheet_name':['Page 1 Generation and Fuel Data', 'Page 4 Generator Data']},
    'rename_vars': rename_cols
}
readin_dict[2007] = {
    'files': [f'{2007}/f906920_{2007}.xls'],
    'excel_params': {'header':7, 'na_values':'.',
                     'sheet_name':['Page 1 Generation and Fuel Data']},
    'rename_vars': rename_cols
}
readin_dict[2006] = {
    'files': [f'{2006}/f906920_{2006}.xls'],
    'excel_params': {'header':7, 'na_values':'.',
                     'sheet_name':['Page 1 Generation and Fuel Data']},
    'rename_vars': rename_cols
}



# %%
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('transform_eia923_ops')

    year_start = snakemake.params.year_start
    year_end = snakemake.params.year_end
    path_raw = snakemake.params.indir
    outfile = snakemake.output.outfile

    # read-in parameters
    print('Reading in data...')
    vars_keep = []
    df_raw = readin_eia_years(path_raw, readin_dict, year_start)
    # drop state-level fuel increments
    df = df_raw.loc[(df_raw.plant_id != 99999) & (df_raw.plant_name != 'State-Fuel Level Increment')].copy()
    print(df.groupby(['year', 'sheet']).agg({'file':'count'}))

    # confirm unique observation ID
    print('Confirming unique observation ids...')
    plant_groupcols = ['year', 'plant_id', 'nuclear_unit_id', 'reported_prime_mover', 
                       'reported_fuel_type_code', 'combined_heat_and_power_plant']
    print('Plant-level duplicates by groupcols:\n', plant_groupcols)
    print(df
        .loc[df.sheet == 'page_1_generation_and_fuel_data']
        .groupby(plant_groupcols, dropna=False)
        .agg({'file':'count'}).value_counts()
    )
    gen_groupcols = ['year', 'plant_id', 'combined_heat_and_power_plant', 'generator_id']
    print('Generator-level duplicates by groupcols:\n', gen_groupcols)
    print(df
        .loc[df.sheet != 'page_1_generation_and_fuel_data']
        .groupby(gen_groupcols, dropna=False)
        .agg({'file':'count'}).value_counts()
    )

    # get annual values for each row
    print('Cleaning up totals columns...')
    vars_tot_dict = {
        'total_elec_fuel_consumption_mmbtu':'elec_mmbtu', 
        'total_elec_fuel_consumption_quantity':'elec_quantity',
        'total_fuel_consumption_mmbtu':'tot_mmbtu',
        'total_fuel_consumption_quantity':'quantity',
        'total_net_gen_mwh':'net_gen'}
    for v, pre in vars_tot_dict.items():
        month_cols = [col for col in df.columns if col.startswith(pre)]
        mask_allna = df[month_cols].isna().all(1)
        mask_totna = df[v].isna()
        df.loc[~mask_allna, f'{pre}_tot_mo'] = df.loc[~mask_allna, month_cols].sum(1)
        df.rename(columns={v:f'{pre}_tot_an'}, inplace=True)

    # write to file
    df = (df
          .drop(columns=df.columns.intersection(['reserved1', 'reserved2', 'reserved_1', 'reserved_2', 'reserved']))
          .astype({'generator_id':str}))
    df.to_parquet(outfile, index=False)
    
    # transpose wide to long
    # print('Transposing wide to long...')
    # vars_id = ['year', 'file', 'sheet', 'operator_id', 'plant_id', 
    #            'generator_id', 'nuclear_unit_id',
    #            'combined_heat_and_power_plant',
    #            'reported_prime_mover', 'reported_fuel_type_code']
    # vars_val_pre = ('net_generation', 'elec_mmbtu', 'netgen')
    # vars_val = [col for col in df.columns if col.startswith(vars_val_pre)]
    # dfl = df.melt(id_vars=vars_id, value_vars=vars_val, var_name='variable_full')
    # dfl['date'] = pd.to_datetime(dfl.year.astype(str) + 
    #                             dfl['variable_full'].str.split('_').str[-1],
    #                             format='%Y%B', errors='coerce')
    # dfl['variable'] = dfl['variable_full'].str.split('_').str[:-1].str.join('_')
    # dfl = dfl.loc[dfl.date.notna()].drop(columns='variable_full')
    # dfl.loc[dfl.variable == 'netgen', 'variable'] = 'net_generation'

    # summarize unique ids over time
    print('Summarizing unique IDs over time...')
    df.loc[df.sheet == 'page_1_generation_and_fuel_data', 'puid'] = (
        df.plant_id.astype(str) + df.reported_prime_mover + 
        df.reported_fuel_type_code + df.combined_heat_and_power_plant)
    df.loc[df.sheet == 'page_4_generator_data', 'gid'] = (
        df.plant_id.astype(str) + df.generator_id.astype(str))
    print(summarize_id_counts_byyear(df.rename(columns={'plant_id':'pid'}), ['pid', 'puid', 'gid']))
