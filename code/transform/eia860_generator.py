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
# GENERATOR DATA
readin_dict = {}
readin_dict[END_YEAR] = {
    'files':        [
            f'{END_YEAR}/3_1_Generator_Y{str(END_YEAR)}.xlsx', 
            f'{END_YEAR}/3_2_Wind_Y{END_YEAR}.xlsx', 
            f'{END_YEAR}/3_3_Solar_Y{END_YEAR}.xlsx', 
            f'{END_YEAR}/3_4_Energy_Storage_Y{END_YEAR}.xlsx',
            f'{END_YEAR}/3_5_Multifuel_Y{END_YEAR}.xlsx'
    ],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  None
}
# cut corner: 2013/2016 is mostly the same
for year in range(2016, END_YEAR+1):
    readin_dict[year] = readin_dict[END_YEAR].copy()
    if year >= 2016:
        readin_dict[year]['files'] = [
                f'{year}/3_1_Generator_Y{str(year)}.xlsx', 
                f'{year}/3_2_Wind_Y{year}.xlsx', 
                f'{year}/3_3_Solar_Y{year}.xlsx', 
                f'{year}/3_4_Energy_Storage_Y{year}.xlsx',
                f'{year}/3_5_Multifuel_Y{year}.xlsx'
        ]
    else:
        readin_dict[year]['files'] = [
                f'{year}/3_1_Generator_Y{str(year)}.xlsx',
                f'{year}/3_2_Wind_Y{year}.xlsx',
                f'{year}/3_3_Solar_Y{year}.xlsx',
                f'{year}/3_4_Multifuel_Y{year}.xlsx']

readin_dict[2012] = {
    'files':        [f'{2012}/GeneratorY{2012}.xlsx', f'{2012}/MultifuelY{2012}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  {'nameplate':'nameplate_capacity_mw', 'sector_number':'sector'}
}
readin_dict[2011] = {
    'files':        [f'{2011}/GeneratorY{2011}.xlsx', f'{2011}/MultifuelY{2011}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  {'nameplate':'nameplate_capacity_mw', 'sector_number':'sector'}
}
readin_dict[2010] = {
    'files':        [f'{2010}/GeneratorsY{2010}.xls', f'{2010}/MultiFuelY{2010}.xls'],
    'excel_params': {'header':0, 'sheet_name':None},
    'rename_vars':  {'nameplate':'nameplate_capacity_mw', 'sector_number':'sector'}
}
readin_dict[2009] = {
    'files':        [f'{2009}/GeneratorY{str(2009)[2:]}.xls', f'{2009}/MultiFuelY{str(2009)[2:]}.xls'],
    'excel_params': {'header':0, 'sheet_name':None},
    'rename_vars':  {'nameplate':'nameplate_capacity_mw', 'sector_number':'sector'}
}
readin_dict[2008] = {
    'files':        [
        f'{2008}/GenY{str(2008)[2:]}.xls', f'{2008}/MFExistY{str(2008)[2:]}.xls',
        f'{2008}/MFPropY{str(2008)[2:]}.xls', f'{2008}/PRGenY{str(2008)[2:]}.xls'],
    'excel_params': {'header':0, 'sheet_name':None},
    'rename_vars':  {
           'utilcode':'utility_id', 'plntcode':'plant_code', 'gencode':'generator_id',
            'owner':'ownership', 'primemover':'prime_mover', 'nameplate':'nameplate_capacity_mw',
            'insvmonth':'operating_month', 'insvyear':'operating_year',
            'retiremonth':'retirement_month', 'retireyear':'retirement_year',
            'prop_cofire_energy_source_1':'cofire_energy_source_1',
            'proposed_nameplate':'nameplate_capacity_mw', 'proposed_energy_source_1':'energy_source_1'
            }
}
readin_dict[2007] = {
    'files':        [
        f'{2007}/GenY{str(2007)[2:]}.xls', f'{2007}/MFExistY{str(2007)[2:]}.xls',
        f'{2007}/MFPropY{str(2007)[2:]}.xls', f'{2007}/PRGenY{str(2007)[2:]}.xls'],
    'excel_params': {'header':0, 'sheet_name':None},
    'rename_vars':  {
           'utilcode':'utility_id', 'plntcode':'plant_code', 'gencode':'generator_id',
            'owner':'ownership', 'primemover':'prime_mover', 'nameplate':'nameplate_capacity_mw',
            'insvmonth':'operating_month', 'insvyear':'operating_year',
            'retiremonth':'retirement_month', 'retireyear':'retirement_year',
            'prop_cofire_energy_source_1':'cofire_energy_source_1',
            'proposed_nameplate':'nameplate_capacity_mw', 'proposed_energy_source_1':'energy_source_1'
            }
}
readin_dict[2006] = {
    'files':        [
        f'{2006}/GenY{str(2006)[2:]}.xls', f'{2006}/MFExistY{str(2006)[2:]}.xls',
        f'{2006}/MFPropY{str(2006)[2:]}.xls', f'{2006}/PRGenY{str(2006)[2:]}.xls'],
    'excel_params': {'header':0, 'sheet_name':None},
    'rename_vars':  {
           'utilcode':'utility_id', 'plntcode':'plant_code', 'gencode':'generator_id',
            'owner':'ownership', 'primemover':'prime_mover', 'nameplate':'nameplate_capacity_mw',
            'insvmonth':'operating_month', 'insvyear':'operating_year',
            'retiremonth':'retirement_month', 'retireyear':'retirement_year',
            'prop_cofire_energy_source_1':'cofire_energy_source_1',
            'proposed_nameplate':'nameplate_capacity_mw', 'proposed_energy_source_1':'energy_source_1'
            }
}

# %%
if __name__ == '__main__':
    vars_date = ['operating_month', 'operating_year',
             'current_month', 'current_year',
             'planned_retirement_month', 'planned_retirement_year',
             'retirement_month', 'retirement_year']
    vars_keep = ['utility_id', 'plant_code', 'generator_id', 'status', 'ownership', 'sector', 
                'energy_source_1', 'cofire_energy_source_1', 
                'prime_mover', 'nameplate_capacity_mw'] + vars_date
    gdf = readin_eia_years(f'{PATH_RAW}eia/f860/', readin_dict, START_YEAR)
    # UPDATE DATA TYPES AND SAVE INTERMEDIATE
    gdf['utility_id'] = pd.to_numeric(gdf.utility_id, errors='coerce').astype('Int64')
    gdf = gdf.loc[gdf.utility_id.notna()]
    gdf['plant_code'] = pd.to_numeric(gdf.plant_code, errors='coerce').astype('Int64')
    gdf['generator_id'] = gdf.generator_id.astype(str)
    gdf['sector'] = pd.to_numeric(gdf.sector, errors='coerce').astype('Int64')
    gdf['nameplate_capacity_mw'] = pd.to_numeric(gdf.nameplate_capacity_mw, errors='coerce').astype('Float64')
    for var in gdf.columns.intersection(vars_date):
        gdf[var] = pd.to_numeric(gdf[var], errors='coerce').astype('Int64')
    gdf = gdf.astype({'status':str, 'prime_mover':str})
    # fix datatypes
    col_str = gdf.columns[gdf.dtypes == 'object']
    gdf[col_str] = gdf[col_str].astype(str)
    gdf.to_parquet(PATH_INTERIM + 'eia860_generator.parquet', index=False)

    # DEDUP GENERATORS
    gdf['dup_key'] = gdf[['utility_id', 'plant_code', 'generator_id', 'year']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    gdf['duplicate'] = gdf.dup_key.isin(gdf.loc[gdf.dup_key.duplicated(), 'dup_key'].drop_duplicates())
    print('total dups:', gdf.duplicate.sum())
    # DECISION: Drop duplicates duplicated across sheets, taking info from the main, "generator" sheet (and "proposed gen" in 06-08)
    gdf['dup_ingen'] = gdf.file.str.startswith(('gen', 'prgen'))
    gdf['dup_anyingen'] = gdf.groupby('dup_key')['dup_ingen'].transform('sum')
    gdf['dup_keep'] = True
    gdf.loc[~gdf.dup_ingen & gdf.dup_anyingen, 'dup_keep'] = False
    gdf_dedup = gdf.loc[gdf.dup_keep]
    gdf_dedup['dup_key'] = gdf_dedup[['utility_id', 'plant_code', 'generator_id', 'year']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    gdf_dedup['duplicate'] = gdf_dedup.dup_key.isin(gdf_dedup.loc[gdf_dedup.dup_key.duplicated(), 'dup_key'].drop_duplicates())
    print('total dups:', gdf_dedup.duplicate.sum())
    gdf_dedup.drop(columns=['dup_key', 'dup_ingen', 'dup_anyingen', 'dup_keep'], inplace=True)
    # WRITE TO FILE
    # write
    gdf_dedup.to_parquet(PATH_PROCESSED + 'eia860_generator.parquet', index=False)
    # summarize
    gdf_dedup['pid'] = gdf_dedup.utility_id.astype(str) + '.' + gdf_dedup.plant_code.astype(str)
    gdf_dedup['gid'] = gdf_dedup.pid + '.' + gdf_dedup.generator_id
    gdf_dedup = gdf_dedup.rename(columns={'utility_id':'uid'})
    print('Generator dataset:')
    print(summarize_id_counts_byyear(gdf_dedup.copy(), ['uid', 'pid', 'gid']))


# %%
