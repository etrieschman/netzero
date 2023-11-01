# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

from utils_transform import (
    PATH_RAW, PATH_INTERIM, PATH_PROCESSED, START_YEAR, END_YEAR)
from utils_transform import readin_eia_years
from utils_summ import summarize_id_counts_byyear

def make_dates(df, vm, vy):
    dates = pd.to_datetime(df[vy].astype(str) + df[vm].astype(str), 
        format='%Y%m', errors='coerce')
    return dates

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
for year in range(2013, END_YEAR+1):
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

vars_date = ['operating_month', 'operating_year',
            'current_month', 'current_year',
            'retirement_month', 'retirement_year',
            'effective_month', 'effective_year',
            'planned_retirement_month', 'planned_retirement_year']
vars_keep = vars_date + [
 'year', 
 'utility_id', 'plant_code', 'generator_id',
 'status', 'ownership',
 'technology', 'prime_mover',
 'sector',
 'nameplate_capacity_mw', 'summer_capacity_mw', 'winter_capacity_mw',
 'multiple_fuels', 'cofire_fuels', 'associated_with_combined_heat_and_power_system',
 'energy_source_1', 'energy_source_2', 'energy_source_3',
 'energy_source_4', 'energy_source_5', 'energy_source_6',
 'cofire_energy_source_1', 'cofire_energy_source_2',
 'cofire_energy_source_3', 'cofire_energy_source_4',
 'cofire_energy_source_5', 'cofire_energy_source_6',
 'startup_source_1', 'startup_source_2',
 'startup_source_3', 'startup_source_4',
 'sheet', 'file']

# %%
if __name__ == '__main__':
    print('Reading in data...')
    gdf = readin_eia_years(f'{PATH_RAW}eia/f860/', readin_dict, START_YEAR, keep_vars=None)
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

    # restrict variables
    new_cols = ['year', 'file', 'sheet']
    gdf = gdf.drop(columns=gdf.columns.difference(vars_keep + new_cols))

    # DEDUP GENERATORS
    print('Deduping generators...')
    gdf['dup_key'] = gdf[['utility_id', 'plant_code', 'generator_id', 'year']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    gdf['duplicate'] = gdf.dup_key.isin(gdf.loc[gdf.dup_key.duplicated(), 'dup_key'].drop_duplicates())
    print('total dups:', gdf.duplicate.sum())
    # DECISION: Drop duplicates duplicated across sheets, taking info from the main, "generator" sheet (and "proposed gen" in 06-08)
    gdf['dup_ingen'] = gdf.file.str.contains(('gen'))
    gdf['dup_anyingen'] = gdf.groupby('dup_key')['dup_ingen'].transform('sum')
    gdf['dup_keep'] = True
    gdf.loc[~gdf.dup_ingen & gdf.dup_anyingen, 'dup_keep'] = False
    gdf_dedup = gdf.loc[gdf.dup_keep]
    gdf_dedup['dup_key'] = gdf_dedup[['utility_id', 'plant_code', 'generator_id', 'year']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    gdf_dedup['duplicate'] = gdf_dedup.dup_key.isin(gdf_dedup.loc[gdf_dedup.dup_key.duplicated(), 'dup_key'].drop_duplicates())
    print('total dups:', gdf_dedup.duplicate.sum())
    gdf_dedup.drop(columns=['dup_key', 'dup_ingen', 'dup_anyingen', 'dup_keep', 'duplicate'], inplace=True)

    print('Cleaning up date columns...')
    # UPDATE COLUMNS: OPERATION STATUS
    # FINDING: MISSING OR UNKOWN STATUS ONLY IN 2006
    print('unknown statuses:')
    print(gdf_dedup.loc[gdf_dedup.status == 'nan'].groupby(['sheet', 'file'])['utility_id'].count())
    # DECISION: These all look like 2006 proposed generators. updating status to `P`
    gdf_dedup.loc[gdf_dedup.status == 'nan', 'status'] = 'IP'
    mask_op = (gdf_dedup.status == 'OP')
    mask_pl = (gdf_dedup.status.isin(('TS', 'P', 'L', 'T', 'U', 'V')))
    mask_sb = (gdf_dedup.status.isin(('SB', 'OA', 'OS', 'BU'))) # BU is from 2006 only, stands for "backup"
    mask_ds = (gdf_dedup.status.isin(('RE', 'CN', 'IP', 'OT')))

    # UPDATE COLUMNS: OPERATION DATES
    vars_month = ['operating_month', 'current_month', 'retirement_month', 'planned_retirement_month']
    vars_year = ['operating_year', 'current_year', 'retirement_year', 'planned_retirement_year']
    for vm in vars_month:
        gdf_dedup.loc[(gdf_dedup[vm] > 12) | (gdf_dedup[vm] < 1), vm] = 1
    for vm, vy in zip(vars_month, vars_year):
        gdf_dedup.loc[gdf_dedup[vy].notna() & gdf_dedup[vm].isna(), vm] = 1
        gdf_dedup.loc[gdf_dedup[vy] == 0, vy] = pd.NA
        gdf_dedup.loc[gdf_dedup[vy].isna(), vm] = pd.NA

    gdf_dedup.loc[mask_op | mask_sb | mask_ds , 'dt_operation_start'] = (
        make_dates(gdf_dedup.loc[mask_op | mask_sb | mask_ds], 'operating_month', 'operating_year'))
    gdf_dedup.loc[mask_op | mask_sb, 'dt_operation_end'] = (
        make_dates(gdf_dedup.loc[mask_op | mask_sb], 'planned_retirement_month', 'planned_retirement_year'))
    gdf_dedup.loc[mask_pl, 'dt_operation_start'] = (
        make_dates(gdf_dedup.loc[mask_pl], 'current_month', 'current_year'))
    gdf_dedup.loc[mask_ds, 'dt_operation_end'] = (
        make_dates(gdf_dedup.loc[mask_ds], 'retirement_month', 'retirement_year'))
    
    # WRITE TO FILE
    print('Writing to file...')
    gdf_dedup.to_parquet(PATH_PROCESSED + 'eia860_generator.parquet', index=False)
    
    # summarize
    print('Summarizing unique IDs over time...')
    gdf_dedup['pid'] = gdf_dedup.utility_id.astype(str) + '.' + gdf_dedup.plant_code.astype(str)
    gdf_dedup['gid'] = gdf_dedup.pid + '.' + gdf_dedup.generator_id
    gdf_dedup = gdf_dedup.rename(columns={'utility_id':'uid'})
    print('Generator dataset:')
    print(summarize_id_counts_byyear(gdf_dedup.copy(), ['uid', 'pid', 'gid']))


# %%
