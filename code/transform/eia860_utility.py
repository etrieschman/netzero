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
# UTILITY DATA
readin_dict = {}
readin_dict[END_YEAR] = {
    'files': [f'{END_YEAR}/1___Utility_Y{END_YEAR}.xlsx'],
    'excel_params':{'header':1, 'sheet_name':None},
    'rename_vars':{'city':'city_util', 'state':'state_util', 'zip':'zip_util'}
}
# cut corner: 2013+ is all the same
for yr in range(2013, END_YEAR+1):
    readin_dict[yr] = readin_dict[END_YEAR].copy()
    readin_dict[yr]['files'] = [f'{yr}/1___Utility_Y{yr}.xlsx']

readin_dict[2012] = {
    'files': [f'{2012}/UtilityY{2012}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars': {'city':'city_util', 'state':'state_util', 'zip5':'zip_util'}
}
readin_dict[2011] = {
    'files': [f'{2011}/UtilityY{2011}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars': {'city':'city_util', 'state':'state_util', 'zip5':'zip_util'}
}
readin_dict[2010] = {
        'files': [f'{2010}/UtilityY{2010}.xls'],
        'excel_params':{'header':0, 'sheet_name':None},
        'rename_vars': {'utility_city':'city_util', 'utility_state':'state_util', 'utility_zip5':'zip_util'}
}
readin_dict[2009] = {
        'files': [f'{2009}/UtilityY{str(2009)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars': {'utility_city':'city_util', 'utility_state':'state_util', 'utility_zip5':'zip_util'}
}
readin_dict[2008] = {
        'files': [f'{2008}/UtilY{str(2008)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars': {'utilcode':'utility_id', 'utilname':'utility_name', 'city':'city_util', 'state':'state_util', 'zipcode':'zip_util'}
}
readin_dict[2007] = {
        'files': [f'{2007}/UtilY{str(2007)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars': {'utilcode':'utility_id', 'utilname':'utility_name', 'city':'city_util', 'state':'state_util', 'zipcode':'zip_util'}
}
readin_dict[2006] = {
        'files': [f'{2006}/UtilY{str(2006)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars': {'utilcode':'utility_id', 'utilname':'utility_name', 'city':'city_util', 'state':'state_util', 'zipcode':'zip_util'}
}

# %%
if __name__ == '__main__':
    # read-in parameters
    vars_keep = ['utility_id', 'utility_name', 'city_util', 'state_util', 'zip_util', 'entity_type']
    df = readin_eia_years(f'{PATH_RAW}eia/f860/', readin_dict, START_YEAR)
    # set datatypes
    df['utility_id'] = pd.to_numeric(df.utility_id).astype('Int64')
    df['zip_util'] = pd.to_numeric(df.zip_util.astype(str).str.strip(), errors='coerce').astype('Int64')
    # save intermediate file
    df.to_parquet(PATH_INTERIM + 'eia860_utility.parquet', index=False)

    # look at duplicates
    df['dup_key'] = df.utility_id.astype(str) + '_' + df.year.astype(str)
    df['duplicate'] = df.dup_key.isin(df.loc[df.dup_key.duplicated(), 'dup_key'].drop_duplicates())
    print('number of dups, by year: ')
    print(df.loc[df.duplicate].groupby('year')['utility_id'].count())
    # NOTE: All the issues are in 2010, and it looks like it's because of differences in naming etc.
    df['num_cols_nan'] = df.isna().sum(axis=1)
    df['min_num_cols_nan'] = df.groupby('dup_key')['num_cols_nan'].transform('min')
    # drop the duplicate row with the most missing info. If both have same missing, doesn't matter
    df['dup_keep'] = True
    df.loc[(df.num_cols_nan != df.min_num_cols_nan) & df.duplicate, 'dup_keep'] = False
    df_dedup = df.loc[df.dup_keep]
    # for remaining, drop duplicate arbitrarily
    df_dedup = df_dedup.loc[~df_dedup.dup_key.duplicated()]
    print('number of dups remaining:', df_dedup.dup_key.duplicated().sum())
    df_dedup.drop(columns=['dup_key', 'duplicate', 'num_cols_nan', 'min_num_cols_nan', 'dup_keep'], inplace=True)
    
    # save final file
    df_dedup.to_parquet(PATH_PROCESSED + 'eia860_utility.parquet', index=False)

    # summarize unique ids
    udf = df.rename(columns={'utility_id':'uid'})
    print('Utility dataset:')
    print(summarize_id_counts_byyear(udf.copy(), ['uid']))
