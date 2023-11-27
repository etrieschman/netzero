# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

from utils_transform import readin_eia_years
from utils_summ import summarize_id_counts_byyear

# %%
# UTILITY DATA
readin_dict = {}
readin_dict[2021] = {
    'files': [f'{2021}/1___Utility_Y{2021}.xlsx'],
    'excel_params':{'header':1, 'sheet_name':None},
    'rename_vars':{'city':'city_util', 'state':'state_util', 'zip':'zip_util'}
}
# cut corner: 2013+ is all the same
for yr in range(2013, 2021+1):
    readin_dict[yr] = readin_dict[2021].copy()
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
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('transform_eia860_utility')

    year_start = snakemake.params.year_start
    year_end = snakemake.params.year_end
    path_raw = snakemake.params.indir
    intfile = snakemake.output.intfile
    outfile = snakemake.output.outfile

    # read-in parameters
    print('Reading in data...')
    vars_keep = ['utility_id', 'utility_name', 'city_util', 'state_util', 'zip_util', 'entity_type']
    df = readin_eia_years(path_raw, readin_dict, year_start)
    # set datatypes
    df['utility_id'] = pd.to_numeric(df.utility_id).astype('Int64')
    df['street_address'] = df.street_address.astype(str)
    df['zip_util'] = pd.to_numeric(df.zip_util.astype(str).str.strip(), errors='coerce').astype('Int64')
    # save intermediate file
    stringcols = df.select_dtypes(include='object').columns
    df[stringcols] = df[stringcols].astype(str)
    df.to_parquet(intfile, index=False)

    # drop variables
    newcols = ['year', 'sheet', 'file']
    df = df[vars_keep + newcols].copy()
    # look at duplicates
    print('Cleaning data...')
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
    print('Writing to file...')
    df_dedup.to_parquet(outfile, index=False)

    # summarize unique ids
    print('Summarizing unique identifiers...')
    udf = df.rename(columns={'utility_id':'uid'})
    print('Utility dataset:')
    print(summarize_id_counts_byyear(udf.copy(), ['uid']))
