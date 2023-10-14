import os
import pandas as pd
import re
from tqdm import tqdm


# set paths
PATH_DATA = '../data/'
PATH_RAW = PATH_DATA + 'raw/'
PATH_INTERIM = PATH_DATA + 'interim/'
PATH_PROCESSED = PATH_DATA + 'processed/'

for path in [PATH_DATA, PATH_INTERIM, PATH_PROCESSED]:
    if not os.path.exists(path):
        os.makedirs(path)

# define years
START_YEAR = 2018
END_YEAR = 2021


# set pandas outputs
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 70)
pd.set_option('mode.chained_assignment', None)


# read-in helper functions
def readin_eia_year(path_folder, path_file, excel_params, rename_vars=None):
    dfs = pd.read_excel(f'{path_folder}{path_file}', sheet_name=None, **excel_params)
    if type(dfs) != dict:
        dfs = {False:dfs}
    sdf = pd.DataFrame({})
    for k, df in dfs.items():
        df.columns = (df.columns.str.lower()
                        .str.replace(' ', '_')
                        .str.replace('(', '').str.replace(')','')
                        .str.replace('?', ''))
        if rename_vars is not None:
            df = df.rename(columns=rename_vars)
            df['sheet'] = k
        sdf = pd.concat([df, sdf], axis=0, ignore_index=True)
    return sdf

def readin_eia_years(path_folder, readin_dict, start_year):
    sdf = pd.DataFrame({})
    for y in readin_dict.keys():
        if y >= start_year:
            for f in tqdm(readin_dict[y]['files'], postfix=f'year: {y}'):
                df = readin_eia_year(path_folder, path_file=f, 
                                     excel_params=readin_dict[y]['excel_params'],
                                     rename_vars=readin_dict[y]['rename_vars'])
                df['year'] = y
                df['file'] = re.sub(r'Y(.*?)\.|xlsx|xls|\_|\d+|\/', '', f).lower()
                sdf = pd.concat([df, sdf.copy()], axis=0, ignore_index=True)
    return sdf


# READIN EPA DATA FROM FILE
def readin_epa(yr_start, yr_end, path_folder, vars_keep=None, vars_coll=None):
    files = [f for f in os.listdir(path_folder)]
    df = pd.DataFrame({})
    for f in tqdm(files):
        df_in = pd.read_csv(path_folder + f, low_memory=False)
        # clean and subset columns
        df_in.columns = (df_in.columns.str.lower()
                         .str.replace(' ', '_').str.replace('/', '_')
                         .str.replace('&', 'and')
                         .str.replace('(', '').str.replace(')',''))
        if vars_keep is not None:
            df_in = df_in[vars_keep]
        if vars_coll is not None:
            df_in['quarter'] = pd.to_datetime(df_in.date).dt.to_period('Q')
            df_in = (df_in.groupby(vars_coll['id']+['quarter'])[vars_coll['val']]
                     .sum().reset_index())
            df_in['year'] = df_in.quarter.dt.year
        # only readin selected range
        if df_in.year[0] not in range(yr_start, yr_end+1):
            continue        
        df = pd.concat([df_in, df], ignore_index=True, axis=0)
    return df


