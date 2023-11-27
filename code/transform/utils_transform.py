import os
import pandas as pd
import re
from tqdm import tqdm

# set pandas outputs
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 70)
pd.set_option('mode.chained_assignment', None)


# read-in helper functions
def readin_eia_year(path_folder, path_file, excel_params, rename_vars=None, keep_vars=None):
    dfs = pd.read_excel(f'{path_folder}{path_file}', **excel_params)
    if type(dfs) != dict:
        dfs = {False:dfs}
    sdf = pd.DataFrame({})
    for k, df in dfs.items():
        df.columns = (df.columns.str.lower()
                        .str.rstrip()
                        .str.replace(' ', '_', regex=False)
                        .str.replace('(', '', regex=False).str.replace(')','', regex=False)
                        .str.replace('?', '', regex=False)
                        .str.replace('\n', '_', regex=False)
                        .str.replace('__', '_', regex=False)
                        .str.replace('.', '', regex=False)
                        .str.replace('&', 'and', regex=False))
        if rename_vars is not None:
            df = df.rename(columns=rename_vars)
        if keep_vars is not None:
            df = df.drop(columns=df.columns.difference(keep_vars))
        df['sheet'] = k.lower().replace(' ', '_')
        sdf = pd.concat([df, sdf], axis=0, ignore_index=True)
    return sdf

def readin_eia_years(path_folder, readin_dict, start_year, keep_vars=None):
    sdf = pd.DataFrame({})
    for y in readin_dict.keys():
        if y >= start_year:
            for f in tqdm(readin_dict[y]['files'], postfix=f'year: {y}'):
                df = readin_eia_year(path_folder, path_file=f, 
                                     excel_params=readin_dict[y]['excel_params'],
                                     rename_vars=readin_dict[y]['rename_vars'],
                                     keep_vars=keep_vars)
                df['year'] = y
                df['file'] = re.sub(r'Y(.*?)\.|xlsx|xls|\_|\.|\/|(\d{4}|\d{6})', '', f).lower()
                sdf = pd.concat([df, sdf.copy()], axis=0, ignore_index=True)
    return sdf


# READIN EPA DATA FROM FILE
def readin_epa(yr_start, yr_end, path_folder, vars_keep=None, vars_coll=None):
    df = pd.DataFrame()
    for yr in range(yr_start, yr_end+1):
        path_folder_yr = f'{path_folder}{yr}/'
        files = [f for f in os.listdir(path_folder_yr)]
        for f in tqdm(files):
            df_in = pd.read_csv(path_folder_yr + f, low_memory=False)
            # clean and subset columns
            df_in.columns = (df_in.columns.str.lower()
                            .str.replace(' ', '_', regex=False).str.replace('/', '_', regex=False)
                            .str.replace('&', 'and', regex=False)
                            .str.replace('(', '', regex=False).str.replace(')','', regex=False))
            if vars_keep is not None:
                df_in = df_in[vars_keep]
            if vars_coll is not None:
                df_in['quarter'] = pd.to_datetime(df_in.date).dt.to_period('Q')
                df_in['year'] = df_in.quarter.dt.year
                df_in = (df_in.groupby(vars_coll['id'])[vars_coll['val']]
                        .sum().reset_index())
            # only readin selected range
            # if df_in.year[0] not in range(yr_start, yr_end+1):
            #     continue        
            df = pd.concat([df_in, df], ignore_index=True, axis=0)
    return df


