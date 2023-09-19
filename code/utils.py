import os
import pandas as pd
from tqdm import tqdm
import re

PATH_RESULTS = '../results/'
PATH_DATA = '../data/'
PATH_EPA = PATH_DATA + 'raw/epa/'
PATH_EIA = PATH_DATA + 'raw/eia/'
PATH_PROCESSED = PATH_DATA + 'processed/'

def readin_eia(path_folder, readin_dict):
    sdf = pd.DataFrame({})
    for y in tqdm(readin_dict.keys()):
        df = pd.read_excel(f'{path_folder}{readin_dict[y]["path_file"]}', 
                           **readin_dict[y]['excel_params'])
        df.columns = (df.columns.str.lower()
                      .str.replace(' ', '_')
                      .str.replace('(', '').str.replace(')',''))
        # print(y)
        # print(df.columns)
        df = df.rename(columns=readin_dict[y]['rename_vars'])
        df = df[df.columns.intersection(readin_dict[y]['vars_keep'])]
        df['year'] = y
        sdf = pd.concat([df, sdf], axis=0, ignore_index=True)
    return sdf

def readin_eia_gen(path_folder, readin_dict):
    sdf = pd.DataFrame({})
    for y in readin_dict.keys():
        for f in tqdm(readin_dict[y]['files'], postfix=f'Year: {y}'):
            dfs = pd.read_excel(f'{path_folder}{f}', sheet_name=None,
                                    **readin_dict[y]['excel_params'])
            if type(dfs) != dict:
                dfs = {False:dfs}
            for k, df in dfs.items():
                df.columns = (df.columns.str.lower()
                            .str.replace(' ', '_')
                            .str.replace('(', '').str.replace(')',''))
                if 'rename_vars' in readin_dict[y].keys():
                    df = df.rename(columns=readin_dict[y]['rename_vars'])
                df = df[df.columns.intersection(readin_dict[y]['vars_keep'])]
                df['year'] = y
                df['sheet'] = k
                df['gen_category'] = re.sub(r'Y(.*?)\.|xlsx|xls|\_|\d+|\\', '', f)
                sdf = pd.concat([df, sdf], axis=0, ignore_index=True)

    return sdf


def readin_epa(yr_start, yr_end, key, vars_keep=None, vars_coll=None):
    files = [f for f in os.listdir(PATH_EPA) if key in f.lower()]
    df = pd.DataFrame({})
    for f in tqdm(files):
        df_in = pd.read_csv(PATH_EPA + f, low_memory=False)
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