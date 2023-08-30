import os
import pandas as pd
from tqdm import tqdm

PATH_DATA = '../data/'
PATH_EPA = PATH_DATA + 'raw/epa/'
PATH_EIA = PATH_DATA + 'raw/eia/'
PATH_PROCESSED = PATH_DATA + 'processed/'

def readin_eia(yr_start, yr_end, path_folder, path_file, vars_keep=None, header=0, suff='.xlsx'):
    years = range(yr_start, yr_end+1)
    sdf = pd.DataFrame({})
    for y in tqdm(years):
        df = pd.read_excel(f'{path_folder}/{y}/{path_file}{y}{suff}', 
                           header=header)
        df.columns = (df.columns.str.lower()
                      .str.replace(' ', '_')
                      .str.replace('(', '').str.replace(')',''))
        if vars_keep is not None:
            df = df[vars_keep]
        df['year'] = y
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