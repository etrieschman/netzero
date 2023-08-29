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