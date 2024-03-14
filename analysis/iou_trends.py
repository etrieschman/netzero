# %%
import numpy as np
import pandas as pd
import re
from tqdm import tqdm
import requests, zipfile, urllib
from io import BytesIO
from bs4 import BeautifulSoup
from pathlib import Path

PATH_DATA = '../data/'
PATH_EIA_F861 = PATH_DATA + 'raw/eia/f861/'
Path(PATH_EIA_F861).mkdir(exist_ok=True, parents=True)
PATH_RESULTS = '../results/analysis/iou_trends/'
Path(PATH_RESULTS).mkdir(exist_ok=True, parents=True)
YEAR_START = 1990

# %%
# SCRAPE ALL FILEPATHS THAT END IN SPECIFIC SUFFIX
def get_eia_files(url, start_year, suffs='.zip'):
    req = requests.get(url)
    soup = BeautifulSoup(req.content, 'html.parser')
    zips = [z['href'] for z in soup.findAll('a', href=True) if z['href'].endswith(suffs)]
    zips_sub = [z for z in zips if z[z.find(suffs)-4:z.find(suffs)-2] in ['19', '20']]
    return [f'{url}{z}' for z in zips_sub if int(z[z.find(suffs)-4:z.find(suffs)]) >= start_year]


# DOWNLOAD AND EXTRACT A ZIP FILE FROM FILEPATH
def download_eia_zip(url, path_save):
    fn = url.split('/')[-1]
    year = int(fn[fn.find('.zip') - 4: fn.find('.zip')])
    # Downloading the file by sending the request to the URL
    req = requests.get(url)
    # extracting the zip file contents
    zip = zipfile.ZipFile(BytesIO(req.content))
    files = zip.namelist()
    downfiles = [f for f in files if f.lower().endswith(('.xls','.xlsx'))]
    path = path_save if downfiles[0].startswith(f'{year}/') else path_save / f'{year}'
    zip.extractall(path, members=downfiles)

# read-in helper functions
def readin_eia_year(path_folder, path_file, excel_params, rename_vars=None, keep_vars=None):
    dfs = pd.read_excel(path_folder + path_file, **excel_params)
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
        # df['sheet'] = k.lower().replace(' ', '_')
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

readin_dict = {}
readin_dict[2022] = {
    'files': [f'{2022}/Sales_Ult_Cust_{2022}.xlsx'],
    'excel_params':{'header':2, 'sheet_name':'States'},
    'rename_vars':{'data_type_o_=_observed_i_=_imputed':'data_type'}
}
for yr in range(2015, 2022):
    readin_dict[yr] = readin_dict[2022].copy()
    readin_dict[yr]['files'] = [f'{yr}/Sales_Ult_Cust_{yr}.xlsx']


# %%
if __name__ == '__main__':
    # DOWNLOAD UTILITY-LEVEL DATA
    print('Downloading EIA 861 data')
    url = 'https://www.eia.gov/electricity/data/eia861/'
    zip_paths = get_eia_files(url, YEAR_START)
    for zp in tqdm(zip_paths):
        download_eia_zip(zp, PATH_EIA_F861)

    # READIN AND STACK DATA
    udf = readin_eia_years(PATH_EIA_F861, readin_dict, YEAR_START)

    # SUMMARIZE
    # NOTE: this fails because we have conflicting dtypes
    udf['sales_sector_total'] = udf.groupby(['data_year'])['megawatthours4'].transform('sum')
    udf['cust_sector_total'] = udf.groupby(['data_year'])['count4'].transform('sum')
    udf['pct_sales_sector_total'] = udf.megawatthours4 / udf.sales_sector_total
    udf['pct_cust_sector_total'] = udf.count4 / udf.cust_sector_total
# %%
