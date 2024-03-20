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

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)

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
    'rename_vars':{'data_type_o_=_observed_i_=_imputed':'data_type',
                   'thousand_dollars':'rev_1k_res', 'megawatthours':'sales_mwh_res', 'count':'cust_res',
                   'thousand_dollars1':'rev_1k_com', 'megawatthours1':'sales_mwh_com', 'count1':'cust_com',
                   'thousand_dollars2':'rev_1k_ind', 'megawatthours2':'sales_mwh_ind', 'count2':'cust_ind',
                   'thousand_dollars3':'rev_1k_trans', 'megawatthours3':'sales_mwh_trans', 'count3':'cust_trans',
                   'thousand_dollars4':'rev_1k_tot', 'megawatthours4':'sales_mwh_tot', 'count4':'cust_tot'}
}
for yr in range(2015, 2022):
    readin_dict[yr] = readin_dict[2022].copy()
    readin_dict[yr]['files'] = [f'{yr}/Sales_Ult_Cust_{yr}.xlsx']
for yr in range(2013, 2015):
    readin_dict[yr] = readin_dict[2022].copy()
    readin_dict[yr]['files'] = [f'{yr}/Sales_Ult_Cust_{yr}.xls']
for yr in range(2001, 2013):
    readin_dict[yr] = readin_dict[2022].copy()
    readin_dict[yr]['files'] = [f'{yr}/Sales_Ult_Cust_{yr}.xlsx']
readin_dict[2012]['rename_vars'] = {
   'data_type_o_=_observed_i_=_imputed':'data_type',
    'thousands_dollars':'rev_1k_res', 'megawatthours':'sales_mwh_res', 'count':'cust_res',
    'thousands_dollars1':'rev_1k_com', 'megawatthours1':'sales_mwh_com', 'count1':'cust_com',
    'thousands_dollars2':'rev_1k_ind', 'megawatthours2':'sales_mwh_ind', 'count2':'cust_ind',
    'thousands_dollars3':'rev_1k_oth', 'megawatthours3':'sales_mwh_oth', 'count3':'cust_oth',
    'thousands_dollars4':'rev_1k_tot', 'megawatthours4':'sales_mwh_tot', 'count4':'cust_tot'}
for yr in range(1990, 2001):
    readin_dict[yr] = readin_dict[2022].copy()
    readin_dict[yr]['files'] = [f'{yr}/Sales_Ult_Cust_{yr}.xlsx']
    readin_dict[yr]['rename_vars'] = {
        'data_type_o_=_observed_i_=_imputed':'data_type',
        'thousand_dollars':'rev_1k_res', 'megawatthours':'sales_mwh_res', 'count':'cust_res',
        'thousand_dollars1':'rev_1k_com', 'megawatthours1':'sales_mwh_com', 'count1':'cust_com',
        'thousand_dollars2':'rev_1k_ind', 'megawatthours2':'sales_mwh_ind', 'count2':'cust_ind',
        'thousand_dollars3':'rev_1k_oth', 'megawatthours3':'sales_mwh_oth', 'count3':'cust_oth',
        'thousand_dollars4':'rev_1k_tot', 'megawatthours4':'sales_mwh_tot', 'count4':'cust_tot'}


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
    vars = [col for col in udf.columns if col.startswith(('rev_1k', 'sales_mwh', 'cust'))]
    for var in vars:
        udf[var] = pd.to_numeric(udf[var], errors='coerce').fillna(0)
    # print missing by year
    print(udf.groupby('year').agg(lambda x: x.isna().sum()))

    # add a filter
    udf['filter_out'] = False
    mask = udf.utility_number.isna() & (udf.year != 2002)
    udf.loc[mask, 'filter_out'] = True
    mask = (udf.utility_number == 88888) | (udf.utility_number == 99999)
    udf.loc[mask, 'filter_out'] = True

    # save to file
    udf.to_csv(PATH_RESULTS + 'df_stacked_f860_sales_ult_cust.csv', index=False)
    
# %%
# NOTE FROM EXCEL
# To calculate a state or the US total, sum Parts (A,B,C & D) for Revenue, but only Parts (A,B & D) for Sales and Customers.\nTo avoid double counting of customers, the aggregated customer counts for the states and US do not include the customer count for respondents with ownership code 'Behind the Meter'.\nThis group consists of Third Party Owners of rooftop solar systems.
