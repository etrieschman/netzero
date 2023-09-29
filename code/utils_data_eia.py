import requests, zipfile, urllib
from io import BytesIO
from bs4 import BeautifulSoup
import pandas as pd
from tqdm import tqdm
import re

# GET ALL FILES THAT END IN SPECIFIC SUFFIX
def get_eia_files(url, start_year, suffs='.zip'):
    req = requests.get(url)
    soup = BeautifulSoup(req.content, 'html.parser')
    zips = [z['href'] for z in soup.findAll('a', href=True) if z['href'].endswith(suffs)]
    zips_sub = [z for z in zips if z[z.find(suffs)-4:z.find(suffs)-2] in ['19', '20']]
    return [f'{url}{z}' for z in zips_sub if int(z[z.find(suffs)-4:z.find(suffs)]) >= start_year]


# DOWNLOAD A ZIP FILE
def download_eia_zip(url, path_save):
    fn = url.split('/')[-1]
    year = int(fn[fn.find('.zip') - 4: fn.find('.zip')])
    # Downloading the file by sending the request to the URL
    req = requests.get(url)
    # extracting the zip file contents
    zip = zipfile.ZipFile(BytesIO(req.content))
    files = zip.namelist()
    downfiles = [f for f in files if f.endswith(('.xls','.xlsx'))]
    path = f'{path_save}' if downfiles[0].startswith(f'{year}/') else f'{path_save}{year}/'
    zip.extractall(path, members=downfiles)


# DOWNLOAD A DIRECT FILE
def download_eia_file(url, path_save):
    fn = url.split('/')[-1]
    req = requests.get(url)
    with open(path_save + fn, 'wb') as output:
        output.write(req.content)


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
                df['gen_category'] = re.sub(r'Y(.*?)\.|xlsx|xls|\_|\d+|\/', '', f).lower()
                sdf = pd.concat([df, sdf], axis=0, ignore_index=True)

    return sdf

