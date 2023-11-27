# %%
import numpy as np
import pandas as pd
from tqdm import tqdm
import requests, zipfile, urllib
from io import BytesIO
from bs4 import BeautifulSoup

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
    path = f'{path_save}' if downfiles[0].startswith(f'{year}/') else f'{path_save}/{year}/'
    zip.extractall(path, members=downfiles)

# DOWNLOAD A DIRECT FILE
def download_eia_file(url, path_save):
    fn = url.split('/')[-1]
    req = requests.get(url)
    with open(path_save + fn, 'wb') as output:
        output.write(req.content)


# %%
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.realpath(__file__))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('extract_eia')

    year_start = snakemake.params.year_start
    year_end = snakemake.params.year_end
    path_eia_f860 = snakemake.output.eia_f860
    path_eia_f861 = snakemake.output.eia_f861
    path_eia_f923 = snakemake.output.eia_f923

    # DOWNLOAD UTILITY-LEVEL DATA
    print('Downloading EIA 861 data')
    url = 'https://www.eia.gov/electricity/data/eia861/'
    zip_paths = get_eia_files(url, year_start)
    for zp in tqdm(zip_paths):
        download_eia_zip(zp, path_eia_f861)

    # DOWNLOAD PLANT-LEVEL DATA
    print('Downloading EIA 860 data')
    url = 'https://www.eia.gov/electricity/data/eia860/'
    zip_paths = get_eia_files(url, year_start)
    for zp in tqdm(zip_paths):
        download_eia_zip(zp, path_eia_f860)

    # DOWNLOAD OPERATIONS-LEVEL DATA
    print('Downloading EIA 923 data')
    url = 'https://www.eia.gov/electricity/data/eia923/'
    zip_paths = get_eia_files(url, year_start)
    for zp in tqdm(zip_paths):
        download_eia_zip(zp, path_eia_f923)

    # DOWNLOAD EMISSIONS DATA
    # print('Downloading EIA estimated emissions data')
    # url = 'https://www.eia.gov/electricity/data/emissions/'
    # file_paths = get_eia_files(url, year_start, suffs='.xlsx')
    # for f in tqdm(file_paths):
    #     if 'region' in f:
    #         continue
    #     download_eia_file(f, PATH_EIA + 'emissions/')

