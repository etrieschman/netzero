# %%
import numpy as np
import pandas as pd
import requests, zipfile, urllib
from io import BytesIO
from bs4 import BeautifulSoup
from tqdm import tqdm

from utils import PATH_EIA, PATH_PROCESSED

# %%
# DOWNLOAD A ZIP FILE
def download_zip(url, path_save):
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


# %%
# GET ALL ZIP FILES
def get_eia_zips(url, start_year):
    req = requests.get(url)
    soup = BeautifulSoup(req.content, 'html.parser')
    zips = [z['href'] for z in soup.findAll('a', href=True) if z['href'].endswith('.zip')]
    zips_sub = [z for z in zips if z[z.find('.zip')-4:z.find('.zip')-2] in ['19', '20']]
    return [f'{url}{z}' for z in zips_sub if int(z[z.find('.zip')-4:z.find('.zip')]) >= start_year]


# %%
START_YEAR = 2018

# %%
# DOWNLOAD UTILITY-LEVEL DATA
url = 'https://www.eia.gov/electricity/data/eia861/'
zip_paths = get_eia_zips(url, START_YEAR)
for zp in tqdm(zip_paths):
    download_zip(zp, PATH_EIA + 'f861/')

# %%
# DOWNLOAD PLANT-LEVEL DATA
url = 'https://www.eia.gov/electricity/data/eia860/'
zip_paths = get_eia_zips(url, START_YEAR)
for zp in tqdm(zip_paths):
    download_zip(zp, PATH_EIA + 'f860/')


# %%
# DOWNLOAD OPERATIONS-LEVEL DATA
url = 'https://www.eia.gov/electricity/data/eia923/'
zip_paths = get_eia_zips(url, START_YEAR)
for zp in tqdm(zip_paths):
    download_zip(zp, PATH_EIA + 'f923/')
# %%
