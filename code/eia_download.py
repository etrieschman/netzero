# %%
import numpy as np
import pandas as pd
from tqdm import tqdm

from utils_data_eia import get_eia_files, download_eia_zip, download_eia_file
from utils import PATH_EIA, START_YEAR


# %%
# DOWNLOAD UTILITY-LEVEL DATA
url = 'https://www.eia.gov/electricity/data/eia861/'
zip_paths = get_eia_files(url, START_YEAR)
for zp in tqdm(zip_paths):
    download_eia_zip(zp, PATH_EIA + 'f861/')

# %%
# DOWNLOAD PLANT-LEVEL DATA
url = 'https://www.eia.gov/electricity/data/eia860/'
zip_paths = get_eia_files(url, START_YEAR)
for zp in tqdm(zip_paths):
    download_eia_zip(zp, PATH_EIA + 'f860/')


# %%
# DOWNLOAD OPERATIONS-LEVEL DATA
url = 'https://www.eia.gov/electricity/data/eia923/'
zip_paths = get_eia_files(url, START_YEAR)
for zp in tqdm(zip_paths):
    download_eia_zip(zp, PATH_EIA + 'f923/')

# %%
# DOWNLOAD EMISSIONS DATA
url = 'https://www.eia.gov/electricity/data/emissions/'
file_paths = get_eia_files(url, START_YEAR, suffs='.xlsx')
for f in tqdm(file_paths):
    if 'region' in f:
        continue
    download_eia_file(f, PATH_EIA + 'emissions/')
# %%
