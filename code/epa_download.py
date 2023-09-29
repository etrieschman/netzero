# Code leveraged from: https://github.com/USEPA/cam-api-examples/blob/main/Python/bulk_data_api_demo.py

# %%
import requests
import json
import sys
import pandas as pd
from datetime import datetime
from datetime import date
from tqdm import tqdm
import os, io

from credentials import EPA_API_KEY # local user: save credentials.py file with api key

from utils import PATH_EPA, START_YEAR
from utils_data_epa import get_epa_bulk_files, download_epa_files

# Set API key

epa_params = {'api_key': EPA_API_KEY}

# %%
# EPA DATA
# The bulk data api allows you to download prepackaged data sets. There are two endpoints for obtaining bulk data.
# The first is the /bulk-files endpoint which returns metadata about files. This metadata includes the path to the
# file.  The second is the /easey/bulk-files endpoint which along with the path, returns the actual file.
URL_BASE = 'https://api.epa.gov/easey/bulk-files/'
URL_METADATA = "https://api.epa.gov/easey/camd-services/bulk-files"
# get files
bulkFiles = get_epa_bulk_files(URL_METADATA, epa_params)
# print out unique data types in the bulk data files
print('Unique data types in the bulk data files:')
print(set([fileObj['metadata']['dataType'] for fileObj in bulkFiles]))


# %%
# DOWNLOAD FACILITY DATA
# facility data
facFiles = [fileObj for fileObj in bulkFiles if 
            (fileObj['metadata']['dataType']=='Facility')]
facFilters = {'minYear': START_YEAR}

filesToDownload = [fileObj for fileObj in facFiles if 
                   ((int(fileObj['metadata']['year']) >= START_YEAR) & 
                   (int(fileObj['metadata']['year']) <= END_YEAR))]
print('Number of files to download: ' + str(len(filesToDownload)))
download_epa_files(URL_BASE, filesToDownload, PATH_EPA + '/facility/', epa_params)



# %%
# DOWNLOAD DAILY EMISSIONS DATA
emFiles = [fileObj for fileObj in bulkFiles if 
            (fileObj['metadata']['dataType']=='Emissions')]
print('Emission granularity:', 
      set(fileObj['metadata']['dataSubType'] for fileObj in emFiles))
emFiles = [fileObj for fileObj in emFiles if
           (fileObj['metadata']['dataSubType']=='Daily')]
filesToDownload = [fileObj for fileObj in emFiles if 
                   (('stateCode' in fileObj['metadata']) & 
                    (int(fileObj['metadata']['year']) >= START_YEAR) & 
                   (int(fileObj['metadata']['year']) <= END_YEAR))]
print('Number of files to download: '+ str(len(filesToDownload)))
download_epa_files(URL_BASE, filesToDownload, PATH_EPA + '/emissions/daily/', epa_params,
               subset_col='Operating Time Count')


# %%
# DOWNLOAD HOURLY EMISSIONS DATA
START_YEAR, END_YEAR = 2013, 2021
emFiles = [fileObj for fileObj in bulkFiles if 
            (fileObj['metadata']['dataType']=='Emissions')]
print('Emission granularity:', 
      set(fileObj['metadata']['dataSubType'] for fileObj in emFiles))
emFiles = [fileObj for fileObj in emFiles if
           (fileObj['metadata']['dataSubType']=='Hourly')]
filesToDownload = [fileObj for fileObj in emFiles if 
                   (('stateCode' in fileObj['metadata']) & 
                    (int(fileObj['metadata']['year']) >= START_YEAR) & 
                   (int(fileObj['metadata']['year']) <= END_YEAR))]
print('Number of files to download: '+ str(len(filesToDownload)))
download_epa_files(URL_BASE, filesToDownload, PATH_EPA + '/emissions/hourly/', epa_params,
               subset_col='Operating Time')
