# Code leveraged from: https://github.com/USEPA/cam-api-examples/blob/main/Python/bulk_data_api_demo.py

# %%
import requests
import json
import sys
from datetime import datetime
from datetime import date
from tqdm import tqdm
import os

PATH_DATA = '../data/'

# Set API key
EPA_API_KEY = 'z6klBZ8QocrDVSam0LpGFfc5CQCgFaoawmH0eUjn'
epa_params = {'api_key': EPA_API_KEY}

# EPA DATA
# The bulk data api allows you to download prepackaged data sets. There are two endpoints for obtaining bulk data.
# The first is the /bulk-files endpoint which returns metadata about files. This metadata includes the path to the
# file.  The second is the /easey/bulk-files endpoint which along with the path, returns the actual file.
URL_BASE = 'https://api.epa.gov/easey/bulk-files/'
URL_METADATA = "https://api.epa.gov/easey/camd-services/bulk-files"



# %%
# GET BULK FILES
def get_bulk_files(parameters):
    # executing get request
    response = requests.get(URL_METADATA, params=parameters)

    # printing the response error message if the response is not successful
    print("Status code: "+str(response.status_code))
    if (int(response.status_code) > 399):
        sys.exit("Error message: "+response.json()['error']['message'])

    # converting the content from json format to a data frame
    resjson = response.content.decode('utf8').replace("'", '"')
    bulkFiles = json.loads(resjson)
    return bulkFiles

# get files
bulkFiles = get_bulk_files(epa_params)
# print out unique data types in the bulk data files
print('Unique data types in the bulk data files:')
print(set([fileObj['metadata']['dataType'] for fileObj in bulkFiles]))

# %%
# DOWNLOAD DATA
def download_files(filesToDownload, parameters):
    if len(filesToDownload) > 0:
        # loop through all files and download them
        for fileObj in tqdm(filesToDownload):
            url = URL_BASE+fileObj['s3Path']
            print('Full path to file on S3: '+url)
            # download and save file
            response = requests.get(url, params=parameters)
            # save file to disk in the data folder
            with open(PATH_DATA + 'epa/' + fileObj['filename'], 'wb') as f:
                f.write(response.content)
    else:
        print('No files to download')


# %%
# DOWNLOAD FACILITY DATA
# facility data
facFiles = [fileObj for fileObj in bulkFiles if 
            (fileObj['metadata']['dataType']=='Facility')]
facFilters = {'minYear': 2018}

filesToDownload = [fileObj for fileObj in facFiles if 
                   (int(fileObj['metadata']['year']) >= facFilters['minYear'])]
print('Number of files to download: ' + str(len(filesToDownload)))
download_files(filesToDownload, epa_params)



# %%
# DOWNLOAD EMISSIONS DATA
# facility data
emFiles = [fileObj for fileObj in bulkFiles if 
            (fileObj['metadata']['dataType']=='Emissions')]
print('Emission granularity:', 
      set(fileObj['metadata']['dataSubType'] for fileObj in emFiles))
emFiles = [fileObj for fileObj in emFiles if
           (fileObj['metadata']['dataSubType']=='Daily')]
emFilters = {'minYear': 2020}

filesToDownload = [fileObj for fileObj in emFiles if 
                   (int(fileObj['metadata']['year']) >= emFilters['minYear'])]
print('Number of files to download: '+str(len(filesToDownload)))
download_files(filesToDownload, epa_params)

