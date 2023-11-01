# %%
# Code leveraged from: https://github.com/USEPA/cam-api-examples/blob/main/Python/bulk_data_api_demo.py
import requests
import json
import sys
import pandas as pd
from datetime import datetime
from datetime import date
from tqdm import tqdm
import os, io

from credentials import EPA_API_KEY # local user: save credentials.py file with api key

# Set API key
START_YEAR = 2015
END_YEAR = 2021
PATH_EPA = '../../data/raw/epa/'
epa_params = {'api_key': EPA_API_KEY}

# GET BULK FILES
def get_bulk_epa_files(url_metadata, parameters):
    # executing get request
    response = requests.get(url_metadata, params=parameters)

    # printing the response error message if the response is not successful
    print("Status code: "+str(response.status_code))
    if (int(response.status_code) > 399):
        sys.exit("Error message: "+response.json()['error']['message'])

    # converting the content from json format to a data frame
    resjson = response.content.decode('utf8').replace("'", '"')
    bulkFiles = json.loads(resjson)
    return bulkFiles

# DOWNLOAD DATA
def download_epa_files(url_base, filesToDownload, filename_prefix, parameters, subset_col=None):
    if len(filesToDownload) > 0:
        # loop through all files and download them
        for fileObj in tqdm(filesToDownload):
            url = url_base+fileObj['s3Path']
            # print('Full path to file on S3: '+url)
            # download and save file
            response = requests.get(url, params=parameters)
            # subset data
            # converters = {col:'Float64' for col in numeric_cols}
            content = pd.read_csv(io.StringIO(response.text), dtype=str)
            if subset_col is not None:
                content = content.astype({subset_col:'Float64'})
                # print(fileObj['s3Path'])
                # print(content.shape)
                content = content.loc[content[subset_col] != 0.]
                # print(content.shape)
                # content = content.astype({subset_col:'Float64'}).loc[content[subset_col] != 0.]
            # save file to disk in the data folder
            content.to_csv(filename_prefix + fileObj['filename'])
            # with open(PATH_EPA + fileObj['filename'], 'wb') as f:
            #     f.write(content)
    else:
        print('No files to download')


# %%
if __name__ == '__main__':
      # The bulk data api allows you to download prepackaged data sets. There are two endpoints for obtaining bulk data.
      # The first is the /bulk-files endpoint which returns metadata about files. This metadata includes the path to the
      # file.  The second is the /easey/bulk-files endpoint which along with the path, returns the actual file.
      URL_BASE = 'https://api.epa.gov/easey/bulk-files/'
      URL_METADATA = "https://api.epa.gov/easey/camd-services/bulk-files"
      # get files
      print('Getting filepaths...')
      bulkFiles = get_bulk_epa_files(URL_METADATA, epa_params)
      # print out unique data types in the bulk data files
      print('Unique data types in the bulk data files:')
      print(set([fileObj['metadata']['dataType'] for fileObj in bulkFiles]))

      # DOWNLOAD FACILITY DATA
      print('Downloading facility data')
      facFiles = [fileObj for fileObj in bulkFiles if 
                  (fileObj['metadata']['dataType']=='Facility')]
      facFilters = {'minYear': START_YEAR}

      filesToDownload = [fileObj for fileObj in facFiles if 
                        ((int(fileObj['metadata']['year']) >= START_YEAR) & 
                        (int(fileObj['metadata']['year']) <= END_YEAR))]
      print('Number of files to download: ' + str(len(filesToDownload)))
      download_epa_files(URL_BASE, filesToDownload, PATH_EPA + '/facility/', epa_params)

      # DOWNLOAD DAILY EMISSIONS DATA
      emFiles = [fileObj for fileObj in bulkFiles if 
                  (fileObj['metadata']['dataType']=='Emissions')]
      print('Emission granularity:', 
            set(fileObj['metadata']['dataSubType'] for fileObj in emFiles))
      print('Downloading daily emissions data...')
      emFiles = [fileObj for fileObj in emFiles if
            (fileObj['metadata']['dataSubType']=='Daily')]
      filesToDownload = [fileObj for fileObj in emFiles if 
                        (('stateCode' in fileObj['metadata']) & 
                        (int(fileObj['metadata']['year']) >= START_YEAR) & 
                        (int(fileObj['metadata']['year']) <= END_YEAR))]
      print('Number of files to download: '+ str(len(filesToDownload)))
      download_epa_files(URL_BASE, filesToDownload, PATH_EPA + '/emissions/daily/', epa_params,
                  subset_col='Operating Time Count')

      # DOWNLOAD HOURLY EMISSIONS DATA
    #   START_YEAR, END_YEAR = 2013, 2021
    #   emFiles = [fileObj for fileObj in bulkFiles if 
    #               (fileObj['metadata']['dataType']=='Emissions')]
    #   print('Emission granularity:', 
    #         set(fileObj['metadata']['dataSubType'] for fileObj in emFiles))
    #   print('Downloading hourly emissions data...')  
    #   emFiles = [fileObj for fileObj in emFiles if
    #         (fileObj['metadata']['dataSubType']=='Hourly')]
    #   filesToDownload = [fileObj for fileObj in emFiles if 
    #                     (('stateCode' in fileObj['metadata']) & 
    #                     (int(fileObj['metadata']['year']) >= START_YEAR) & 
    #                     (int(fileObj['metadata']['year']) <= END_YEAR))]
    #   print('Number of files to download: '+ str(len(filesToDownload)))
    #   download_epa_files(URL_BASE, filesToDownload, PATH_EPA + '/emissions/hourly/', epa_params,
    #               subset_col='Operating Time')

# %%
