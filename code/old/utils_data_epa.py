import sys, json, io
import requests, zipfile, urllib
from tqdm import tqdm
import pandas as pd

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


# READIN EPA DATA FROM FILE
def readin_epa(yr_start, yr_end, path_folder, vars_keep=None, vars_coll=None):
    files = [f for f in os.listdir(path_folder)]
    df = pd.DataFrame({})
    for f in tqdm(files):
        df_in = pd.read_csv(path_folder + f, low_memory=False)
        # clean and subset columns
        df_in.columns = (df_in.columns.str.lower()
                         .str.replace(' ', '_').str.replace('/', '_')
                         .str.replace('&', 'and')
                         .str.replace('(', '').str.replace(')',''))
        if vars_keep is not None:
            df_in = df_in[vars_keep]
        if vars_coll is not None:
            df_in['quarter'] = pd.to_datetime(df_in.date).dt.to_period('Q')
            df_in = (df_in.groupby(vars_coll['id']+['quarter'])[vars_coll['val']]
                     .sum().reset_index())
            df_in['year'] = df_in.quarter.dt.year
        # only readin selected range
        if df_in.year[0] not in range(yr_start, yr_end+1):
            continue        
        df = pd.concat([df_in, df], ignore_index=True, axis=0)
    return df

