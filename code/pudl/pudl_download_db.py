# %%
# SETUP
import subprocess, os
import tarfile 


# PICK PUDL DATABASE VERSION
# ACCESSED ON: 2023.10.11
# ACCESSED AT: https://zenodo.org/record/7472137
ZENODO_DOI = '10.5281/zenodo.7472137'
# SETUP PATH TO SAVE DATA
PATH_DATA = '../../data/pudl/'
if not os.path.exists(PATH_DATA):
   os.makedirs(PATH_DATA)

# %%
# IMPORT PUDL DATABASE
# NOTE: THIS IS A 6GB FILE THAT TAKES >15 MIN TO DOWNLOAD
subprocess.run(['zenodo_get', 
                f'--doi={ZENODO_DOI}',
                f'--output-dir={PATH_DATA}'])

# %%
# UNZIP INTO FINAL LOCATION
filenames = [fn for fn in os.listdir(PATH_DATA) if fn.endswith('.tgz')]
for fn in filenames:
   file = tarfile.open(PATH_DATA + fn) 
   print(file.getnames()) 
   file.extractall(PATH_DATA) 
   file.close() 


# %%
# ===== OLD =====
# PREVIOUS ATTEMPT WHERE I DOWNLOAD INDIVIDUAL FILES AND TRY TO RUN THE ETL PIPELINE
# I ABANDONED THIS APPROACH BECAUSE A LOT OF THE FUNCTIONS I WAS USING WERE DEPRECATED
# PUDL SEEMS TO BE MOVING TOWARD ALWAYS USING ALL OF THE DATA, WHICH REQUIRES MORE DATA STORAGE

# # Step 1: Setup datastore
# subprocess.run(["pudl_setup", datastore_path])

# # Step 2: pull data
# subprocess.run(['pudl_datastore', '--dataset', 'eia860'])
# subprocess.run(['pudl_datastore', '--dataset', 'eia923'])
# subprocess.run(['pudl_datastore', '--dataset', 'epacems'])

# # Step 3. Extract Transform Load (ETL) data
# settings = 'etl_custom.yml'
# subprocess.run(['pudl_etl', f'{datastore_path}/{settings}'])