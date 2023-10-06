# %%
import subprocess
import os
import pandas as pd

# %%

# %%
# Step 1: Setup datastore
os.getcwd()
datastore_folder = 'data_pudl'
datastore_path = f'../../{datastore_folder}'
if not os.path.exists(datastore_path):
   os.makedirs(datastore_path)
subprocess.run(["pudl_setup", datastore_path])
# Step 2: pull data
subprocess.run(['pudl_datastore', '--dataset', 'eia860'])
subprocess.run(['pudl_datastore', '--dataset', 'eia923'])
subprocess.run(['pudl_datastore', '--dataset', 'epacems'])
# Step 3. Extract Transform Load (ETL) data
settings = 'etl_custom.yml'
subprocess.run(['pudl_etl', f'{datastore_path}/{settings}'])