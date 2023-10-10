# %%
import subprocess
import numpy as np
import os
import pandas as pd
import pudl.workspace.datastore

os.getcwd()
datastore_folder = 'data_pudl'
datastore_path = f'../../{datastore_folder}'
if not os.path.exists(datastore_path):
   os.makedirs(datastore_path)

# %%
# Step 1: Setup datastore
subprocess.run(["pudl_setup", datastore_path])

# %%
# Step 2: pull data
subprocess.run(['pudl_datastore', '--dataset', 'eia860'])
subprocess.run(['pudl_datastore', '--dataset', 'eia923'])
subprocess.run(['pudl_datastore', '--dataset', 'epacems'])
# year_range = np.arange(2020, 2022)
# for year in year_range:
#    print(year)
#    subprocess.run(['pudl_datastore', '--dataset', 'eia860', '--loglevel', 'INFO', '--partition', f'year={year}'])
#    subprocess.run(['pudl_datastore', '--dataset', 'eia923', '--loglevel', 'INFO', '--partition', f'year={year}'])
#    subprocess.run(['pudl_datastore', '--dataset', 'epacems', '--loglevel', 'INFO', '--partition', f'year={year}'])

# %%
# Step 3. Extract Transform Load (ETL) data
settings = 'etl_custom.yml'
subprocess.run(['pudl_etl', f'{datastore_path}/{settings}'])