# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

from utils_transform import (
    PATH_RAW, PATH_INTERIM, PATH_PROCESSED, START_YEAR, END_YEAR)
from utils_transform import readin_eia_years
from utils_summ import summarize_id_counts_byyear

# %%
# OWNERSHIP DATA
readin_dict = {}
readin_dict[END_YEAR] = {
    'files': [f'{END_YEAR}/4___Owner_Y{END_YEAR}.xlsx'],
    'excel_params':{'header':1},
    'rename_vars':{'city_owner':'owner_city', 'owner_state':'state_owner', 'owner_zip':'zip_owner'}
}
# cut corner: 2013+ is all the same
for yr in range(2013, END_YEAR+1):
    readin_dict[yr] = readin_dict[END_YEAR].copy()
    readin_dict[yr]['files'] = [f'{yr}/4___Owner_Y{yr}.xlsx']

readin_dict[2012] = {
    'files':    [f'{2012}/OwnerY{2012}.xlsx'],
    'excel_params': {'header':1},
    'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2011] = {
    'files':    [f'{2011}/OwnershipY{2011}.xlsx'],
    'excel_params': {'header':1},
    'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2010] = {
        'files':    [f'{2010}/OwnerY{2010}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2009] = {
        'files':    [f'{2009}/OwnerY{str(2009)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2008] = {
        'files':    [f'{2008}/OwnerY{str(2008)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plantcode':'plant_code', 'gencode':'generator_id', 'owner_id':'ownership_id'}
}
readin_dict[2007] = {
        'files':    [f'{2007}/OwnerY{str(2007)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plantcode':'plant_code', 'gencode':'generator_id', 'owner_id':'ownership_id'}
}
readin_dict[2006] = {
        'files':    [f'{2006}/OwnerY{str(2006)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plantcode':'plant_code', 'gencode':'generator_id', 'owner_id':'ownership_id'}
}


# %%
if __name__ == '__main__':
    # read-in parameters
    vars_keep = ['utility_id', 'plant_code', 'generator_id', 'ownership_id', 'status', 
             'owner_name','city_owner', 'state_owner', 'zip_owner', 'percent_owned']
    odf = readin_eia_years(f'{PATH_RAW}eia/f860/', readin_dict, START_YEAR)
    odf['utility_id'] = pd.to_numeric(odf.utility_id).astype('Int64')
    odf['plant_code'] = pd.to_numeric(odf.plant_code).astype('Int64')
    odf['ownership_id'] = pd.to_numeric(odf.ownership_id).astype('Int64')
    odf['zip_owner'] = pd.to_numeric(odf.zip_owner.astype(str).str.strip(), errors='coerce').astype('Int64')
    odf['percent_owned'] = pd.to_numeric(odf.percent_owned.astype(str).str.strip(), errors='coerce').astype('Float64').round(3)
    odf = odf.astype({'generator_id':str, 'status':str, 'owner_name':str, 'state_owner':str})
    odf.to_parquet(PATH_INTERIM + 'eia860_ownership.parquet', index=False)
    odf.to_parquet(PATH_PROCESSED + 'eia860_ownership.parquet', index=False)
    odf['pid'] = odf.utility_id.astype(str) + '.' + odf.plant_code.astype(str)
    odf['gid'] = odf.pid + '.' + odf.generator_id
    odf['oid'] = odf.gid + '.' + odf.ownership_id.astype(str)
    odf = odf.rename(columns={'utility_id':'uid'})
    print('Ownership dataset:')
    summarize_id_counts_byyear(odf.copy(), ['uid', 'pid', 'gid', 'oid', 'ownership_id'])


