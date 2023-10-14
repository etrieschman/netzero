# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

from utils_transform import (
    PATH_RAW, PATH_INTERIM, PATH_PROCESSED, START_YEAR, END_YEAR)
from utils_eia import readin_eia_years
from utils_summ import summarize_id_counts_byyear

# %%
# OWNERSHIP DATA
readin_dict = {}
readin_dict[END_YEAR] = {
    'files':    [f'{END_YEAR}/2___Plant_Y{END_YEAR}.xlsx'],
    'excel_params': {'header':1},
    'rename_vars':  {'state':'state_plant', 'zip':'zip_plant', 'primary_purpose_naics_code':'naics_primary'}
}
# cut corner: 2013+ is all the same
for yr in range(2013, END_YEAR+1):
    readin_dict[yr] = readin_dict[END_YEAR].copy()
    readin_dict[yr]['files'] = [f'{yr}/2___Plant_Y{yr}.xlsx']

readin_dict[2012] = {
    'files':    [f'{2012}/PlantY{2012}.xlsx'],
    'excel_params': {'header':1},
    'rename_vars':  {'state':'state_plant', 'zip':'zip_plant', 'primary_purpose_naics_code':'naics_primary'}
}
readin_dict[2011] = {
    'files':    [f'{2011}/Plant.xlsx'],
    'excel_params': {'header':1},
    'rename_vars':  {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'}
}
readin_dict[2010] = {
        'files':    [f'{2010}/PlantY{2010}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'},
}
readin_dict[2009] = {
        'files':    [f'{2009}/PlantY{str(2009)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'}
}
readin_dict[2008] = {
        'files':    [f'{2008}/PlantY{str(2008)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plntname':'plant_name', 'state':'state_plant', 'plntzip':'zip_plant', 'zip5':'zip_plant', 'naics':'naics_primary', 'primary_purpose':'naics_primary'}
}
readin_dict[2007] = {
        'files':    [f'{2007}/PlantY{str(2007)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plntname':'plant_name', 'state':'state_plant', 'plntzip':'zip_plant', 'zip5':'zip_plant', 'naics':'naics_primary', 'primary_purpose':'naics_primary'}
}
readin_dict[2006] = {
        'files':    [f'{2006}/PlantY{str(2006)[2:]}.xls'],
        'excel_params': {'header':0},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plntname':'plant_name', 'state':'state_plant', 'plntzip':'zip_plant', 'zip5':'zip_plant', 'naics':'naics_primary', 'primary_purpose':'naics_primary'}
}

# %%
if __name__ == '__main__':
    # read-in parameters
    vars_keep = ['utility_id', 'plant_code', 'plant_name', 'state_plant', 
                 'zip_plant', 'latitude', 'longitude', 'naics_primary']
    pdf = readin_eia_years(f'{PATH_RAW}eia/f860/', readin_dict, START_YEAR)
    pdf['utility_id'] = pd.to_numeric(pdf.utility_id).astype('Int64')
    pdf['plant_code'] = pd.to_numeric(pdf.plant_code).astype('Int64')
    pdf['zip_plant'] = pd.to_numeric(pdf.zip_plant.astype(str).str.strip(), errors='coerce').astype('Int64')
    pdf = pdf.astype({'plant_name':str, 'state_plant':str, 'latitude':str, 'longitude':str})
    col_str = pdf.columns[pdf.dtypes == 'object']
    pdf[col_str] = pdf[col_str].astype(str)
    pdf.to_parquet(PATH_INTERIM + 'eia860_plant.parquet', index=False)
    pdf.to_parquet(PATH_PROCESSED + 'eia860_plant.parquet', index=False)
    pdf['pid'] = pdf.utility_id.astype(str) + '.' + pdf.plant_code.astype(str)
    pdf = pdf.rename(columns={'utility_id':'uid'})
    print('Plant dataset:')
    print(summarize_id_counts_byyear(pdf.copy(), ['uid', 'pid']))



# %%