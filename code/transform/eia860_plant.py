# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from pathlib import Path

from utils_transform import readin_eia_years
from utils_summ import summarize_id_counts_byyear

# %%
# OWNERSHIP DATA
readin_dict = {}
readin_dict[2021] = {
    'files':    [f'{2021}/2___Plant_Y{2021}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  {'state':'state_plant', 'zip':'zip_plant', 'primary_purpose_naics_code':'naics_primary'}
}
# cut corner: 2013+ is all the same
for yr in range(2013, 2021+1):
    readin_dict[yr] = readin_dict[2021].copy()
    readin_dict[yr]['files'] = [f'{yr}/2___Plant_Y{yr}.xlsx']

readin_dict[2012] = {
    'files':    [f'{2012}/PlantY{2012}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  {'state':'state_plant', 'zip':'zip_plant', 'primary_purpose_naics_code':'naics_primary'}
}
readin_dict[2011] = {
    'files':    [f'{2011}/Plant.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'}
}
readin_dict[2010] = {
        'files':    [f'{2010}/PlantY{2010}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'},
}
readin_dict[2009] = {
        'files':    [f'{2009}/PlantY{str(2009)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'}
}
readin_dict[2008] = {
        'files':    [f'{2008}/PlantY{str(2008)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plntname':'plant_name', 'state':'state_plant', 'plntzip':'zip_plant', 'zip5':'zip_plant', 'naics':'naics_primary', 'primary_purpose':'naics_primary'}
}
readin_dict[2007] = {
        'files':    [f'{2007}/PlantY{str(2007)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plntname':'plant_name', 'state':'state_plant', 'plntzip':'zip_plant', 'zip5':'zip_plant', 'naics':'naics_primary', 'primary_purpose':'naics_primary'}
}
readin_dict[2006] = {
        'files':    [f'{2006}/PlantY{str(2006)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plntname':'plant_name', 'state':'state_plant', 'plntzip':'zip_plant', 'zip5':'zip_plant', 'naics':'naics_primary', 'primary_purpose':'naics_primary'}
}

vars_keep = [
    'utility_id', 
    'plant_code', 'plant_name',
    'street_address', 'city', 'state_plant', 'zip_plant',
    'latitude', 'longitude', 'nerc_region',
    'balancing_authority_code', 'balancing_authority_name',
    'name_of_water_source',
    'naics_primary', 'regulatory_status',
    'sector', 'sector_name',
    'transmission_or_distribution_system_owner',
    'transmission_or_distribution_system_owner_id',
    'transmission_or_distribution_system_owner_state']

# %%
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('transform_eia860_plant')

    year_start = snakemake.params.year_start
    year_end = snakemake.params.year_end
    path_raw = Path(snakemake.params.indir)
    path_results = Path(snakemake.params.resultsdir)
    path_results.mkdir(parents=True, exist_ok=True)
    intfile = Path(snakemake.output.intfile)
    intfile.parent.mkdir(parents=True, exist_ok=True)
    outfile = Path(snakemake.output.outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    # read-in parameters
    print('Reading in data...')
    pdf = readin_eia_years(path_raw, readin_dict, year_start, keep_vars=None)
    pdf['utility_id'] = pd.to_numeric(pdf.utility_id).astype('Int64')
    pdf['plant_code'] = pd.to_numeric(pdf.plant_code).astype('Int64')
    pdf['zip_plant'] = pd.to_numeric(pdf.zip_plant.astype(str).str.strip(), errors='coerce').astype('Int64')
    pdf['sector'] = pd.to_numeric(pdf.sector, errors='coerce').astype('Int64')
    pdf = pdf.astype({'plant_name':str, 'state_plant':str, 'latitude':str, 'longitude':str})
    col_str = pdf.columns[pdf.dtypes == 'object']
    pdf[col_str] = pdf[col_str].astype(str)
    print('Writing to file...')
    pdf.to_parquet(intfile, index=False)

    # drop variables
    new_cols = ['year', 'file', 'sheet']
    pdf = pdf.drop(columns=pdf.columns.difference(vars_keep + new_cols))
    pdf.to_parquet(outfile, index=False)
    pdf['pid'] = pdf.utility_id.astype(str) + '.' + pdf.plant_code.astype(str)
    pdf = pdf.rename(columns={'utility_id':'uid'})
    print('Summarizing unique identifiers...')
    print('Plant dataset:')
    out = summarize_id_counts_byyear(pdf.copy(), ['uid', 'pid'])
    out.to_csv(path_results / 'df_summ_final.csv')



# %%
