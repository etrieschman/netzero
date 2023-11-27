# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

from utils_transform import readin_eia_years
from utils_summ import summarize_id_counts_byyear

# %%
# OWNERSHIP DATA
readin_dict = {}
readin_dict[2021] = {
    'files': [f'{2021}/4___Owner_Y{2021}.xlsx'],
    'excel_params':{'header':1, 'sheet_name':None},
    'rename_vars':{'city_owner':'owner_city', 'owner_state':'state_owner', 'owner_zip':'zip_owner'}
}
# cut corner: 2013+ is all the same
for yr in range(2013, 2021+1):
    readin_dict[yr] = readin_dict[2021].copy()
    readin_dict[yr]['files'] = [f'{yr}/4___Owner_Y{yr}.xlsx']

readin_dict[2012] = {
    'files':    [f'{2012}/OwnerY{2012}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2011] = {
    'files':    [f'{2011}/OwnershipY{2011}.xlsx'],
    'excel_params': {'header':1, 'sheet_name':None},
    'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2010] = {
        'files':    [f'{2010}/OwnerY{2010}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2009] = {
        'files':    [f'{2009}/OwnerY{str(2009)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'owner_state':'state_owner'}
}
readin_dict[2008] = {
        'files':    [f'{2008}/OwnerY{str(2008)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plantcode':'plant_code', 'gencode':'generator_id', 'owner_id':'ownership_id'}
}
readin_dict[2007] = {
        'files':    [f'{2007}/OwnerY{str(2007)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plantcode':'plant_code', 'gencode':'generator_id', 'owner_id':'ownership_id'}
}
readin_dict[2006] = {
        'files':    [f'{2006}/OwnerY{str(2006)[2:]}.xls'],
        'excel_params': {'header':0, 'sheet_name':None},
        'rename_vars':  {'utilcode':'utility_id', 'plntcode':'plant_code', 'plantcode':'plant_code', 'gencode':'generator_id', 'owner_id':'ownership_id'}
}


# %%
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('transform_eia860_ownership')

    year_start = snakemake.params.year_start
    year_end = snakemake.params.year_end
    path_raw = snakemake.params.indir
    intfile = snakemake.output.intfile
    outfile = snakemake.output.outfile

    # read-in parameters
    print('Reading in data...')
    vars_keep = ['utility_id', 'plant_code', 'generator_id', 'ownership_id', 'status', 
             'owner_name','city_owner', 'state_owner', 'zip_owner', 'percent_owned']
    odf = readin_eia_years(path_raw, readin_dict, year_start)
    odf['utility_id'] = pd.to_numeric(odf.utility_id).astype('Int64')
    odf['plant_code'] = pd.to_numeric(odf.plant_code).astype('Int64')
    odf['ownership_id'] = pd.to_numeric(odf.ownership_id).astype('Int64')
    odf['zip_owner'] = pd.to_numeric(odf.zip_owner.astype(str).str.strip(), errors='coerce').astype('Int64')
    odf['percent_owned'] = pd.to_numeric(odf.percent_owned.astype(str).str.strip(), errors='coerce').astype('Float64').round(3)
    odf = odf.astype({'generator_id':str, 'status':str, 'owner_name':str, 'state_owner':str})
    odf.to_parquet(intfile, index=False)

    # drop unnecessary columns
    newcols = ['year', 'sheet', 'file']
    odf = odf.drop(columns=odf.columns.difference(vars_keep + newcols))

    print('Cleaning key variables...')
    summ_dict = {}
    # FIRST: OWNERSHIP PERCENT
    # TODO: 2010 is a messy year that needs to be reconciled
    odf['pct_ownership_total'] = (
        odf.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)
        ['percent_owned'].transform('sum').round(5))
    odf['owner_count'] = (
        odf.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)
        ['percent_owned'].transform('count'))
    print("Before cleaning: count of generators where ownership doesn't add up:")
    print(odf.loc[np.abs(odf.pct_ownership_total - 1.) > 0.05].groupby('year')['generator_id'].agg(['count', 'nunique']))
    # 0. Divide by 100 in settings where decimal is off
    mask_decimal = odf.pct_ownership_total >= 50.
    summ_dict['decimal off'] = mask_decimal.sum()
    odf.loc[mask_decimal, 'percent_owned'] /= 100.
    # 1. Scale if close to 100
    mask_almost_100pct = np.abs(odf.pct_ownership_total - 1) <= 0.06
    summ_dict['close_to_1'] = mask_almost_100pct.sum()
    odf.loc[mask_almost_100pct, 'percent_owned'] /= odf.loc[mask_almost_100pct, 'pct_ownership_total']
    # 2. If ownership is missing, but no one else
    mask_miss_only_owner = (odf.owner_count == 0)
    mask_miss_other_owners = (odf.percent_owned.isna()) & (odf.owner_count > 0)
    odf.loc[mask_miss_only_owner, 'percent_owned'] = 1.
    odf.loc[mask_miss_other_owners, 'percent_owned'] = 0.
    summ_dict['missing_only_owner'] = mask_miss_only_owner.sum()
    summ_dict['missing_other_owner'] = mask_miss_other_owners.sum()
    odf['pct_ownership_total'] = (
        odf.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)
        ['percent_owned'].transform('sum').round(5))
    odf['owner_count'] = (
        odf.groupby(['utility_id', 'plant_code', 'generator_id', 'year'], dropna=False)
        ['percent_owned'].transform('count'))
    print('Cleaning summary:\n', summ_dict)
    print("After cleaning: count of generators where ownership doesn't add up:")
    print(odf.loc[np.abs(odf.pct_ownership_total - 1.) > 0.05].groupby('year')['generator_id'].agg(['count', 'nunique']))
    # TODO: need to handle 2010, it's really messy
    
    # SECOND: OWNERSHIP STATUS
    mask_o_is_u = odf.ownership_id == odf.utility_id
    odf.loc[mask_o_is_u & (odf.percent_owned == 1.), 'ownership'] = 'S' # Single ownership by respondent
    odf.loc[~mask_o_is_u & (odf.percent_owned == 1.), 'ownership'] = 'W' # Wholly owned by an entity other than respondent
    odf.loc[~mask_o_is_u & (odf.percent_owned != 1.), 'ownership'] = 'J' # Jointly owned with another entity
    odf = odf.drop(columns=['pct_ownership_total', 'owner_count'])
    odf.to_parquet(outfile, index=False)

    print('Summarizing unique identifiers...')
    odf['pid'] = odf.utility_id.astype(str) + '.' + odf.plant_code.astype(str)
    odf['gid'] = odf.pid + '.' + odf.generator_id
    odf['oid'] = odf.gid + '.' + odf.ownership_id.astype(str)
    odf = odf.rename(columns={'utility_id':'uid'})
    print('Ownership dataset:')
    print(summarize_id_counts_byyear(odf.copy(), ['uid', 'pid', 'gid', 'oid', 'ownership_id']))