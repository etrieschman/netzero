# %%
import pandas as pd
import numpy as np
from tqdm import tqdm
import os


from utils import PATH_EIA, PATH_PROCESSED
from utils import readin_eia, readin_eia_gen, summarize_id_counts_byyear

YR_START, YR_END = 2013, 2021

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 70)

# %%
# UTILITY DATA
vars_keep = ['utility_id', 'utility_name', 'city_util', 'state_util', 'zip_util', 'entity_type']
u_readin_dict={year:{} for year in range(YR_START, YR_END+1)}
for year in u_readin_dict.keys():
    u_readin_dict[year]['vars_keep'] = vars_keep
    if year >= 2013:
        u_readin_dict[year]['path_file'] = f'{year}/1___Utility_Y{year}.xlsx'
        u_readin_dict[year]['excel_params'] = {'header':1}
        u_readin_dict[year]['rename_vars'] = {'city':'city_util', 'state':'state_util', 'zip':'zip_util'}
    elif year >= 2011:
        u_readin_dict[year]['path_file'] = f'{year}/UtilityY{year}.xlsx'
        u_readin_dict[year]['excel_params'] = {'header':1}
        u_readin_dict[year]['rename_vars'] = {'city':'city_util', 'state':'state_util', 'zip5':'zip_util'}
    elif year >= 2010:
        u_readin_dict[year]['path_file'] = f'{year}/UtilityY{year}.xls'
        u_readin_dict[year]['excel_params'] = {'header':0}
        u_readin_dict[year]['rename_vars'] = {'utility_city':'city_util', 'utility_state':'state_util', 'utility_zip5':'zip_util'}
    elif year >= 2009:
        u_readin_dict[year]['path_file'] = f'{year}/UtilityY{str(year)[2:]}.xls'
        u_readin_dict[year]['excel_params'] = {'header':0}
        u_readin_dict[year]['rename_vars'] = {'utility_city':'city_util', 'utility_state':'state_util', 'utility_zip5':'zip_util'}
    elif year >= 2006:
        u_readin_dict[year]['path_file'] = f'{year}/UtilY{str(year)[2:]}.xls'
        u_readin_dict[year]['excel_params'] = {'header':0}
        u_readin_dict[year]['rename_vars'] = {'utilcode':'utility_id', 'utilname':'utility_name', 'city':'city_util', 'state':'state_util', 'zipcode':'zip_util'}

udf = readin_eia(path_folder=f'{PATH_EIA}f860/', readin_dict=u_readin_dict)
# 2010 is a messy year with repeat rows. drop these
# udf = udf.loc[~((udf.year == 2010) & udf.duplicated())]
udf['utility_id'] = pd.to_numeric(udf.utility_id).astype('Int64')
udf['zip_util'] = pd.to_numeric(udf.zip_util.astype(str).str.strip(), errors='coerce').astype('Int64')

# %%
# LOOK AT DUPLICATES
udf['dup_key'] = udf.utility_id.astype(str) + '_' + udf.year.astype(str)
udf['duplicate'] = udf.dup_key.isin(udf.loc[udf.dup_key.duplicated(), 'dup_key'].drop_duplicates())
print('number of dups, by year: ')
display(udf.loc[udf.duplicate].groupby('year')['utility_id'].count())
# NOTE: All the issues are in 2010, and it looks like it's because of differences in naming etc.
udf['num_cols_nan'] = udf.isna().sum(axis=1)
udf['min_num_cols_nan'] = udf.groupby('dup_key')['num_cols_nan'].transform('min')
# drop the duplicate row with the most missing info. If both have smae missing, doesn't matter
udf['dup_keep'] = True
udf.loc[(udf.num_cols_nan != udf.min_num_cols_nan) & udf.duplicate, 'dup_keep'] = False
udf_dedup = udf.loc[udf.dup_keep]
# for remaining, drop duplicate arbitrarily
udf_dedup = udf_dedup.loc[~udf.dup_key.duplicated()]
udf_dedup.drop(columns=['num_cols_nan', 'min_num_cols_nan', 'dup_keep'])
print('number of dups remaining:', udf_dedup.dup_key.duplicated().sum())
# WRITE TO FILE
udf_dedup = udf_dedup.astype({'utility_name':str, 'city_util':str, 'state_util':str})
udf_dedup.to_parquet(PATH_PROCESSED + 'eia_f860_utility.parquet', index=False)

# %%
# PLANT DATA
vars_keep = ['utility_id', 'plant_code', 'plant_name', 'state_plant', 'zip_plant', 'latitude', 'longitude', 'naics_primary']
p_readin_dict={year:{} for year in range(YR_START, YR_END+1)}
for year in p_readin_dict.keys():
    p_readin_dict[year]['vars_keep'] = vars_keep
    if year >= 2013:
        p_readin_dict[year]['path_file'] = f'{year}/2___Plant_Y{year}.xlsx'
        p_readin_dict[year]['excel_params'] = {'header':1}
        p_readin_dict[year]['rename_vars'] = {'state':'state_plant', 'zip':'zip_plant', 'primary_purpose_naics_code':'naics_primary'}
    elif year >= 2012:
        p_readin_dict[year]['path_file'] = f'{year}/PlantY{year}.xlsx'
        p_readin_dict[year]['excel_params'] = {'header':1}
        p_readin_dict[year]['rename_vars'] = {'state':'state_plant', 'zip':'zip_plant', 'primary_purpose_naics_code':'naics_primary'}
    elif year >= 2011:
        p_readin_dict[year]['path_file'] = f'{year}/Plant.xlsx'
        p_readin_dict[year]['excel_params'] = {'header':1}
        p_readin_dict[year]['rename_vars'] = {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'}
    elif year >= 2010:
        p_readin_dict[year]['path_file'] = f'{year}/PlantY{year}.xls'
        p_readin_dict[year]['excel_params'] = {'header':0}
        p_readin_dict[year]['rename_vars'] = {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'}
    elif year >= 2009:
        p_readin_dict[year]['path_file'] = f'{year}/PlantY{str(year)[2:]}.xls'
        p_readin_dict[year]['excel_params'] = {'header':0}
        p_readin_dict[year]['rename_vars'] = {'state':'state_plant', 'zip5':'zip_plant', 'primary_purpose':'naics_primary'}
    elif year >= 2006:
        p_readin_dict[year]['path_file'] = f'{year}/PlantY{str(year)[2:]}.xls'
        p_readin_dict[year]['excel_params'] = {'header':0}
        p_readin_dict[year]['rename_vars'] = {'utilcode':'utility_id', 'plntcode':'plant_code', 'plntname':'plant_name', 'state':'state_plant', 'plntzip':'zip_plant', 'zip5':'zip_plant', 'naics':'naics_primary', 'primary_purpose':'naics_primary'}


pdf = readin_eia(f'{PATH_EIA}f860/', p_readin_dict)
pdf['utility_id'] = pd.to_numeric(pdf.utility_id).astype('Int64')
pdf['plant_code'] = pd.to_numeric(pdf.plant_code).astype('Int64')
pdf['zip_plant'] = pd.to_numeric(pdf.zip_plant.astype(str).str.strip(), errors='coerce').astype('Int64')
pdf = pdf.astype({'plant_name':str, 'state_plant':str, 'latitude':str, 'longitude':str})
pdf.to_parquet(PATH_PROCESSED + 'eia_f860_plant.parquet', index=False)

# %%
# GENERATOR DATA
vars_date = ['operating_month', 'operating_year',
             'current_month', 'current_year',
             'planned_retirement_month', 'planned_retirement_year',
             'retirement_month', 'retirement_year']
vars_keep = ['utility_id', 'plant_code', 'generator_id', 'status', 'ownership', 'sector', 
             'energy_source_1', 'cofire_energy_source_1', 
             'prime_mover', 'nameplate_capacity_mw'] + vars_date

g_readin_dict={year:{} for year in range(YR_START, YR_END+1)}
for year in g_readin_dict.keys():
    g_readin_dict[year]['vars_keep'] = vars_keep
    if year >= 2016:
        g_readin_dict[year]['files'] = [
            f'{year}/3_1_Generator_Y{str(year)}.xlsx', 
            f'{year}/3_2_Wind_Y{year}.xlsx', 
            f'{year}/3_3_Solar_Y{year}.xlsx', 
            f'{year}/3_4_Energy_Storage_Y{year}.xlsx',
            f'{year}/3_5_Multifuel_Y{year}.xlsx']
        g_readin_dict[year]['excel_params'] = {'header':1}
        g_readin_dict[year]['rename_vars'] = {}
    elif year >= 2013:
        g_readin_dict[year]['files'] = [
            f'{year}/3_1_Generator_Y{str(year)}.xlsx',
            f'{year}/3_2_Wind_Y{year}.xlsx',
            f'{year}/3_3_Solar_Y{year}.xlsx',
            f'{year}/3_4_Multifuel_Y{year}.xlsx']
        g_readin_dict[year]['excel_params'] = {'header':1}
        g_readin_dict[year]['rename_vars'] = {}
    elif year >= 2011:
        g_readin_dict[year]['files'] = [f'{year}/GeneratorY{year}.xlsx', f'{year}/MultifuelY{year}.xlsx']
        g_readin_dict[year]['excel_params'] = {'header':1}
        g_readin_dict[year]['rename_vars'] = {'nameplate':'nameplate_capacity_mw', 'sector_number':'sector'}
    elif year >= 2010:
        g_readin_dict[year]['files'] = [f'{year}/GeneratorsY{year}.xls', f'{year}/MultiFuelY{year}.xls']
        g_readin_dict[year]['excel_params'] = {'header':0}
        g_readin_dict[year]['rename_vars'] = {'nameplate':'nameplate_capacity_mw', 'sector_number':'sector'}
    elif year >= 2009:
        g_readin_dict[year]['files'] = [f'{year}/GeneratorY{str(year)[2:]}.xls', f'{year}/MultiFuelY{str(year)[2:]}.xls']
        g_readin_dict[year]['excel_params'] = {'header':0}
        g_readin_dict[year]['rename_vars'] = {'nameplate':'nameplate_capacity_mw', 'sector_number':'sector'}
    elif year >= 2006:
        g_readin_dict[year]['files'] = [f'{year}/GenY{str(year)[2:]}.xls',
                                        f'{year}/MFExistY{str(year)[2:]}.xls',
                                        f'{year}/MFPropY{str(year)[2:]}.xls',
                                        f'{year}/PRGenY{str(year)[2:]}.xls']
        g_readin_dict[year]['excel_params'] = {'header':0}
        g_readin_dict[year]['rename_vars'] = {'utilcode':'utility_id', 'plntcode':'plant_code', 'gencode':'generator_id',
                                              'owner':'ownership', 'primemover':'prime_mover', 'nameplate':'nameplate_capacity_mw',
                                              'insvmonth':'operating_month', 'insvyear':'operating_year',
                                              'retiremonth':'retirement_month', 'retireyear':'retirement_year',
                                              'prop_cofire_energy_source_1':'cofire_energy_source_1',
                                              'proposed_nameplate':'nameplate_capacity_mw', 'proposed_energy_source_1':'energy_source_1'}

gdf = readin_eia_gen(f'{PATH_EIA}f860/', g_readin_dict)
# %%
gdf['utility_id'] = pd.to_numeric(gdf.utility_id, errors='coerce').astype('Int64')
gdf = gdf.loc[gdf.utility_id.notna()]
gdf['plant_code'] = pd.to_numeric(gdf.plant_code, errors='coerce').astype('Int64')
gdf['generator_id'] = gdf.generator_id.astype(str)
gdf['sector'] = pd.to_numeric(gdf.sector, errors='coerce').astype('Int64')
gdf['nameplate_capacity_mw'] = pd.to_numeric(gdf.nameplate_capacity_mw, errors='coerce').astype('Float64')
for var in gdf.columns.intersection(vars_date):
    gdf[var] = pd.to_numeric(gdf[var], errors='coerce').astype('Int64')

# %%
# LOOK AT DUPLICATES
gdf['dup_key'] = gdf[['utility_id', 'plant_code', 'generator_id', 'year']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
gdf['duplicate'] = gdf.dup_key.isin(gdf.loc[gdf.dup_key.duplicated(), 'dup_key'].drop_duplicates())
# NOTE: Duplicates appear across sheets!
gdf_summ = gdf.loc[gdf.duplicate].sort_values(['utility_id', 'plant_code', 'generator_id', 'year'])
display(gdf_summ.groupby(['year', 'gen_category'])['utility_id'].count())
print('total dups:', gdf.duplicate.sum())
# DECISION: Drop duplicates When this is the case
# Take the info from the main, "generator" sheet (and "proposed gen" in 06-08)
gdf['dup_ingen'] = gdf.gen_category.str.startswith(('Gen', 'PRGen'))
gdf['dup_anyingen'] = gdf.groupby('dup_key')['dup_ingen'].transform('sum')
gdf['dup_keep'] = True
gdf.loc[~gdf.dup_ingen & gdf.dup_anyingen, 'dup_keep'] = False
print('remaining dups:', gdf.loc[gdf.dup_keep, 'dup_key'].duplicated().sum())
gdf_dedup = gdf.loc[gdf.dup_keep]
gdf_dedup.drop(columns=['dup_key', 'dup_ingen', 'dup_anyingen', 'dup_keep'], inplace=True)
# WRITE TO FILE
gdf_dedup = gdf_dedup.astype({'status':str, 'prime_mover':str})
gdf_dedup.to_parquet(PATH_PROCESSED + 'eia_f860_generator.parquet', index=False)

# %%
# OWNER DATA
vars_keep = ['utility_id', 'plant_code', 'generator_id', 'ownership_id', 'status', 
             'owner_name','city_owner', 'state_owner', 'zip_owner', 'percent_owned']
o_readin_dict={year:{} for year in range(YR_START, YR_END+1)}
for year in o_readin_dict.keys():
    o_readin_dict[year]['vars_keep'] = vars_keep
    if year >= 2013:
        o_readin_dict[year]['path_file'] = f'{year}/4___Owner_Y{year}.xlsx'
        o_readin_dict[year]['excel_params'] = {'header':1}
        o_readin_dict[year]['rename_vars'] = {'city_owner':'owner_city', 'owner_state':'state_owner', 'owner_zip':'zip_owner'}
    elif year >= 2012:
        o_readin_dict[year]['path_file'] = f'{year}/OwnerY{year}.xlsx'
        o_readin_dict[year]['excel_params'] = {'header':1}
        o_readin_dict[year]['rename_vars'] = {'owner_state':'state_owner'}
    elif year >= 2011:
        o_readin_dict[year]['path_file'] = f'{year}/OwnershipY{year}.xlsx'
        o_readin_dict[year]['excel_params'] = {'header':1}
        o_readin_dict[year]['rename_vars'] = {'owner_state':'state_owner'}
    elif year >= 2010:
        o_readin_dict[year]['path_file'] = f'{year}/OwnerY{year}.xls'
        o_readin_dict[year]['excel_params'] = {'header':0}
        o_readin_dict[year]['rename_vars'] = {'owner_state':'state_owner'}
    elif year >= 2009:
        o_readin_dict[year]['path_file'] = f'{year}/OwnerY{str(year)[2:]}.xls'
        o_readin_dict[year]['excel_params'] = {'header':0}
        o_readin_dict[year]['rename_vars'] = {'owner_state':'state_owner'}
    elif year >= 2006:
        o_readin_dict[year]['path_file'] = f'{year}/OwnerY{str(year)[2:]}.xls'
        o_readin_dict[year]['excel_params'] = {'header':0}
        o_readin_dict[year]['rename_vars'] = {'utilcode':'utility_id', 'plntcode':'plant_code', 'plantcode':'plant_code', 'gencode':'generator_id', 'owner_id':'ownership_id'}


odf = readin_eia(f'{PATH_EIA}f860/', o_readin_dict)
odf['utility_id'] = pd.to_numeric(odf.utility_id).astype('Int64')
odf['plant_code'] = pd.to_numeric(odf.plant_code).astype('Int64')
odf['ownership_id'] = pd.to_numeric(odf.ownership_id).astype('Int64')
odf['zip_owner'] = pd.to_numeric(odf.zip_owner.astype(str).str.strip(), errors='coerce').astype('Int64')
odf['percent_owned'] = pd.to_numeric(odf.percent_owned.astype(str).str.strip(), errors='coerce').astype('Float64').round(3)
odf = odf.astype({'generator_id':str, 'status':str, 'owner_name':str, 'state_owner':str})
odf.to_parquet(PATH_PROCESSED + 'eia_f860_ownership.parquet', index=False)



# %%
# GENERATE UTILITY SUMMARY
# utility data
udf = udf.rename(columns={'utility_id':'uid'})
print('Utility dataset:')
summarize_id_counts_byyear(udf.copy(), ['uid'])

# %%
# GENERATE PLANT SUMMARY
# plant data
pdf['pid'] = pdf.utility_id.astype(str) + '.' + pdf.plant_code.astype(str)
pdf = pdf.rename(columns={'utility_id':'uid'})
print('Plant dataset:')
summarize_id_counts_byyear(pdf.copy(), ['uid', 'pid'])

# %%
# GENERATE GENERATOR SUMMARY
# unit data
gdf['pid'] = gdf.utility_id.astype(str) + '.' + gdf.plant_code.astype(str)
gdf['gid'] = gdf.pid + '.' + gdf.generator_id
gdf = gdf.rename(columns={'utility_id':'uid'})
print('Generator dataset:')
summarize_id_counts_byyear(gdf.copy(), ['uid', 'pid', 'gid'])

# %%
# GENERATE OWNER SUMMARY
# unit data
odf['pid'] = odf.utility_id.astype(str) + '.' + odf.plant_code.astype(str)
odf['gid'] = odf.pid + '.' + odf.generator_id
odf['oid'] = odf.gid + '.' + odf.ownership_id.astype(str)
odf = odf.rename(columns={'utility_id':'uid'})
print('Ownership dataset:')
summarize_id_counts_byyear(odf.copy(), ['uid', 'pid', 'gid', 'oid', 'ownership_id'])
# %%
