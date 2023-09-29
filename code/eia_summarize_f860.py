# %%
import pandas as pd

from utils import PATH_INTERIM, PATH_PROCESSED
from utils_summ import summarize_id_counts_byyear


# %%
# GENERATE UTILITY SUMMARY
# utility data
udf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_utility.parquet')
udf = udf.rename(columns={'utility_id':'uid'})
print('Utility dataset:')
summarize_id_counts_byyear(udf.copy(), ['uid'])

# %%
# GENERATE PLANT SUMMARY
# plant data
pdf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_plant.parquet')
pdf['pid'] = pdf.utility_id.astype(str) + '.' + pdf.plant_code.astype(str)
pdf = pdf.rename(columns={'utility_id':'uid'})
print('Plant dataset:')
summarize_id_counts_byyear(pdf.copy(), ['uid', 'pid'])

# %%
# GENERATE GENERATOR SUMMARY
# unit data
gdf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_generator.parquet')
gdf['pid'] = gdf.utility_id.astype(str) + '.' + gdf.plant_code.astype(str)
gdf['gid'] = gdf.pid + '.' + gdf.generator_id
gdf = gdf.rename(columns={'utility_id':'uid'})
print('Generator dataset:')
summarize_id_counts_byyear(gdf.copy(), ['uid', 'pid', 'gid'])

# %%
# GENERATE OWNER SUMMARY
# unit data
odf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_ownership.parquet')
odf['pid'] = odf.utility_id.astype(str) + '.' + odf.plant_code.astype(str)
odf['gid'] = odf.pid + '.' + odf.generator_id
odf['oid'] = odf.gid + '.' + odf.ownership_id.astype(str)
odf = odf.rename(columns={'utility_id':'uid'})
print('Ownership dataset:')
summarize_id_counts_byyear(odf.copy(), ['uid', 'pid', 'gid', 'oid', 'ownership_id'])




# %%
# EXPLORE: UNDERSTAND UTILITY DUPLICATES
udf = pd.read_parquet(PATH_INTERIM + 'eia_f860_utility.parquet')
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

# %%
# EXPLORE: UNDERSTAND GENERATOR DUPLICATES
gdf = pd.read_parquet(PATH_INTERIM + 'eia_f860_generator.parquet')
gdf['dup_key'] = gdf[['utility_id', 'plant_code', 'generator_id', 'year']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
gdf['duplicate'] = gdf.dup_key.isin(gdf.loc[gdf.dup_key.duplicated(), 'dup_key'].drop_duplicates())
# NOTE: Duplicates appear across sheets!
gdf_summ = gdf.loc[gdf.duplicate].sort_values(['utility_id', 'plant_code', 'generator_id', 'year'])
display(gdf_summ.groupby(['year', 'gen_category'])['utility_id'].count())
print('total dups:', gdf.duplicate.sum())
# DECISION: Drop duplicates When this is the case
# Take the info from the main, "generator" sheet (and "proposed gen" in 06-08)
gdf['dup_ingen'] = gdf.gen_category.str.startswith(('gen', 'prgen'))
gdf['dup_anyingen'] = gdf.groupby('dup_key')['dup_ingen'].transform('sum')
gdf['dup_keep'] = True
gdf.loc[~gdf.dup_ingen & gdf.dup_anyingen, 'dup_keep'] = False
print('remaining dups:', gdf.loc[gdf.dup_keep, 'dup_key'].duplicated().sum())
gdf_dedup = gdf.loc[gdf.dup_keep]
