# %%
import pandas as pd
import numpy as np
from tqdm import tqdm

# from process_epa_data import readin_data

PATH_DATA = '../data/'
PATH_EIA = PATH_DATA + 'eia/'
PATH_EPA = PATH_DATA + 'epa/'

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 25)


# %%
# READIN DATA
years = range(2018, 2022)

# %%
# UTILITY DATA
cols_keep = ['utility_id', 'utility_name', 'city', 'state', 'entity_type']
udf = pd.DataFrame({})
for y in tqdm(years):
    df = pd.read_excel(f'{PATH_EIA}f860/{y}/1___Utility_Y{y}.xlsx', header=1)
    df.columns = df.columns.str.lower().str.replace(' ', '_')
    df = df[cols_keep]
    df['year'] = y
    udf = pd.concat([df, udf], axis=0, ignore_index=True)
udf

# %%
# PLANT DATA
cols_keep = ['utility_id', 'utility_name', 'plant_code', 'plant_name',
       'city', 'state', 'zip', 'county', 'latitude',
       'longitude', 'primary_purpose_naics_code']
pdf = pd.DataFrame({})
for y in tqdm(years):
    df = pd.read_excel(f'{PATH_EIA}f860/{y}/2___Plant_Y{y}.xlsx', header=1)
    df.columns = df.columns.str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')','')
    df = df[cols_keep]
    df['year'] = y
    pdf = pd.concat([df, pdf], axis=0, ignore_index=True)
pdf

# %%
# GENERATOR DATA
gen_tags = ['1_Generator', '2_Wind', '3_Solar', '4_Energy_Storage', '5_Multifuel']
cols_keep = [
    'utility_id', 'utility_name', 'plant_code', 'plant_name', 'state',
       'county', 'generator_id', 'technology', 'prime_mover',
       'nameplate_capacity_mw',
]
gdf = pd.DataFrame({})
for y in tqdm(years):
    for g in tqdm(gen_tags):
        df = pd.read_excel(f'{PATH_EIA}f860/{y}/3_{g}_Y{y}.xlsx', header=1)
        df.columns = df.columns.str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')','')
        df = df[cols_keep]
        df['year'] = y
        gdf = pd.concat([df, gdf], axis=0, ignore_index=True)
gdf


# %%
# SUMMARIZE BY COUNTS
def summarize_id_counts(df, ids):
    for id in ids:
        prev = set()
        counts = {}
        counts['n'], counts['n_new'], counts['n_drop'] = [], [], []
        for i, row in df.iterrows():
            curr = row[id]
            # Calculate overlap and new IDs
            counts['n'] += [len(curr)]
            counts['n_new'] += [len(curr.difference(prev))]
            counts['n_drop'] += [len(prev.difference(curr))]
            # Update prev_ids for the next iteration
            prev = curr

        # Add new columns to the DataFrame
        for k, v in counts.items():
            df[f'{id}_{k}'] = v
        df = df.drop(columns=[id])
    return df

# %%
# GENERATE UTILITY SUMMARY
# utility data
ludf = (udf.groupby('year')[['utility_id']]
           .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
ludf = ludf.rename(columns={'utility_id':'uid'})
print('Utility dataset:')
summarize_id_counts(ludf.copy(), ['uid'])

# %%
# GENERATE PLANT SUMMARY
# plant data
pdf['pid'] = pdf.utility_id.astype(str) + '.' + pdf.plant_code.astype(str)
lpdf = (pdf.groupby('year')[['utility_id', 'pid']]
        .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
lpdf = lpdf.rename(columns={'utility_id':'uid'})
print('Plant dataset:')
summarize_id_counts(lpdf.copy(), ['uid', 'pid'])

# %%
# GENERATE GENERATOR SUMMARY
# unit data
gdf['pid'] = gdf.utility_id.astype(str) + '.' + gdf.plant_code.astype(str)
gdf['gid'] = gdf.pid + '.' + gdf.generator_id
lgdf = (gdf.groupby(['year'])[['utility_id', 'pid', 'gid']]
        .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
lgdf = lgdf.rename(columns={'utility_id':'uid'})
# .reset_index(drop=True).sort_values(['technology', 'year'])
print('Generator dataset:')
summarize_id_counts(lgdf.copy(), ['uid', 'pid', 'gid'])

# %%
