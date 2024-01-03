# %%
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from utils import format_sumpct_col, format_npct_col, format_mstd_col, format_iqr_col
from utils import sample_summ

# global variables
PATH_DATA = '../data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
PATH_RESOURCES = PATH_DATA + 'resources/'
PATH_RESULTS = '../results/analysis/cdp/'
DENOM, ROUND = 1e6, 2
rpath = Path(PATH_RESULTS)
rpath.mkdir(parents=True, exist_ok=True)

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 500)

# read in data
gdf = pd.read_parquet(PATH_PROCESSED + 'df_generators.parquet')
gendf = pd.read_parquet(PATH_PROCESSED + 'df_generation.parquet')
emdf = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'df_plants.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'df_utilities.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'df_owners.parquet')
cdp = pd.read_csv(PATH_RESOURCES + 'cdp_elec.csv', encoding = 'ISO-8859-1')

# %%
# GET RELEVANT IDS
conids = cdp.gvkey.drop_duplicates()
conids = conids.astype(str).str.zfill(6)
conids.to_csv(PATH_RESULTS + 'in_gvkey_list.csv', index=False)

# %%
# READIN PARENT-SUBSIDIARY MAPPING
# This was pulled from WRDS Company Subsidiary data using gvkeys in CDP dataset
csids = pd.read_csv(PATH_RESOURCES + 'out_parent_subsidiary_mapping.csv')
csids['year'] = csids.fdate.str[:4].astype(int)

# %%
# PREPROCESS NAMES
def preprocess_names(x:pd.Series) -> pd.Series:
    x_pp = (x
        .str.lower()
        .str.replace('\.|\,', '', regex=True)
        .str.replace('-', ' ', regex=True)
        .str.replace('\s\/[a-zA-Z]{2}\/?', '', regex=True)
        .str.replace('\\(.*?\\)', '', regex=True)
        .str.strip()
        .str.replace(r'\s\b(corp|company|corporation|co)\b(\s|$)', ' co ', regex=True).str.strip()
        .str.replace(r'\s\b(elec|electric)\b(\s|$)', ' electric ', regex=True).str.strip()
        .str.replace(r'\s\b(inc|incorporated)\b(\s|$)', ' inc ', regex=True).str.strip()
        .str.replace(r'\s\b(grp|group)\b(\s|$)', ' group ', regex=True).str.strip()
        .str.replace(r'\s\b(int|international)\b(\s|$)', ' int ', regex=True).str.strip()
        .str.replace(r'\s\b(ltd|limited)\b(\s|$)', ' ltd ', regex=True).str.strip()
    )

    # one-offs
    x_pp_oo = (x_pp
               .str.replace('pacific gas & electric', 'pg&e', regex=True))
    return x_pp_oo

# cdp
names_cdp = cdp[['gvkey', 'conm']].drop_duplicates().reset_index(drop=True)
names_cdp['name_cdp'] = preprocess_names(names_cdp.conm)

# subsidiary (includes year)
names_cs = csids[['gvkey', 'year', 'coname', 'clean_company']].drop_duplicates().reset_index(drop=True)
names_cs['pp_cdp_parname'] = preprocess_names(names_cs.coname)
names_cs['pp_cdp_subname'] = preprocess_names(names_cs.clean_company)

# EIA datasets (includes year)
udf['pp_utility_name'] = preprocess_names(udf.utility_name)
uids = udf[['utility_id', 'utility_name', 'pp_utility_name', 'year']].drop_duplicates().reset_index(drop=True)
odf['pp_owner_name'] = preprocess_names(odf.owner_name)
oids = odf[['utility_id', 'plant_code', 'generator_id', 'ownership_id', 'owner_name', 'pp_owner_name', 'year']].drop_duplicates().reset_index(drop=True)

# %%
# FLAG EIA IDENTIFIERS (ID-YEAR) WITH CDP COMPANIES
from thefuzz import fuzz, process
import concurrent
import time

def fuzzy_search(df_target, df_search, key_target, key_search, 
                 scorer, num_matches, threshold):
    # Step 0: fuzzy search
    searches = df_search[key_search].drop_duplicates()
    searches = searches.tolist()
    targets = df_target[key_target].drop_duplicates()
    target_map = targets.apply(lambda x: process.extract(x, searches, scorer=scorer, limit=num_matches))

    # Step 1: Expand the list of tuples into separate columns
    df_expanded = pd.DataFrame(target_map.tolist(), index=targets)
    df_expanded.columns = [f'pair_{i}' for i in range(num_matches)]
    
    # Step 2: Transform to Long Form
    df_long = df_expanded.melt(var_name='pair', value_name='name_score', ignore_index=False)
    
    # Step 3: Split the Tuples into Two Columns
    df_long[[key_search, f'{key_search}_score']] = pd.DataFrame(df_long['name_score'].tolist(), index=df_long.index)
    
    # Drop columns and rows
    df_long = df_long.drop(['name_score', 'pair'], axis=1)
    df_long_sub = df_long.loc[df_long[f'{key_search}_score'] >= threshold].reset_index()

    return df_long_sub

def fuzzy_search_years(df_target, df_search, key_target, key_search, 
                       scorer=fuzz.WRatio, num_matches=1, threshold=95):
    dfs = pd.DataFrame()
    for y in tqdm(df_target.year.unique()):
        df = fuzzy_search(df_target.loc[df_target.year == y],
                        df_search.loc[df_search.year == y],
                        key_target, key_search, scorer, num_matches, threshold)
        df['year'] = y
        dfs = pd.concat([df, dfs], ignore_index=True)
    return dfs

# %%
# CREATE MAPPINGS
scorer = fuzz.WRatio
num_matches = 1
threshold = 95

updf = fuzzy_search_years(uids, names_cs, 'pp_utility_name', 'pp_cdp_parname',
                        scorer, num_matches, threshold)
usdf = fuzzy_search_years(uids, names_cs, 'pp_utility_name', 'pp_cdp_subname',
                        scorer, num_matches, threshold)

# %%
uspdf = pd.merge(
    left=usdf, how='left',
    right=names_cs[['year', 'pp_cdp_parname', 'pp_cdp_subname', 'gvkey']].drop_duplicates(),
    on=['year', 'pp_cdp_subname'])

# %%
opdf = fuzzy_search_years(oids, names_cs, 'pp_utility_name', 'pp_cdp_parname',
                        scorer, num_matches, threshold)
osdf = fuzzy_search_years(oids, names_cs, 'pp_utility_name', 'pp_cdp_subname',
                        scorer, num_matches, threshold)


# %%
# MERGE CDP FLAGS ONTO UTILITIES AND GENERATORS
# merge on cdp parent name
udf_p = pd.merge(left=udf, right=updf, how='outer', 
                 on=['year', 'pp_utility_name'])
assert len(udf) == len(udf_p)
# merge on cdp subsidiary name (this will create dups b/c of partial ownership)
# NOTE: because subsidiaries are reported for 50% ownership, we expect dups of only magnitute 2
# i.e., can't have more than 2 companies that own 50%
udf_ps = pd.merge(left=udf_p, right=uspdf, how='outer',
                  on=['year', 'pp_utility_name'],
                  suffixes=['', '_s'])







# %%
# MAKE ANALYSIS DATASET
gedf = pd.merge(left=gdf.drop(columns='gid'), right=emdf.drop(columns='nameplate_capacity_mw'), how='left', on=['year', 'plant_code', 'generator_id'])
assert len(gdf) == len(gdf)
gepdf = pd.merge(left=gedf, right=pdf[['year', 'plant_code', 'nerc_region']], 
               how='inner', on=['year', 'plant_code'])
gepgdf = pd.merge(left=gepdf, right=gendf, how='left',
                  on=['year', 'utility_id', 'plant_code', 'generator_id'])
gepgdf['in_cdp'] = gepdf.utility_id.isin(uid)


# %% 
# SUMMARIZE
ROUND = 3
summdict = {
    'age':{
        'params':{'round':ROUND, 'denom':1},
        'aggfns':[format_iqr_col, format_mstd_col]
    },
    'nameplate_capacity_mw':{
        'params':{'round':ROUND, 'denom':1},
        'aggfns':[format_iqr_col, format_mstd_col, format_sumpct_col]
    },
    'net_gen_tot_an':{
        'params':{'round':ROUND, 'denom':1e6},
        'aggfns':[format_iqr_col, format_mstd_col, format_sumpct_col]
    },
    'co2_mass_tons_gen':{
        'params':{'round':ROUND, 'denom':1e6},
        'aggfns':[format_iqr_col, format_mstd_col, format_sumpct_col]
    }, 
    'co2_mass_tons_gen_923':{
        'params':{'round':ROUND, 'denom':1e6},
        'aggfns':[format_iqr_col, format_mstd_col, format_sumpct_col]
    },
    'utility_id':{
        'params':{'round':ROUND, 'denom':1},
        'aggfns':[format_npct_col]
    },
    'plant_code':{
        'params':{'round':ROUND, 'denom':1},
        'aggfns':[format_npct_col]
    }, 
    'gid':{
        'params':{'round':ROUND, 'denom':1},
        'aggfns':[format_npct_col]
    }
}


# %%
# SUMMARIZE BY CATEGORY
cats = [['year'], ['energy_source_1_cat', 'energy_source_1_subcat', 'year'], ['nerc_region', 'year']]
cats.reverse()
summ = pd.DataFrame()

for cat in cats:
    summ_cat = sample_summ(gepgdf, cat=cat+['in_cdp'], summdict=summdict)
    summ = pd.concat([summ_cat.reset_index(), summ], axis=0)
summ = summ.fillna('total')
cols_leading = ['year', 'energy_source_1_cat', 'energy_source_1_subcat', 'nerc_region']
cols_tailing = [col for col in summ.columns if col not in cols_leading]
summ = summ[cols_leading + cols_tailing]
summ.to_csv(PATH_RESULTS + 'df_compare_cdp_in_out.csv', index=False)
# %%
