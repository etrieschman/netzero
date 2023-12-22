# %%
import pandas as pd
from pathlib import Path

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
        .str.replace('\s\/[a-zA-Z]{2}\/?', '', regex=True)
        .str.replace('\\(the\\)', '', regex=True)
        .str.replace(' inc$| corp$| company$| corporation$| corp$| co$', '', regex=True)
        .str.replace(' grp$| industries$| international$', '', regex=True)
        .str.replace(' inc$| corp$| company$| corporation$| corp$| co$', '', regex=True)
        .str.replace(' grp$| industries$| international$', '', regex=True)
        .str.replace(' llc| ltd', '', regex=True)
        .str.strip()
    )
    return x_pp

# cdp
names_cdp = cdp[['gvkey', 'conm']].drop_duplicates().reset_index(drop=True)
names_cdp['name_cdp'] = preprocess_names(names_cdp.conm)

# subsidiary (includes year)
names_cs = csids[['gvkey', 'year', 'coname', 'clean_company']].drop_duplicates().reset_index(drop=True)
names_cs['name_cs_parent'] = preprocess_names(names_cs.coname)
names_cs['name_cs_subsid'] = preprocess_names(names_cs.clean_company)

# EIA datasets (includes year)
uids = udf[['utility_id', 'utility_name', 'year']].drop_duplicates().reset_index(drop=True)
uids['name_utility'] = preprocess_names(uids.utility_name)
oids = odf[['utility_id', 'plant_code', 'generator_id', 'ownership_id', 'owner_name', 'year']].drop_duplicates().reset_index(drop=True)
oids['name_owner'] = preprocess_names(oids.owner_name)
pids = pdf[['utility_id', 'plant_code', 'plant_name', 'year']].drop_duplicates().reset_index(drop=True)
pids['name_plant'] = preprocess_names(pids.plant_name)

# %%
# FLAG EIA IDENTIFIERS (ID-YEAR) WITH CDP COMPANIES
from thefuzz import fuzz
from thefuzz import process
import jellyfish

df_target = uids
target_id = 'name_utility'
df_search = names_cs.copy()
search_id = 'name_cs_parent'
num_matches = 5
threshold = 95

# fuzzy search
searches = df_search[search_id].drop_duplicates()
search = searches.tolist()
targets = df_target[target_id].drop_duplicates()
target_map = targets.apply(lambda x: process.extract(x, search, scorer=fuzz.WRatio, limit=num_matches))

# %%
# Step 1: Expand the list of tuples into separate columns
df_expanded = pd.DataFrame(target_map.tolist(), index=targets)
df_expanded.columns = [f'pair_{i}' for i in range(num_matches)]
# Step 2: Transform to Long Form
df_long = df_expanded.melt(var_name='pair', value_name='name_score', ignore_index=False)
# Step 3: Split the Tuples into Two Columns
df_long[[search_id, 'lscore']] = pd.DataFrame(df_long['name_score'].tolist(), index=df_long.index)
# Drop the original tuple column
df_long = df_long.drop('name_score', axis=1).reset_index()
# drop pairs where score is low
df_long_sub = df_long.loc[df_long.lscore >= threshold]
df_long_sub.sort_values(by=['name_cs_parent', 'lscore'], ascending=False)







# %%
# get operated IDs
uid_op = udf.loc[udf.utility_name.str.lower().str.contains('|'.join(conms)),
                 ['utility_id', 'utility_name']].drop_duplicates().reset_index(drop=True)
uid_op = set(uid_op.utility_id)

# get owned IDs
uid_own = odf.loc[odf.owner_name.str.lower().str.contains('|'.join(conms)),
                 ['utility_id', 'ownership_id', 'owner_name']].drop_duplicates().reset_index(drop=True)
uid_own = set(uid_own.utility_id)

# get any IDs
uid = uid_op.union(uid_own)

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
