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
pd.set_option('display.max_rows', 100)

# read in data
gdf = pd.read_parquet(PATH_PROCESSED + 'df_generators.parquet')
gendf = pd.read_parquet(PATH_PROCESSED + 'df_generation.parquet')
emdf = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'df_plants.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'df_utilities.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'df_owners.parquet')
cdp = pd.read_csv(PATH_RESOURCES + 'cdp_elec.csv', encoding = 'ISO-8859-1')

# %%
# GET RELEVANT IDS (to be updated)
conms = (cdp.conm.drop_duplicates()
         .str.lower()
         .str.replace(' inc| corp| co| grp| group', '', regex=True)
         .str.replace('\\(the\\)', '', regex=True)
         .values
)

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
