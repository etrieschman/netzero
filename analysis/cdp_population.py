# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm
from thefuzz import fuzz, process

from utils import format_sumpct_col, format_npct_col, format_mstd_col, format_iqr_col
from utils import sample_summ

# global variables
PATH_DATA = '../data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
PATH_INTERIM = PATH_DATA + 'interim/'
PATH_RESOURCES = PATH_DATA + 'resources/'
PATH_RESULTS = '../results/analysis/cdp/'
DENOM, ROUND = 1e6, 2
rpath = Path(PATH_RESULTS)
rpath.mkdir(parents=True, exist_ok=True)

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)

# READ IN DATA
suff = '_llm'
# INTERIM DATA
i_ucdp = pd.read_csv(PATH_INTERIM + 'map_utility_cdp.csv')
i_ocdp = pd.read_csv(PATH_INTERIM + 'map_owner_cdp.csv')
i_ucdpllm = pd.read_csv(PATH_INTERIM + 'map_utility_cdp_llm.csv')
i_ocdpllm = pd.read_csv(PATH_INTERIM + 'map_owner_cdp_llm.csv')
# FINAL DATA
gdf = pd.read_parquet(PATH_PROCESSED + 'df_generators.parquet')
gendf = pd.read_parquet(PATH_PROCESSED + 'df_generation.parquet')
emdf = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'df_plants.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'df_utilities.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'df_owners.parquet')
flag_cdp_gdf = pd.read_parquet(PATH_PROCESSED + 'flag_cdp_generators.parquet')
flag_cdp_udf = pd.read_parquet(PATH_PROCESSED + 'flag_cdp_utilities.parquet')
flag_cdp_gdfllm = pd.read_parquet(PATH_PROCESSED + 'flag_cdp_generators_llm.parquet')
flag_cdp_udfllm = pd.read_parquet(PATH_PROCESSED + 'flag_cdp_utilities_llm.parquet')

# %%
# MAKE ANALYSIS DATASET
# generation
ggdf = pd.merge(left=gdf, right=gendf, how='outer',
                on=['year', 'utility_id', 'plant_code', 'generator_id'])
assert len(gdf) == len(ggdf)
# emissions
ggedf = pd.merge(left=ggdf, how='left',
                 right=emdf[['year', 'plant_code', 'generator_id', 
                             'co2_mass_tons_gen', 'co2_mass_tons_gen_923', 
                             'has_emissions_epa', 'has_emissions_eia']], 
                 on=['year', 'plant_code', 'generator_id'])
assert len(ggdf) == len(ggedf)
# plants
ggepdf = pd.merge(left=ggedf, how='left',
                 right=pdf[['year', 'utility_id', 'plant_code', 'nerc_region']],
                 on=['year', 'utility_id', 'plant_code'])
assert len(ggedf) == len(ggepdf)
# CDP flag
ggepcdf = pd.merge(left=ggepdf, right=flag_cdp_gdf, how='left',
                 on=['year', 'utility_id', 'plant_code', 'generator_id'])
assert len(ggdf) == len(ggepcdf)

ggepcdf['in_cdp'] = ggepcdf.in_cdp | ggepcdf.in_cdp_own

# subset to sample
ggepuodf_samp = ggepcdf.loc[ggepcdf.in_sample]

# %%
# SUMMARIZE FUZZY MATCH RESULTS
sample_cdp_udf = flag_cdp_udf.loc[flag_cdp_udf.utility_id.isin(ggepuodf_samp.utility_id.values)]
sample_cdp_udfllm = flag_cdp_udfllm.loc[flag_cdp_udfllm.utility_id.isin(ggepuodf_samp.utility_id.values)]
# histogram
bins=30
threshold = 90
thresholdllm = 0.98
fig, ax = plt.subplots(nrows=2, sharex=True, sharey=False, figsize=(4,5))
ax[0].hist(sample_cdp_udf.pp_cdp_parname_score, bins=bins, color='C0', density=True)
ax[1].hist(sample_cdp_udf.pp_cdp_subname_score, bins=bins, color='C1', density=True)
ax[0].axvline(x=threshold, color='black', linestyle=':')
ax[1].axvline(x=threshold, color='black', linestyle=':')
ax[0].set_title('Parent name match score')
ax[1].set_title('Subsidiary name match score')
plt.savefig(PATH_RESULTS + 'fig_fuzzyscore_hist.png', dpi=300, bbox_inches='tight')

# %%
# scatterplot
years = [2006, 2012, 2018, 2021]
for year in years:
    plt.figure(figsize=(4,4))
    plt.axvline(x=threshold, color='C1', linestyle='-')
    plt.axhline(y=threshold, color='C1', linestyle='-')
    plt.plot(sample_cdp_udf.loc[sample_cdp_udf.year == year, 'pp_cdp_parname_score'],
            sample_cdp_udf.loc[sample_cdp_udf.year == year, 'pp_cdp_subname_score'],
                marker='.', markersize=1, linestyle='')
    plt.ylabel('Subsidiary name match score')
    plt.xlabel('Parent name match score')
    plt.title(f'EIA utility name matches, {year}')
    plt.savefig(PATH_RESULTS + f'fig_fuzzyscore_scatter_{year}.png', dpi=300, bbox_inches='tight')

# %%
# scatterplot llm vs fuzzy
years = [2006, 2012, 2018, 2021]
for year in years:
    plt.figure(figsize=(4,4))
    plt.axvline(x=threshold, color='C1', linestyle='-')
    plt.axhline(y=thresholdllm, color='C1', linestyle='-')
    plt.plot(sample_cdp_udf.loc[sample_cdp_udf.year == year, 'pp_cdp_subname_score'],
            sample_cdp_udfllm.loc[sample_cdp_udf.year == year, 'pp_cdp_subname_score'],
                marker='.', markersize=1, linestyle='')
    plt.xlabel('Subsidiary name match fuzzy score')
    plt.ylabel('Subsidiary name match LLM score')
    plt.title(f'EIA utility name matches, {year}')
    plt.savefig(PATH_RESULTS + f'fig_fuzzyscore_scatter_{year}.png', dpi=300, bbox_inches='tight')


# %%
# summary table
summ = (sample_cdp_udf
        .groupby(['year', 'in_cdp', 'in_cdp_as_subsid'])
        [['utility_id']]
        .nunique()
        .reset_index())
summ.to_csv(PATH_RESULTS + 'df_utilities_in_cdp.csv', index=False)

# %%
# score of 86
sample_cdp_udf.loc[(sample_cdp_udf.pp_cdp_subname_score == 85.5) & 
                   (sample_cdp_udf.in_cdp)]



# %% 
# SUMMARIZE BY CATEGORY
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

# SUMMARIZE BY CATEGORY
cats = [['year'], ['energy_source_1_cat', 'energy_source_1_subcat', 'year'], ['nerc_region', 'year']]
cats.reverse()
summ = pd.DataFrame()

for cat in cats:
    summ_cat = sample_summ(ggepuodf_samp, cat=cat+['in_cdp'], summdict=summdict)
    summ = pd.concat([summ_cat.reset_index(), summ], axis=0)
summ = summ.fillna('total')
cols_leading = ['year', 'energy_source_1_cat', 'energy_source_1_subcat', 'nerc_region']
cols_leading = [(col, '') for col in cols_leading]
cols_tailing = [col for col in summ.columns if col not in cols_leading]
summ = summ[cols_leading + cols_tailing]
summ.to_csv(PATH_RESULTS + 'df_compare_cdp_in_out.csv', index=False)
# %%
# SUMMARIZE IN A TIGHTER FORMAT
summsub = summ.loc[summ.year.isin([2015, 2021]) & 
                (summ.energy_source_1_cat == 'total') & 
                (summ.nerc_region == 'total')]
summsub.drop(columns=['energy_source_1_cat', 'energy_source_1_subcat', 'nerc_region'], inplace=True)

summsub = summsub.pivot(index='in_cdp', columns='year').T
summsub.to_csv(PATH_RESULTS + 'df_compare_cdp_abbrev.csv', index=True)

# %%
# SUMMARIZE CDP COUNTS
cdp = pd.read_csv(PATH_DATA + 'resources/cdp_elec.csv', encoding='ISO-8859-1')
print('N gvkeys CDP:', len(cdp.gvkey.unique()))
print('N gvkeys mapped:', len(flag_cdp_udf.loc[flag_cdp_udf.in_cdp].gvkey.unique()))
summ = (flag_cdp_udf.loc[flag_cdp_udf.in_cdp]
        .groupby(['year'])
        .agg(unique_gvkey=('gvkey','nunique'),
             unique_utilityid=('utility_id','nunique')))
summ.to_csv(PATH_RESULTS + 'df_unique_gvkeys.csv', index=True)

# %%
summ = (flag_cdp_udf.loc[flag_cdp_udf.in_cdp]
        .groupby(['year', 'gvkey'])['utility_id'].nunique()
        .reset_index())
# Pivot the DataFrame
pivot_summ = summ.pivot(index='year', columns='gvkey', values='utility_id')

# Plotting
pivot_summ.plot(kind='line', linewidth=0.5, legend=False)
plt.title('Count of unique utility_id per GVKEY')
plt.savefig(PATH_RESULTS + 'fig_utility_id_counts.png', dpi=300, bbox_inches='tight')
# %%
