# %%
# SETUP
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import format_sumpct_col, format_npct_col, format_mstd_col, format_iqr_col
from utils import sample_summ

# global variables
PATH_PROCESSED = '../data/processed/'
PATH_RESULTS = '../results/analysis/'
ROUND, DENOM = 3, 1e3

pd.set_option('display.max_columns', None)

# %%
# READ IN DATA
g = pd.read_parquet(PATH_PROCESSED + 'df_generators.parquet')
p = pd.read_parquet(PATH_PROCESSED + 'df_plants.parquet')
df = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
gdf = pd.merge(left=g.drop(columns='gid'), right=df.drop(columns='nameplate_capacity_mw'), how='left', on=['year', 'plant_code', 'generator_id'])
assert len(g) == len(gdf)
gdf = pd.merge(left=gdf, right=p[['year', 'plant_code', 'nerc_region']], 
               how='inner', on=['year', 'plant_code'])
gdf.loc[:,'em_cat'] = 'none'
gdf[['has_emissions_epa', 'has_emissions_eia']] = gdf[['has_emissions_epa', 'has_emissions_eia']].astype(bool)
gdf[['co2_mass_tons_gen', 'co2_mass_tons_gen_923']] = gdf[['co2_mass_tons_gen', 'co2_mass_tons_gen_923']].astype(float)
gdf.loc[gdf.has_emissions_epa & ~gdf.has_emissions_eia, 'em_cat'] = 'epa_only'
gdf.loc[~gdf.has_emissions_epa & gdf.has_emissions_eia, 'em_cat'] = 'eia_only'
gdf.loc[gdf.has_emissions_epa & gdf.has_emissions_eia, 'em_cat'] = 'both'
gdf.loc[:,'co2_mass_tons_diff'] = gdf.co2_mass_tons_gen - gdf.co2_mass_tons_gen_923


# SUBSET SAMPLE
gdf_sub = gdf.loc[gdf.in_sample].copy()

# %%
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
    'co2_mass_tons_gen':{
        'params':{'round':ROUND, 'denom':1e6},
        'aggfns':[format_iqr_col, format_mstd_col, format_sumpct_col]
    }, 
    'co2_mass_tons_gen_923':{
        'params':{'round':ROUND, 'denom':1e6},
        'aggfns':[format_iqr_col, format_mstd_col, format_sumpct_col]
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
    summ_cat = sample_summ(gdf_sub, cat=cat+['em_cat'], summdict=summdict)
    summ = pd.concat([summ_cat.reset_index(), summ], axis=0)
summ = summ.fillna('total')
cols_leading = ['year', 'energy_source_1_cat', 'energy_source_1_subcat', 'nerc_region']
cols_tailing = [col for col in summ.columns if col not in cols_leading]
summ = summ[cols_leading + cols_tailing]
summ.to_csv(PATH_RESULTS + 'df_compare_emissions_eia_epa.csv', index=False)


# %%
# PLOT
DENOM = 1e6
FIGSIZE = 2
var_cat = 'energy_source_1_subcat'
cat_sub =['coal', 'natural_gas', 'waste_and_tires', 'biomass']
cat = gdf_sub.loc[(gdf_sub.year == 2021) & 
                  (gdf_sub[var_cat].isin(cat_sub)), var_cat].drop_duplicates().to_list()
fig, ax = plt.subplots(nrows=len(cat), ncols=3, sharey=True, sharex=False, figsize=(FIGSIZE*3, FIGSIZE*len(cat)))
for i, c in enumerate(cat):
    ax[i,0].hist(gdf_sub.loc[(gdf_sub[var_cat] == c) & (gdf_sub.em_cat == 'epa_only'), 'co2_mass_tons_gen']/DENOM, 
                 bins=50, color='C0')
    ax[i,1].hist(gdf_sub.loc[(gdf_sub[var_cat] == c) & (gdf_sub.em_cat == 'both'), 'co2_mass_tons_diff']/DENOM, 
                 bins=50, color='C1')
    ax[i,2].hist(gdf_sub.loc[(gdf_sub[var_cat] == c) & (gdf_sub.em_cat == 'eia_only'), 'co2_mass_tons_gen_923']/DENOM, 
                 bins=50, color='C2')
    ax[i,0].set_yscale('log')
    ax[i,1].set_yscale('log')
    ax[i,2].set_yscale('log')
    ax[i,0].set_ylabel(c)
ax[0,0].set_title('EPA only')
ax[0,1].set_title('Histogram of 2021 CO2 emissions (M metric tons)\nEIA - EPA')
ax[0,2].set_title('EIA only')
plt.savefig(PATH_RESULTS + 'fig_compare_emissions_eia_epa.png', dpi=300)
plt.show()
# %%
