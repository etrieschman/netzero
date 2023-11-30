# %%
# SETUP
import pandas as pd
import numpy as np

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
gdf.loc[:,'em_cat'] = 'both'
gdf[['has_emissions_epa', 'has_emissions_eia']] = gdf[['has_emissions_epa', 'has_emissions_eia']].astype(bool)
gdf[['co2_mass_tons_gen', 'co2_mass_tons_gen_923']] = gdf[['co2_mass_tons_gen', 'co2_mass_tons_gen_923']].astype(float)
gdf.loc[gdf.has_emissions_epa & ~gdf.has_emissions_eia, 'em_cat'] = 'epa_only'
gdf.loc[~gdf.has_emissions_epa & gdf.has_emissions_eia, 'em_cat'] = 'eia_only'
gdf.loc[:,'co2_mass_diff'] = gdf.co2_mass_tons_gen - gdf.co2_mass_tons_gen_923

# SUBSET SAMPLE
gdf_sub = gdf.loc[gdf.in_sample].copy()



# %%
def format_aggpct_col(grouper, col, agg, params):
    summ = grouper[[col]].agg(agg) / params['denom']
    pct = (
        summ[col] / 
        summ.reset_index().groupby(cat, dropna=False)[col].transform('sum').values * 100 )
    out = pd.Series(
        summ[col].round(params['round']).astype(str).str.replace('nan', '-') + ' (' + 
        pct.round(params['round']).astype(str).str.replace('nan', '-') + '%)', name=f'{col}_n_pct')
    return out.to_frame()

def format_mstd_col(grouper, col, params):
    summ_m = grouper[[col]].mean() / params['denom']
    summ_s = grouper[[col]].std() / params['denom']
    out = pd.Series(
        summ_m[col].round(params['round']).astype(str).str.replace('nan', '-') + '±' + 
        summ_s[col].round(params['round']).astype(str).str.replace('nan', '-'), name=f'{col}_mean_std')
    return out.to_frame()

def format_iqr_col(grouper, col, params):
    summ_l = grouper[[col]].quantile(0.25) / params['denom']
    summ_h = grouper[[col]].quantile(0.75) / params['denom']
    out = pd.Series(
        '[' + summ_l[col].round(params['round']).astype(str).str.replace('nan', '') + '-' + 
        summ_h[col].round(params['round']).astype(str).str.replace('nan', '') + ']', name=f'{col}_iqr')
    return out.to_frame()
    


# %% 
# SUMMARIZE BY CATEGORY
cats = [['year'], ['year', 'energy_source_1_cat', 'energy_source_1_subcat'], ['year', 'nerc_region']]
cats.reverse()
# cats = [['year', 'energy_source_1_cat']]
summ = pd.DataFrame()

for cat in cats:
    grouper = gdf_sub.groupby(cat + ['em_cat'], dropna=False)
    summ_cat = pd.DataFrame()
    
    # get age
    cols = ['age']
    params = {'round':ROUND, 'denom':1}
    for col in cols:
        summ_cat = pd.concat([format_mstd_col(grouper, col, params), summ_cat], axis=1)
        summ_cat = pd.concat([format_iqr_col(grouper, col, params), summ_cat], axis=1)
    
    # get capacity
    cols = ['nameplate_capacity_mw']
    params = {'round':ROUND, 'denom':1}
    for col in cols:
        summ_cat = pd.concat([format_aggpct_col(grouper, col, 'sum', params), summ_cat], axis=1)
        summ_cat = pd.concat([format_mstd_col(grouper, col, params), summ_cat], axis=1)
        summ_cat = pd.concat([format_iqr_col(grouper, col, params), summ_cat], axis=1)

    # get emissions
    cols = ['co2_mass_tons_gen', 'co2_mass_tons_gen_923', ]
    params = {'round':ROUND, 'denom':1e6}
    for col in cols:
        summ_cat = pd.concat([format_aggpct_col(grouper, col, 'sum', params), summ_cat], axis=1)
        summ_cat = pd.concat([format_mstd_col(grouper, col, params), summ_cat], axis=1)
        summ_cat = pd.concat([format_iqr_col(grouper, col, params), summ_cat], axis=1)

    # get counts
    cols = ['plant_code', 'gid']
    params = {'round':ROUND, 'denom':1}
    for col in cols:
        summ_cat = pd.concat([format_aggpct_col(grouper, col, 'nunique', params), summ_cat], axis=1)

    summ = pd.concat([summ_cat.reset_index(), summ], axis=0)
summ = summ.fillna('total')
cols_leading = ['year', 'energy_source_1_cat', 'energy_source_1_subcat', 'nerc_region']
cols_tailing = [col for col in summ.columns if col not in cols_leading]
summ = summ[cols_leading + cols_tailing]
summ.to_csv(PATH_RESULTS + 'df_compare_emissions_eia_epa.csv', index=False)
# %%

# %%
