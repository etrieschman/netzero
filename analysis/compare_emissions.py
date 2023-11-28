# %%
# SETUP
import pandas as pd
import numpy as np

# global variables
PATH_PROCESSED = '../data/processed/'
ROUND, DENOM = 3, 1e3

# %%
# READ IN DATA
g = pd.read_parquet(PATH_PROCESSED + 'df_generators.parquet')
df = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
gdf = pd.merge(left=g.drop(columns='gid'), right=df, how='left', on=['year', 'plant_code', 'generator_id'])
assert len(g) == len(gdf)
gdf.loc[:,'em_cat'] = 'both'
gdf[['has_emissions_epa', 'has_emissions_eia']] = gdf[['has_emissions_epa', 'has_emissions_eia']].astype(bool)
gdf[['co2_mass_tons_gen', 'co2_mass_tons_gen_923']] = gdf[['co2_mass_tons_gen', 'co2_mass_tons_gen_923']].astype(float)
gdf.loc[gdf.has_emissions_epa & ~gdf.has_emissions_eia, 'em_cat'] = 'epa_only'
gdf.loc[~gdf.has_emissions_epa & gdf.has_emissions_eia, 'em_cat'] = 'eia_only'
gdf.loc[:,'co2_mass_diff'] = gdf.co2_mass_tons_gen - gdf.co2_mass_tons_gen_923

# %%
gdf.has_emissions_epa.mean()


# %% 
# SUMMARIZE BY CATEGORY
cats = [['year'], ['year', 'energy_source_1_cat']]

summ_cats = pd.DataFrame()

for cat in cats:
    grouper = gdf.groupby(cat + ['em_cat'], dropna=False)

    # get counts
    summ_n = grouper['plant_code', 'gid'].nunique()
    for col in summ_n.columns:
        summ_n.loc[:, f'{col}_pct'] = (
            summ_n[col].values / 
            summ_n.reset_index().groupby(cat)[col].transform('sum').values * 100).round(ROUND).astype(str)
        summ_n[col] = summ_n[col].astype(str) + ' (' + summ_n[f'{col}_pct'] + '%)'
        summ_n.drop(columns=f'{col}_pct', inplace=True)
    summ_n

    # get totals
    co2s = ['co2_mass_tons_gen', 'co2_mass_tons_gen_923']
    summ_em = grouper[co2s].sum() / DENOM
    for col in summ_em.columns:
        summ_em.loc[:, f'{col}_pct'] = (
            summ_em[col].values / 
            summ_em.reset_index().groupby(cat)[col].transform('sum').values * 100).round(ROUND).astype(str)
        summ_em[col] = summ_em[col].round(ROUND).astype(str) + ' (' + summ_em[f'{col}_pct'] + '%)'
        summ_em.drop(columns=f'{col}_pct', inplace=True)
    summ_em

    # get means
    co2s = ['co2_mass_tons_gen', 'co2_mass_tons_gen_923']
    summ_m = grouper[co2s].mean() / DENOM
    summ_s = grouper[co2s].std() / DENOM
    for col in summ_m.columns:
        summ_m[f'{col}_ms'] = summ_m[col].round(ROUND).astype(str) + '±' + summ_s[col].round(ROUND).astype(str)
        summ_m.drop(columns=col, inplace=True)
    summ_m

    # get range
    co2s = ['co2_mass_tons_gen', 'co2_mass_tons_gen_923']
    summ_l = grouper[co2s].quantile(0.25) / DENOM
    summ_h = grouper[co2s].quantile(0.75) / DENOM
    for col in summ_l.columns:
        summ_l[f'{col}_range'] = '[' + summ_l[col].round(ROUND).astype(str) + '-' + summ_h[col].round(ROUND).astype(str) + ']'
        summ_l.drop(columns=col, inplace=True)
    summ_l

    summ = pd.concat([summ_n, summ_em, summ_m, summ_l], axis=1)
    summ_cats = pd.concat([summ.reset_index(), summ_cats])

summ_cats


# %%
