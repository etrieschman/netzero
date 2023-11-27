# %%
# SETUP
import pandas as pd
import numpy as np

# global variables
PATH_PROCESSED = '../data/processed/'
denom = 1e6

# %%
# READ IN DATA
df = pd.read_parquet('../data/processed/df_emissions.parquet')
df.loc[:,'em_cat'] = 'both'
df[['has_emissions_epa', 'has_emissions_eia']] = df[['has_emissions_epa', 'has_emissions_eia']].astype(bool)
df[['co2_mass_short_tons_gen', 'co2_mass_short_tons_gen_923']] = df[['co2_mass_short_tons_gen', 'co2_mass_short_tons_gen_923']].astype(float)
df.loc[df.has_emissions_epa & ~df.has_emissions_eia, 'em_cat'] = 'epa_only'
df.loc[~df.has_emissions_epa & df.has_emissions_eia, 'em_cat'] = 'eia_only'
df.loc[:,'co2_mass_diff'] = df.co2_mass_short_tons_gen - df.co2_mass_short_tons_gen_923


# %%
# SUMMARIZE BY CATEGORY
cats = [['year']]
cat = cats[0]

summ_cats = pd.DataFrame()

for cat in cats:
    grouper = df.groupby(cat + ['em_cat'], dropna=False)

    # get counts
    summ_n = grouper['plant_code', 'gid'].nunique()
    for col in summ_n.columns:
        summ_n.loc[:, f'{col}_pct'] = (
            summ_n[col].values / 
            summ_n.reset_index().groupby(cat)[col].transform('sum').values * 100).round(2).astype(str)
        summ_n[col] = summ_n[col].astype(str) + ' (' + summ_n[f'{col}_pct'] + '%)'
        summ_n.drop(columns=f'{col}_pct', inplace=True)
    summ_n

    # get totals
    co2s = ['co2_mass_short_tons_gen', 'co2_mass_short_tons_gen_923']
    summ_em = grouper[co2s].sum() / denom
    for col in summ_em.columns:
        summ_em.loc[:, f'{col}_pct'] = (
            summ_em[col].values / 
            summ_em.reset_index().groupby(cat)[col].transform('sum').values * 100).round(2).astype(str)
        summ_em[col] = summ_em[col].round(2).astype(str) + ' (' + summ_em[f'{col}_pct'] + '%)'
        summ_em.drop(columns=f'{col}_pct', inplace=True)
    summ_em

    # get means
    co2s = ['co2_mass_short_tons_gen', 'co2_mass_short_tons_gen_923', 'co2_mass_diff']
    summ_m = grouper[co2s].mean() / denom
    summ_s = grouper[co2s].std() / denom
    for col in summ_m.columns:
        summ_m[f'{col}_ms'] = summ_m[col].round(2).astype(str) + '±' + summ_s[col].round(2).astype(str)
        summ_m.drop(columns=col, inplace=True)
    summ_m

    # get range
    co2s = ['co2_mass_short_tons_gen', 'co2_mass_short_tons_gen_923', 'co2_mass_diff']
    summ_l = grouper[co2s].agg(lambda x: np.percentile(x, 25)) / denom
    summ_h = grouper[co2s].agg(lambda x: np.percentile(x, 75)) / denom
    for col in summ_l.columns:
        summ_l[f'{col}_range'] = '[' + summ_l[col].round(2).astype(str) + '-' + summ_h[col].round(2).astype(str) + ']'
        summ_l.drop(columns=col, inplace=True)
    summ_l

    summ = pd.concat([summ_n, summ_em, summ_m, summ_l], axis=1)
    summ_cats = pd.concat([summ, summ_cats])


# %%
summ_cats
# %%

