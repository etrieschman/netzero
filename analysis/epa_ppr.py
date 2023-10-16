# %%
# SETUP CODE
# Data downloaded from:           https://www.epa.gov/power-sector-modeling/analysis-proposed-greenhouse-gas-standards-and-guidelines
# Data documentation available:   https://www.epa.gov/system/files/documents/2023-04/IPM%20InputOutputGuide_2023.pdf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

pd.set_option('display.max_columns', None)

PATH_DATA = '../data/epa_ppr/'
scen_base = 'Post-IRA 2022 Reference Case'
scen_base_updated = 'Updated Baseline'
scen_prop = 'Proposal'
scen_prop_updated = 'Integrated Proposal'

cols_keep = [
    'year', 'region_group_acronym', 'unitid', 'unit_long_name',
    'prior_stage', 'base_unit', 'fuel_type',
    'co2_content_lbs_per_mmbtu', 'hcl_content_lbs_per_mmbtu',
    'mer_content_lbs_per_mmbtu', 'so2_content_lbs_per_mmbtu',
    'fuel_consumption_total_tbtu', 'generation_total_gwh',
    'co2_emissions_total_thousand_metric_tons',
    'co2_emissions_total_thousand_tons',
    'average_capacity_mw',
    'dispatchable_capacity_mw', 'capacity_reporting_type',
    'fom_usd_per_kwyear', 'vom_usd_per_mwh', 'fuel_usd_per_mmbtu',
    'capital_cost_usd_per_kwyear',
    'transportation_and_storage_costs_usd_per_mwh'
]

# %%
# READIN DATA
def download_ipm_gendata(scenario, cols_keep):
    # readin data
    df = pd.read_excel(f'{PATH_DATA}{scenario}/{scenario} RPE File.xlsx', sheet_name='RPE Report-1')
    df.columns = (df.columns
                .str.lower()
                .str.replace(' ', '_').str.replace('/', '_per_')
                .str.replace('$', 'usd').str.replace('-', ''))
    # clean data
    df = df[cols_keep].copy()
    df['scenario'] = scenario.lower().replace(' ', '')
    return df

def clean_ipm_gendata(df):
    df['dispatchable_capacity_gw'] = df.dispatchable_capacity_mw / 1000
    df['capacity_factor'] = (df.generation_total_gwh) / (df.dispatchable_capacity_gw*8760)
    # categorize fuel types
    df.loc[df.fuel_type.isin(
        ['Wind', 'Solar', 'Solar-DER', 'Geothermal', 'EnerStor']), 'fuel_cat'] = 'renewables'
    df.loc[df.fuel_type.isin(['Nuclear']), 'fuel_cat'] = 'nuclear'
    df.loc[df.fuel_type.isin(['Hydro', 'Pumps']), 'fuel_cat'] = 'hydro'
    df.loc[df.fuel_type.isin(
        ['Coal', 'Waste Coal', 'Pet. Coke']), 'fuel_cat'] = 'coal'
    df.loc[df.fuel_type.isin(['NaturalGas']), 'fuel_cat'] = 'natural_gas'
    df.loc[df.fuel_type.isin(['Oil']), 'fuel_cat'] = 'oil'
    df.loc[df.fuel_type.isin(
        ['LF Gas', 'Non-Fossil', 'MSW', 'Fwaste', 'Biomass', 'Tires']), 'fuel_cat'] = 'nonfossil_nonren'
    df.loc[df.fuel_cat.isna() & 
            df.capacity_reporting_type.str.contains('Retirement'), 'fuel_cat'] = 'retire'
    df.loc[df.fuel_cat.isna(), ['fuel_cat', 'capacity_reporting_type']].value_counts(dropna=False)
    # bucket capacity
    bucket_cap = {0:'0000-0149', 149:'0150-0299', 299:'0300-0599', 599:'0600-1199', 1199:'1200-2399', 
                2399:'2400-3599', 3599:'3600+', 1e4:'>10k'}
    df['capacity_mw_bucket'] = pd.cut(df.dispatchable_capacity_mw, bins=list(bucket_cap.keys()), 
                                 labels=list(bucket_cap.values())[:-1], right=True, include_lowest=True).astype(str)
    df.loc[df.capacity_reporting_type.str.contains('New'), 'capacity_mw_bucket'] = 'New'
    return df



# %%
# GET DATASETS
dfb = download_ipm_gendata(scen_base, cols_keep)
dfbu = download_ipm_gendata(scen_base_updated, cols_keep)
dfp = download_ipm_gendata(scen_prop, cols_keep)
dfpu = download_ipm_gendata(scen_prop_updated, cols_keep)

# %%
# CLEAN DATASETS
dfbc = clean_ipm_gendata(dfb.copy())
dfbuc = clean_ipm_gendata(dfbu.copy())
dfpc = clean_ipm_gendata(dfp.copy())
dfpuc = clean_ipm_gendata(dfpu.copy())

# %%
# PLOT SUMMARY
fig, ax = plt.subplots(ncols=2, sharey=True)

groups = ['year', 'fuel_cat']
# agg_var, agg_fn = 'dispatchable_capacity_gw', 'sum'
# agg_var, agg_fn = 'capacity_factor', 'mean'
agg_var, agg_fn = 'generation_total_gwh', 'sum'
agg_var, agg_fn = 'co2_emissions_total_thousand_tons', 'sum'
dfs = [dfbuc, dfpuc]
df_titles = ['baseline', 'integrated proposal']
for i, (df, df_title) in enumerate(zip(dfs, df_titles)):
    # df = df.loc[df.fuel_cat == 'natural_gas']
    s = df.groupby(groups)[agg_var].agg(agg_fn).reset_index()
    legend = (i+1 == len(dfs))
    sns.lineplot(data=s, x='year', y=agg_var, hue=groups[-1], ax=ax[i], legend=legend)
    ax[i].set_title(df_title)
    if legend:
        ax[i].legend(bbox_to_anchor=(1,1), loc='upper left', title=groups[-1])
plt.suptitle(f'Electricity production from natural gas ({agg_var})')

# %%