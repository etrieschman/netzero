# %%
# SETUP CODE
# Data downloaded from:           https://www.epa.gov/power-sector-modeling/analysis-proposed-greenhouse-gas-standards-and-guidelines
# Data documentation available:   https://www.epa.gov/system/files/documents/2023-04/IPM%20InputOutputGuide_2023.pdf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 70)

PATH_DATA = '../data/raw/ipm/'
PATH_RESULTS = '../results/epa_ppr/'
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
    df = pd.read_excel(f'{PATH_DATA}{scenario}/{scenario} RPE File.xlsx', sheet_name='RPE Report-2 with States')
    df.columns = (df.columns
                .str.lower()
                .str.replace(' ', '_').str.replace('/', '_per_')
                .str.replace('$', 'usd').str.replace('-', ''))
    # clean data
    df = df[cols_keep].copy()
    df['scenario'] = scenario.lower().replace(' ', '_')
    return df

def clean_ipm_gendata(df):
    df['dispatchable_capacity_gw'] = df.dispatchable_capacity_mw / 1000
    df['co2_emissions_total_million_metric_tons'] = df.co2_emissions_total_thousand_metric_tons / 1000
    df['co2_emissions_per_generation_total_tmt_per_gwh'] = df.co2_emissions_total_thousand_metric_tons / df.generation_total_gwh
    df['capacity_factor'] = (df.generation_total_gwh) / (df.dispatchable_capacity_gw*8760)
    df['in_canada'] = df.region_group_acronym.str.startswith('CN')
    # categorize fuel types
    df.loc[df.fuel_type.isin(
        ['Wind', 'Solar', 'Solar-DER', 'Geothermal', 'EnerStor']), 'fuel_cat'] = 'renewables'
    df.loc[df.fuel_type.isin(['Nuclear']), 'fuel_cat'] = 'nuclear'
    df.loc[df.fuel_type.isin(['Hydro', 'Pumps']), 'fuel_cat'] = 'hydro'
    df.loc[df.fuel_type.isin(
        ['Coal', 'Waste Coal', 'Pet. Coke']), 'fuel_cat'] = 'coal'
    df.loc[df.fuel_type.isin(['NaturalGas']), 'fuel_cat'] = 'natural_gas'
    df.loc[df.fuel_type.isin(['Oil']), 'fuel_cat'] = 'oil'
    df.loc[df.fuel_type.isin(['Hydrogen']), 'fuel_cat'] = 'hydrogen'
    df.loc[df.fuel_type.isin(
        ['LF Gas', 'Non-Fossil', 'MSW', 'Fwaste', 'Biomass', 'Tires']), 'fuel_cat'] = 'nonfossil_nonren'
    # bucket capacity
    bucket_cap = {-1:'0000-0149', 149:'0150-0299', 299:'0300-0599', 599:'0600-1199', 1199:'1200-2399', 
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
# ANALYSIS DATASETS
dfs_all = pd.concat([dfbuc, dfpc, dfpuc], axis=0).reset_index()
dfs_all = dfs_all.loc[dfs_all.fuel_type.notna() & dfs_all.capacity_reporting_type.notna()]
dfs = dfs_all.loc[~dfs_all.in_canada]
dfs_gas = dfs.loc[dfs.fuel_cat == 'natural_gas']
dfs_canada = dfs_all.loc[dfs_all.in_canada]
dfs_canada_gas = dfs_canada.loc[dfs_canada.fuel_cat == 'natural_gas']

# %%
# MAKE SUMMARY TABLE
def make_pctchg_table(sdf, groups, agg, varbase, valbase):
    sdf_pct = sdf.copy()
    for k in agg.keys():
        baseline_mask = sdf_pct[varbase] == valbase
        sdf_pct['baseline'] = np.nan
        sdf_pct.loc[baseline_mask, 'baseline'] = sdf_pct.loc[baseline_mask, k]
        sdf_pct['baseline'] = sdf_pct.groupby(groups)['baseline'].transform('max')
        sdf_pct[k] = sdf_pct[k] / sdf_pct.baseline - 1
        sdf_pct.drop(columns='baseline', inplace=True)
    return sdf_pct

def make_table_summs(dfs, group, agg, outfn=None):
    # summary table
    dfs_total = dfs.copy()
    dfs_total[group] = 'total'
    dfs_all = pd.concat([dfs, dfs_total], axis=0)
    s = (dfs_all
        .groupby([group, 'year', 'scenario'], dropna=False)
        .agg(agg)
        .reset_index()
        )
    # get results as a percent of baseline
    spct_scen = make_pctchg_table(s, [group, 'year'], agg, 'scenario', 'updated_baseline')
    spct_scen = spct_scen.pivot(index=[group, 'year'], columns='scenario')
    spct_scen.drop(columns=[col for col in spct_scen.columns if 'updated_baseline' in col], inplace=True)
    # get results as a percent of 2028
    spct_yr = make_pctchg_table(s, [group, 'scenario'], agg, 'year', 2028)
    spct_yr = spct_yr.pivot(index=[group, 'year'], columns='scenario')
    # pivot
    s = s.pivot(index=[group, 'year'], columns='scenario')
    if outfn is not None:
        s.to_csv(outfn+'.csv')
        spct_scen.to_csv(outfn + '_pctbase.csv')
        spct_yr.to_csv(outfn + '_pct2028.csv')
    return {'value':s, 'pct2028':spct_yr, 'pctbase':spct_scen}

# %%
# SUMMARIZE
agg = {
    # 'dispatchable_capacity_gw':'sum',
    'co2_emissions_total_million_metric_tons':'sum',
    'generation_total_gwh':'sum',
    'capacity_factor':'mean',
    }
s = make_table_summs(dfs, 'fuel_cat', agg, PATH_RESULTS + 'df_fuelcat')
sgas = make_table_summs(dfs_gas, 'capacity_mw_bucket', agg, PATH_RESULTS + 'df_gascap')



# %%
# MAKE PLOT SUMMARY
def plot_summ(dfs:pd.DataFrame, groupcol:str, 
              agg_var:str, agg_fn:str, 
              fmt_scenario:str='style', logscale:bool=True):
    fmt_scenario = {fmt_scenario:'scenario'}
    s = dfs.groupby(['year', groupcol, 'scenario'])[agg_var].agg(agg_fn).reset_index()
    # plot
    sns.relplot(kind='line', data=s, x='year', y=agg_var, hue=groupcol, **fmt_scenario, legend=True)
    if logscale:
        plt.yscale('log')
    plt.show()

# plot fuel category
logscales = [False, True, False]
for (agg_var, agg_fn), logscale in zip(agg.items(), logscales):
    plot_summ(dfs=dfs, groupcol='fuel_cat', 
              agg_var=agg_var, agg_fn=agg_fn, fmt_scenario='style', logscale=logscale)
    

# plot natgas capacity
logscales = [False, True, False]
for (agg_var, agg_fn), logscale in zip(agg.items(), logscales):
    plot_summ(dfs=dfs_gas, groupcol='capacity_mw_bucket', 
              agg_var=agg_var, agg_fn=agg_fn, fmt_scenario='style', logscale=logscale)
    
# %%
# CANADA ANALYSIS
scan = make_table_summs(dfs_canada, 'fuel_cat', agg)
scangas = make_table_summs(dfs_canada_gas, 'capacity_mw_bucket', agg)

# %%
logscales = [True, False, True, True]
for (agg_var, agg_fn), logscale in zip(agg.items(), logscales):
    plot_summ(dfs=dfs_canada, groupcol='fuel_cat', 
              agg_var=agg_var, agg_fn=agg_fn, fmt_scenario='style', logscale=logscale)
# plot natgas capacity
logscales = [True, False, False, False]
for (agg_var, agg_fn), logscale in zip(agg.items(), logscales):
    plot_summ(dfs=dfs_canada_gas, groupcol='capacity_mw_bucket', 
              agg_var=agg_var, agg_fn=agg_fn, fmt_scenario='style', logscale=logscale)

# %%
