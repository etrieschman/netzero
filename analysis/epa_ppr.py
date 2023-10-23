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
    'transportation_and_storage_costs_usd_per_mwh', 'state', 'scenario'
]

# %%
# READIN DATA
def download_ipm_gendata(scenario):
    # readin data
    df = pd.read_excel(f'{PATH_DATA}{scenario}/{scenario} RPE File.xlsx', sheet_name='RPE Report-2 with States')
    df.columns = (df.columns
                .str.lower()
                .str.replace(' ', '_', regex=False)
                .str.replace('/', '_per_', regex=False)
                .str.replace('$', 'usd', regex=False)
                .str.replace('-', '', regex=False))
    # clean data
    df['scenario'] = scenario.lower().replace(' ', '_')
    return df

# CLEAN DATA
def clean_ipm_gendata(df, cols_keep):
    df = df[cols_keep].copy()
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
    # bucket capacity factor
    step = 0.05
    bucket_capfac_bins = np.arange(0, 1+step, step)
    bucket_capfac_labs = [f'{l:0.2f}-{h:0.2f}' for l, h in 
                        zip(bucket_capfac_bins[0:-1], bucket_capfac_bins[1:])]
    df['capacity_factor_bucket'] = pd.cut(df.capacity_factor, bins=bucket_capfac_bins, 
                                          labels=bucket_capfac_labs, right=True, include_lowest=True)
    return df

# MAKE SUMMARY TABLES
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
    dfs_all = pd.concat([dfs, dfs_total], axis=0, ignore_index=True)
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

# MAKE PLOT SUMMARY
def plot_summ(dfs:pd.DataFrame, groupcol:str, 
              agg_var:str, agg_fn:str, 
              fmt_scenario:str='style', logscale:bool=True, outfn=None):
    fmt_scenario = {fmt_scenario:'scenario'}
    s = dfs.groupby(['year', groupcol, 'scenario'])[agg_var].agg(agg_fn).reset_index()
    # plot
    sns.relplot(kind='line', data=s, x='year', y=agg_var, hue=groupcol, **fmt_scenario, legend=True)
    if logscale:
        plt.yscale('log')
    
    if outfn is not None:
        plt.savefig(f'{outfn}_{agg_var[:3]}.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()

def plot_all_summ(dfs, groupcol, agg, fmt_scenario, logscales, outfn=None):
    for (agg_var, agg_fn), ls in zip(agg.items(), logscales):
        plot_summ(dfs, groupcol, agg_var, agg_fn, fmt_scenario, ls, outfn)


# %%
if __name__ == '__main__':
    # GET DATASETS
    print('reading in data...')
    dfb = download_ipm_gendata(scen_base)
    dfbu = download_ipm_gendata(scen_base_updated)
    dfp = download_ipm_gendata(scen_prop)
    dfpu = download_ipm_gendata(scen_prop_updated)

    # CLEAN DATASETS
    print('cleaning data...')
    dfbc = clean_ipm_gendata(dfb.copy(), cols_keep)
    dfbuc = clean_ipm_gendata(dfbu.copy(), cols_keep)
    dfpc = clean_ipm_gendata(dfp.copy(), cols_keep)
    dfpuc = clean_ipm_gendata(dfpu.copy(), cols_keep)

    # MAKE ANALYSIS DATASETS
    dfs_all = pd.concat([dfbuc, dfpc, dfpuc], axis=0).reset_index()
    dfs_all = dfs_all.loc[dfs_all.fuel_type.notna() & dfs_all.capacity_reporting_type.notna()]
    dfs = dfs_all.loc[~dfs_all.in_canada]
    dfs_gas = dfs.loc[dfs.fuel_cat == 'natural_gas']
    dfs_canada = dfs_all.loc[dfs_all.in_canada]
    dfs_canada_gas = dfs_canada.loc[dfs_canada.fuel_cat == 'natural_gas']

    # SUMMARIZE
    print('summarizing data...')
    agg = {
        # 'dispatchable_capacity_gw':'sum',
        'co2_emissions_total_million_metric_tons':'sum',
        'generation_total_gwh':'sum',
        'capacity_factor':'mean',
        }
    s = make_table_summs(dfs, 'fuel_cat', agg, PATH_RESULTS + 'df_fuelcat')
    sgas = make_table_summs(dfs_gas, 'capacity_mw_bucket', agg, PATH_RESULTS + 'df_gascap')
    # canada analysis
    scan = make_table_summs(dfs_canada, 'fuel_cat', agg)
    scangas = make_table_summs(dfs_canada_gas, 'capacity_mw_bucket', agg)

    # PLOT
    logscales = [False, True, False]
    assert(len(logscales) == len(agg.keys()))
    plot_all_summ(dfs, 'fuel_cat', agg, 'style', logscales, PATH_RESULTS + 'fig_fuelcat')
    plot_all_summ(dfs_gas, 'capacity_mw_bucket', agg, 'style', logscales, PATH_RESULTS + 'fig_gascap')
    plot_all_summ(dfs_canada, 'fuel_cat', agg, 'style', logscales)
    plot_all_summ(dfs_canada_gas, 'capacity_mw_bucket', agg, 'style', logscales)


# %%
# ANALYSIS: quantify CCS
dfs['ccs'] = dfs.capacity_reporting_type.str.lower().str.contains('ccs')
dfs_gascoal = dfs.loc[
    dfs.scenario.isin(['updated_baseline', 'integrated_proposal']) & 
    dfs.fuel_cat.isin(['coal', 'natural_gas'])]
agg = {'unitid':'count',
    'generation_total_gwh':'sum'}

# full summary dataset
display(dfs_gascoal.groupby(['scenario', 'fuel_cat', 'ccs', 'year']).agg(agg))


# PLOT
summ = (dfs_gascoal
        .loc[dfs_gascoal.ccs]
        .groupby(['scenario', 'year'])
        .agg(agg).reset_index())


fig, ax = plt.subplots()
sns.barplot(data=summ, x='year', y='generation_total_gwh', 
            hue='scenario', ax=ax)
ax.set_ylabel('Generation (GWh)')
ax.set_title('Generation from NG and coal with CCS')
plt.show()

fig, ax = plt.subplots()
sns.barplot(data=summ, x='year', y='unitid', 
            hue='scenario', ax=ax)
ax.set_ylabel('Units')
ax.set_title('NG and coal units with CCS')
plt.show()




# %%
# ANALYSIS: shift in capacity factor
summ = (dfs_gas
        .loc[dfs_gas.scenario == 'integrated_proposal']
        .groupby(['year', 'capacity_mw_bucket', 'capacity_factor_bucket'])
        .agg({'unitid':'count',
              'generation_total_gwh':'sum',
              'co2_emissions_total_million_metric_tons':'sum'})
        .reset_index())

summ_w = summ.pivot(index=['capacity_mw_bucket', 'year'],
            columns=['capacity_factor_bucket'],
            values=['co2_emissions_total_million_metric_tons'])
summ_w.to_csv(PATH_RESULTS + 'cf_heatmap.csv')


# %%
# ANALYSIS: read-out for Hannah Dobie
# MOTIVATING QUESTION: What state is getting hit the hardest?
dfs['regulated_gas_cap'] = False
dfs.loc[(dfs.fuel_cat == 'natural_gas') & (dfs.dispatchable_capacity_mw >= 300), 'regulated_gas_cap'] = True
dfs['regulated_gas_cf'] = False
dfs.loc[(dfs.fuel_cat == 'natural_gas') & (dfs.capacity_factor >= 0.5), 'regulated_gas_cf'] = True
dfs['regulated_gas'] = dfs.regulated_gas_cap & dfs.regulated_gas_cf
dfs['regulated_coal'] = dfs.fuel_cat == 'coal'
dfs['regulated'] = dfs.regulated_gas | dfs.regulated_coal

state_readout = (dfs.loc[dfs.fuel_cat.isin(['natural_gas', 'coal'])]
                 .groupby(['scenario', 'year', 'state', 'fuel_cat', 
                           'regulated', 'regulated_gas_cap'], dropna=False)
                 .agg({'unitid':'count',
                        'co2_emissions_total_million_metric_tons':['sum'], 
                        'generation_total_gwh':['sum']})
                 .reset_index()
)
state_readout.to_csv(PATH_RESULTS + 'temp_statereadout.csv', index=False)
# %%

(dfbuc.loc[(dfbuc.year == 2028) ]
 .groupby(['fuel_cat']).agg({'unitid':'count'}))
# %%
