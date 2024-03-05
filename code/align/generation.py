# %%
# SET UP
import numpy as np
import pandas as pd
from pathlib import Path

# global variables
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 150)

# %%
# HELPER FUNCTIONS
def get_annual_counts(df, groupvars='year', countvars='utility_id'):
    return df.groupby(groupvars)[countvars].count()

def subset_from_ids(df, fromdf, idcol):
    return df.loc[df[idcol].isin(fromdf[idcol].drop_duplicates().values)]

def get_annual_generation(df, groupvars='year', sumvars='net_gen_tot_an', 
                          denom=DENOM, round=ROUND):
    return (df.groupby(groupvars)[sumvars].sum() / denom).round(round)


# %%
# CLEAN DATASETS
def clean_eia923(opsdf):
    # drop generation columns to keep things cleaner
    dropvar_keys = ('_jan', '_feb', '_mar', '_apr', '_may', '_jun', '_jul', 
                    '_aug', '_sep', '_oct', '_nov', '_dec',
                    'file', '_per_unit', 'balancing_authority_code', 'aer_fuel_type_code',
                    'physical_unit_label', 'plant_name', 'operator_name',
                    'plant_state', 'census_region', 'nerc_region','naics_code',
                    'sector_number','sector_name', 'respondent_frequency')
    dropvar = opsdf.columns[opsdf.columns.str.contains('|'.join(dropvar_keys))]
    # rename key variables
    opsdf = (opsdf.drop(columns=dropvar)
            .rename(columns={'plant_id':'plant_code', 'operator_id':'utility_id'})
            .astype({'utility_id':'Int64'}))
    # some generator ids have leading 0s
    opsdf['generator_id'] = opsdf.generator_id.str.lstrip('0')

    # summarize
    print('Total generation by aggregation and year (TWh):')
    print(opsdf.groupby(['year', 'sheet'])['net_gen_tot_an'].sum())

    # split df into generator and plant part data
    opgendf = (opsdf
              .loc[opsdf.sheet == 'page_4_generator_data']
              .drop(columns=['reported_fuel_type_code', 'sheet', 'quantity_tot_an', 
                             'elec_quantity_tot_an','tot_mmbtu_tot_an', 'elec_mmbtu_tot_an']))
    opppdf = (opsdf
             .loc[opsdf.sheet == 'page_1_generation_and_fuel_data']
             .drop(columns=['generator_id','sheet']))
    summ = {}
    summ['n_id_generators'] = get_annual_counts(opgendf)
    summ['n_id_plantpart'] = get_annual_counts(opppdf)
    summ['twh_generators'] = get_annual_generation(opgendf)
    summ['twh_plantpart'] = get_annual_generation(opppdf)

    # collapse on primary unit ids (across combined heat and power and nuclear unit id)
    opgendf = (opgendf.groupby(['year', 'utility_id', 'plant_code', 'generator_id'], dropna=False)
               ['net_gen_tot_an'].sum().reset_index())
    opppdf = (opppdf.groupby(['year', 'utility_id', 'plant_code', 
                              'reported_prime_mover', 'reported_fuel_type_code'], dropna=False)
              [opppdf.columns[opppdf.columns.str.endswith(('_an'))]].sum().reset_index())
    summ['n_id_generators_postcol'] = get_annual_counts(opgendf)
    summ['n_id_plantpart_postcol'] = get_annual_counts(opppdf)
    summ = pd.DataFrame(summ).T

    # return
    return opgendf, opppdf, summ



# %% 
# ALIGN DATASETS
def align_eia923_to_eia860(gdf, opsgendf, opsppdf):
    # First, merge generator generation data onto EIA860 generator data
    gg = pd.merge(
        right=opsgendf, how='outer',
        left=gdf,
        on=['year', 'utility_id', 'plant_code', 'generator_id'],
        suffixes=('_g', '_go'), indicator='merge_go')

    # Second, merge plant part generation data onto EIA860 generator data
    pp_on = ['year', 'utility_id', 'plant_code', 'reported_prime_mover', 'reported_fuel_type_code']
    gg_on = ['year','utility_id', 'plant_code', 'prime_mover', 'energy_source_1']
    ggp = pd.merge(
        left=gg, how='outer',
        right=opsppdf,
        left_on=gg_on,
        right_on=pp_on,
        suffixes=('', '_po'), indicator='merge_po')

    # summarize
    summ_g = gg.groupby(['year', 'merge_go']).agg({'utility_id':'count', 'net_gen_tot_an':'sum'})
    summ_g.columns = ['n', 'twh']
    summ_g.twh /= DENOM
    summ_g['cat'] = 'egu'
    summ_pp = (ggp[gg_on + ['net_gen_tot_an_po', 'merge_po']].drop_duplicates()
               .groupby(['year', 'merge_po']).agg({'utility_id':'count', 'net_gen_tot_an_po':'sum'}))
    summ_pp.columns = ['n', 'twh']
    summ_pp.twh /= DENOM
    summ_pp['cat'] = 'plant_part'
    summ = pd.concat([summ_g, summ_pp], names=['year', 'merge'], axis=0)

    # drop right-only merges
    ggp = ggp.loc[(ggp.merge_po != 'right_only') & (ggp.merge_go != 'right_only')]

    # return
    return ggp, summ


# %%
# ALLOCATE OPERATIONS TO GENERATORS
def allocate_eia923_operations(ggp):
    # flag source of generation data (gen or plant plart)
    groupvars = ['year', 'utility_id', 'plant_code', 'prime_mover', 'energy_source_1']
    ggp['flag_gen_from_go'] = (ggp.merge_go == 'both')
    ggp['net_gen_tot_an_go'] = ggp.groupby(groupvars)['net_gen_tot_an'].transform('sum')
    ggp['0_gen_from_egu'] = ggp.flag_gen_from_go
    
    # get plant-part remaining generation, and summarize
    ggp['net_gen_tot_an_po_remaining'] = (ggp.net_gen_tot_an_po.fillna(0) - ggp.net_gen_tot_an_go.fillna(0)).clip(0)
         
    # which generators should we allocate to?
    print('Generation by status among EGUs with generation, 2013/2021 (N, MWh):')
    print(ggp.loc[(ggp.year.isin([2013, 2021])) & (ggp.net_gen_tot_an > 0)].groupby(['year', 'status']).agg({'utility_id':'count', 'net_gen_tot_an':lambda x:(x.sum()/DENOM).round(ROUND)}))
    # DECISION: allocate only to OP generators if available. The only to OA, OS, and SB if available. Then to All others if not
    ggp['group_has_opstatus'] = ggp.status == 'OP'
    ggp['group_has_opstatus'] = ggp.groupby(groupvars)['group_has_opstatus'].transform('max')
    ggp['group_has_otheractivestatus'] = ggp.status.isin(['OA', 'OS', 'SB'])
    ggp['group_has_otheractivestatus'] = ggp.groupby(groupvars)['group_has_otheractivestatus'].transform('max')
    
    # allocate to group with operation status
    ggp.loc[ggp.flag_gen_from_go, 'capacity_denominator'] = ggp.loc[ggp.flag_gen_from_go, 'nameplate_capacity_mw']
    mask_ops = (~ggp.flag_gen_from_go) & (ggp.status == 'OP')
    ggp.loc[mask_ops, 'capacity_denominator'] = ggp.loc[mask_ops].groupby(groupvars)['nameplate_capacity_mw'].transform('sum')
    ggp['1_plus_gen_from_opstatus'] = ggp.flag_gen_from_go | mask_ops
       
    # allocate to group with other active status
    mask_otops = (~ggp.flag_gen_from_go) & (ggp.status.isin(['OA', 'OS', 'SB'])) & (~ggp.group_has_opstatus)
    ggp.loc[mask_otops, 'capacity_denominator'] = ggp.loc[mask_otops].groupby(groupvars)['nameplate_capacity_mw'].transform('sum')
    ggp['2_plus_gen_from_otopstatus'] = ggp.flag_gen_from_go | mask_ops | mask_otops
          
    # allocate to groups with no operation or active status
    mask_res = (~ggp.flag_gen_from_go) & (~ggp.group_has_opstatus) & (~ggp.group_has_otheractivestatus)
    ggp.loc[mask_res, 'capacity_denominator'] = ggp.loc[mask_res].groupby(groupvars)['nameplate_capacity_mw'].transform('sum')
    ggp['3_plus_gen_from_reststatus'] = ggp.flag_gen_from_go | mask_ops | mask_otops | mask_res
         
    # final calculation
    ggp['pct_generation'] = ggp.nameplate_capacity_mw / ggp.capacity_denominator
    ggp['net_gen_tot_an_go'] = ggp.net_gen_tot_an.copy()
    ggp.loc[~ggp.flag_gen_from_go, 'net_gen_tot_an'] = ggp.pct_generation.astype(float) * ggp.net_gen_tot_an_po_remaining
    
    # allocate other operational data
    ggp['pct_allocation'] = ggp.net_gen_tot_an / ggp.groupby(groupvars)['net_gen_tot_an'].transform('sum')
    num_vars = ['quantity_tot_an', 'elec_quantity_tot_an', 'tot_mmbtu_tot_an', 'elec_mmbtu_tot_an']
    for var in num_vars:
        ggp[f'{var}_po'] = ggp[var].copy()
        ggp[var] *= ggp['pct_allocation']

    # summary
    print('Allocation by status, 2013/2021 (N, TWh):\n', 
          (ggp.loc[(ggp.year.isin([2013, 2021]))]
           .groupby(['year', 'status'])
           .agg({'utility_id':'count', 'net_gen_tot_an':lambda x:(x.sum()/DENOM).round(ROUND)})))

    # summarize
    steps = ['0_gen_from_egu', '1_plus_gen_from_opstatus', '2_plus_gen_from_otopstatus', '3_plus_gen_from_reststatus']
    summ = pd.DataFrame()
    for step in steps:
        summ_g = (ggp[groupvars + ['generator_id', 'net_gen_tot_an', step]]
                  .drop_duplicates()
                  .groupby([step, 'year'])
                  .agg({'utility_id':'count', 'net_gen_tot_an':'sum'}))
        summ_pp = (ggp[groupvars + ['net_gen_tot_an_po_remaining', step]]
                  .drop_duplicates()
                  .groupby([step, 'year'])
                  .agg({'net_gen_tot_an_po_remaining':'sum'}))
        summ_s = pd.merge(summ_g, summ_pp, left_index=True, right_index=True)
        summ_s.columns = ['n', 'twh_tot', 'twh_pp_remaining']
        summ_s.twh_tot /= DENOM
        summ_s.twh_pp_remaining /= DENOM
        summ_s['step'] = step
        summ = pd.concat([summ_s, summ], axis=0)
    summ = summ.pivot(columns='step')
    summ.index.names = ['allocated', 'year']
    summ.loc[summ.index.get_level_values(0), [col for col in summ.columns if 'remaining' in col[0]]] = pd.NA
    summ.loc[summ.index.get_level_values(0) == False, [col for col in summ.columns if 'remaining' not in col[0]]] = pd.NA
    summ = summ.groupby('year').sum()

    return ggp, summ



# %%
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('align_generation')
    
    results_dir = Path(snakemake.params.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    final_csv_dir = Path(snakemake.params.final_csv_dir)
    final_csv_dir.mkdir(parents=True, exist_ok=True)
    
    # load data
    gdf = pd.read_parquet(snakemake.input.infile_gen)
    opsdf = pd.read_parquet(snakemake.input.infile_ops)
    
    # clean datasets
    opsgendf, opsppdf, opssumm = clean_eia923(opsdf)
    opssumm = pd.DataFrame(opssumm).T
    opssumm.to_csv(results_dir / 'df_summ_generation.csv')

    # align datasets
    ggp, alignsumm = align_eia923_to_eia860(gdf, opsgendf, opsppdf)
    alignsumm.to_csv(results_dir / 'df_summ_align.csv')
    
    # allocate operations data
    ggp, allsumm = allocate_eia923_operations(ggp)
    allsumm.to_csv(results_dir / 'df_summ_allocate.csv')

    # calculate heatrate
    ggp['tot_hr_mmbtu_per_unit_an'] = np.where(
        ggp.quantity_tot_an.notna() & (ggp.quantity_tot_an != 0),
        ggp.tot_mmbtu_tot_an / ggp.quantity_tot_an, np.nan)

    # write to file
    vars_keep = ['year', 'utility_id', 'plant_code', 'generator_id', 
                #  'nameplate_capacity_mw', 'status', 
                #  'energy_source_1', 'energy_source_1_subcat', 'energy_source_1_cat',
                # 'reported_prime_mover', 'reported_fuel_type_code',
                'net_gen_tot_an', 'quantity_tot_an', 'elec_quantity_tot_an', 
                'tot_mmbtu_tot_an', 'elec_mmbtu_tot_an', 'tot_hr_mmbtu_per_unit_an',
                'pct_generation', 'pct_allocation']
    vars_raw = ['net_gen_tot_an_go', 'quantity_tot_an_po', 
                'elec_quantity_tot_an_po', 'tot_mmbtu_tot_an_po', 
                'elec_mmbtu_tot_an']
    ggp[vars_keep].to_parquet(snakemake.output.outfile)

    # save to csv
    ggp[vars_keep].to_csv(final_csv_dir / snakemake.output.outfile.split('/')[-1].replace('.parquet', '.csv'))


# %%
# READOUT FOR ALLSUMM
unique_cats = allsumm.columns.get_level_values(0).unique()
for cat in unique_cats:
    print(allsumm[allsumm.columns[allsumm.columns.get_level_values(0) == cat]])
# %%
# READOUT FOR ALIGNSUMM
temp = alignsumm.reset_index()
temp.loc[temp.level_1 != 'left_only'].pivot(index=['year', 'level_1'], columns='cat')
# %%
# summarize 2021 by category to spotcheck
# COMPARE GENRATION TO EIA REPORTS: https://www.eia.gov/tools/faqs/faq.php?id=427&t=3
# year = 2021
# catsumm = ggp.loc[ggp.year == year].groupby(['energy_source_1_cat', 'energy_source_1_subcat'])[['net_gen_tot_an']].sum() / 1e6
# catsumm['ratio'] = catsumm.net_gen_tot_an / catsumm.net_gen_tot_an.sum()
# catsumm.to_csv(results_dir / f'df_summ_{year}_totals.csv')
