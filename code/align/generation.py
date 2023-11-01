# %%
# SET UP
import pandas as pd

# global variables
PATH_PROCESSED = '../data/processed/'
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 150)

# %%
# HELPER FUNCTIONS
def get_annual_counts(df, groupvars='year', countvars='utility_id'):
    return df.groupby(groupvars)[countvars].count().values

def subset_from_ids(df, fromdf, idcol):
    return df.loc[df[idcol].isin(fromdf[idcol].drop_duplicates().values)]

def get_annual_generation(df, groupvars='year', sumvars='net_gen_tot_an', 
                          denom=DENOM, round=ROUND):
    return (df.groupby(groupvars)[sumvars].sum().values / denom).round(round)


# %%
# LOAD DATASETS
def load_datasets(path=PATH_PROCESSED):
    gdf = pd.read_parquet(path + 'eia860_generator.parquet')
    pdf = pd.read_parquet(path + 'eia860_plant.parquet')
    udf = pd.read_parquet(path + 'eia860_utility.parquet')
    odf = pd.read_parquet(path + 'eia860_ownership.parquet')
    opsdf = pd.read_parquet(path + 'eia923_ops.parquet')
    return gdf, pdf, udf, odf, opsdf


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
    print(get_annual_generation(opsdf, groupvars=['year', 'sheet'], sumvars=['net_gen_tot_mo', 'net_gen_tot_an']))

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
    opgendf = opgendf.groupby(['year', 'utility_id', 'plant_code', 'generator_id'], dropna=False)['net_gen_tot_an'].sum().reset_index()
    opppdf = opppdf.groupby(['year', 'utility_id', 'plant_code', 'reported_prime_mover', 'reported_fuel_type_code'], dropna=False)[opppdf.columns[opppdf.columns.str.endswith(('_an'))]].sum().reset_index()
    summ['n_id_generators_postcol'] = get_annual_counts(opgendf)
    summ['n_id_plantpart_postcol'] = get_annual_counts(opppdf)
    summ = pd.DataFrame(summ).T

    # return
    return opgendf, opppdf, summ



def clean_eia860(gdf, pdf, udf, odf):
    # 0. Sample selection on generators
    summ = {'gen':{}, 'pla':{}, 'uti':{}, 'own':{}}
    summ['gen']['0_n_total'] = get_annual_counts(gdf)
    summ['pla']['0_n_total'] = get_annual_counts(pdf)
    summ['uti']['0_n_total'] = get_annual_counts(udf)
    summ['own']['0_n_total'] = get_annual_counts(odf)
    # 1. drop non-contig states
    pdf_sub = pdf.loc[~pdf.state_plant.isin(['AK', 'HI', 'PR'])]
    gdf_sub = subset_from_ids(df=gdf, fromdf=pdf, idcol='plant_code')
    udf_sub = subset_from_ids(df=udf, fromdf=gdf, idcol='utility_id')
    odf_sub = subset_from_ids(df=odf, fromdf=gdf, idcol='generator_id')
    summ['gen']['1_n_drop_noncontig_state'] = get_annual_counts(gdf_sub)
    summ['pla']['1_n_drop_noncontig_state'] = get_annual_counts(pdf_sub)
    summ['uti']['1_n_drop_noncontig_state'] = get_annual_counts(udf_sub)
    summ['own']['1_n_drop_noncontig_state'] = get_annual_counts(odf_sub)
    # 2. drop CHP plants
    pdf_sub = pdf_sub.loc[pdf_sub.sector.isin([1,2,4,6])]
    gdf_sub = subset_from_ids(df=gdf_sub, fromdf=pdf_sub, idcol='plant_code')
    udf_sub = subset_from_ids(df=udf_sub, fromdf=gdf_sub, idcol='utility_id')
    odf_sub = subset_from_ids(df=odf_sub, fromdf=gdf_sub, idcol='generator_id')
    summ['gen']['2_n_drop_chp'] = get_annual_counts(gdf_sub)
    summ['pla']['2_n_drop_chp'] = get_annual_counts(pdf_sub)
    summ['uti']['2_n_drop_chp'] = get_annual_counts(udf_sub)
    summ['own']['2_n_drop_chp'] = get_annual_counts(odf_sub)
    # 3. drop generation units with no nameplate capacity
    gdf_sub = gdf_sub.loc[(gdf_sub.nameplate_capacity_mw > 0) & gdf_sub.nameplate_capacity_mw.notna()]
    pdf_sub = subset_from_ids(df=pdf_sub, fromdf=gdf_sub, idcol='plant_code')
    udf_sub = subset_from_ids(df=udf_sub, fromdf=gdf_sub, idcol='utility_id')
    odf_sub = subset_from_ids(df=odf_sub, fromdf=gdf_sub, idcol='generator_id')
    summ['gen']['3_n_drop_nocapacity'] = get_annual_counts(gdf_sub)
    summ['pla']['3_n_drop_nocapacity'] = get_annual_counts(pdf_sub)
    summ['uti']['3_n_drop_nocapacity'] = get_annual_counts(udf_sub)
    summ['own']['3_n_drop_nocapacity'] = get_annual_counts(odf_sub)
    # summary
    summ = pd.DataFrame(summ)
    return gdf_sub, pdf_sub, udf_sub, odf_sub, summ



# %% 
# ALIGN DATASETS
def align_eia923_to_eia860(gdf, opgendf, opppdf):
    # First, merge generator generation data onto EIA860 generator data
    summ = {'gen':{}, 'pp':{}}
    summ['gen']['n_ids'] = gdf.groupby('year')['utility_id'].count().values
    gg = pd.merge(
        right=opgendf, how='outer',
        left=gdf,
        on=['year', 'utility_id', 'plant_code', 'generator_id'],
        suffixes=('_g', '_go'), indicator=True)
    # summarize
    summ['gen']['n_ids_postmerge'] = get_annual_counts(gg.loc[gg._merge != 'right_only'])
    summ['gen']['n_ids_matched_generation'] = get_annual_counts(gg.loc[gg._merge == 'both'])
    summ['gen']['twh_matched_generation'] = get_annual_generation(gg.loc[gg._merge == 'both'])
    summ['gen']['n_unmatched_generation'] = get_annual_counts(gg.loc[gg._merge == 'right_only'])
    summ['gen']['twh_unmatched_generation'] = get_annual_generation(gg.loc[gg._merge == 'right_only'])
    gg = gg.loc[gg._merge != 'right_only'].rename(columns={'_merge':'merge_go'})

    # Second, merge plant part generation data onto EIA860 generator data
    ggp = pd.merge(
        left=gg, how='outer',
        right=opppdf,
        left_on=['year','utility_id', 'plant_code', 'prime_mover', 'energy_source_1'],
        right_on=['year', 'utility_id', 'plant_code', 'reported_prime_mover', 'reported_fuel_type_code'],
        suffixes=('', '_po'), indicator=True)
    # summarize
    summ['pp']['n_ids_postmerge'] = get_annual_counts(ggp.loc[ggp._merge != 'right_only'])
    summ['pp']['n_ids_matched_generation'] = get_annual_counts(ggp.loc[ggp._merge == 'both'])
    summ['pp']['twh_matched_generation'] = get_annual_generation(ggp.loc[ggp._merge == 'both'])
    summ['pp']['n_unmatched_generation'] = get_annual_counts(ggp.loc[ggp._merge == 'right_only'])
    summ['pp']['twh_unmatched_generation'] = get_annual_generation(ggp.loc[ggp._merge == 'right_only'])
    ggp = ggp.loc[ggp._merge != 'right_only'].rename(columns={'_merge':'merge_po'})

    # return
    summ = pd.DataFrame(summ)
    return ggp, summ



# ALLOCATE OPERATIONS TO GENERATORS
def allocate_eia923_operations(ggp):
    groupvars = ['year', 'plant_code', 'prime_mover', 'energy_source_1']
    summ = {}
    summ['n_generators'] = get_annual_counts(ggp)
    
    # flag source of generation data (gen or plant plart)
    ggp['flag_genfgo'] = (ggp.merge_go == 'both')
    summ['nid_w_generation_matchedto_generator'] = get_annual_counts(ggp.loc[ggp.flag_genfgo])
    summ['twh_from_generation_matchedto_generator'] = get_annual_generation(ggp.loc[ggp.flag_genfgo])
    ggp['net_gen_tot_an_go'] = ggp.groupby(groupvars)['net_gen_tot_an'].transform('sum')
    
    # get plant-part remaining generation
    ggp['net_gen_tot_an_po_remaining'] = (ggp.net_gen_tot_an_po - ggp.net_gen_tot_an_go).clip(0)
    summ['nid_w_generation_from_remaining_plantpart_gen'] = get_annual_counts(ggp.loc[(ggp.net_gen_tot_an_po_remaining > 0)])
    summ['twh_from_generation_from_remaining_plantpart_gen'] = get_annual_generation(ggp[groupvars + ['net_gen_tot_an_po_remaining']].drop_duplicates(), sumvars='net_gen_tot_an_po_remaining')
    
    # allocate net generation by capacity share
    print('Status of generators with reported generation:')
    print(ggp.loc[ggp.net_gen_tot_an > 0].groupby('status')['net_gen_tot_an'].agg(['sum', 'count']))
    # DECISION: allocate only to OP generators if available. The only to OA, OS, and SB if available. Then to All others if not
    ggp['group_has_opstatus'] = ggp.status == 'OP'
    ggp['group_has_opstatus'] = ggp.groupby(groupvars)['group_has_opstatus'].transform('max')
    ggp['group_has_otheractivestatus'] = ggp.status.isin(['OA', 'OS', 'SB'])
    ggp['group_has_otheractivestatus'] = ggp.groupby(groupvars)['group_has_otheractivestatus'].transform('max')
    # allocate to group with operation status
    ggp.loc[ggp.flag_genfgo, 'capacity_denominator'] = ggp.loc[ggp.flag_genfgo, 'nameplate_capacity_mw']
    mask = (~ggp.flag_genfgo) & (ggp.status == 'OP')
    ggp.loc[mask, 'capacity_denominator'] = ggp.loc[mask].groupby(groupvars)['nameplate_capacity_mw'].transform('sum')
    summ['nid_w_remaining_gen_added_to_opstatus'] = get_annual_counts(ggp.loc[mask])
    summ['twh_from_remaining_gen_added_to_opstatus'] = get_annual_generation(ggp.loc[mask, groupvars + ['net_gen_tot_an_po_remaining']].drop_duplicates(), sumvars='net_gen_tot_an_po_remaining')
    # allocate to group with other active status
    mask = (~ggp.flag_genfgo) & (ggp.status.isin(['OA', 'OS', 'SB'])) & (~ggp.group_has_opstatus)
    ggp.loc[mask, 'capacity_denominator'] = ggp.loc[mask].groupby(groupvars)['nameplate_capacity_mw'].transform('sum')
    summ['nid_w_remaining_gen_added_to_otheractivestatus'] = get_annual_counts(ggp.loc[mask])
    summ['twh_from_remaining_gen_added_to_otheractivestatus'] = get_annual_generation(ggp.loc[mask, groupvars + ['net_gen_tot_an_po_remaining']].drop_duplicates(), sumvars='net_gen_tot_an_po_remaining')
    # allocate to groups with no operation or active status
    mask = (~ggp.flag_genfgo) & (~ggp.group_has_opstatus) & (~ggp.group_has_otheractivestatus)
    ggp.loc[mask, 'capacity_denominator'] = ggp.loc[mask].groupby(groupvars)['nameplate_capacity_mw'].transform('sum')
    summ['nid_w_remaining_gen_added_to_remainingstatus'] = get_annual_counts(ggp.loc[mask])
    summ['twh_from_remaining_gen_added_to_remainingstatus'] = get_annual_generation(ggp.loc[mask, groupvars + ['net_gen_tot_an_po_remaining']].drop_duplicates(), sumvars='net_gen_tot_an_po_remaining')
    # final calculation
    ggp['pct_generation'] = ggp.nameplate_capacity_mw / ggp.capacity_denominator
    ggp.loc[~ggp.flag_genfgo, 'net_gen_tot_an'] = ggp.pct_generation.astype(float) * ggp.net_gen_tot_an_po_remaining
    
    # allocate other operational data
    ggp['pct_allocation'] = ggp.net_gen_tot_an / ggp.groupby(groupvars)['net_gen_tot_an'].transform('sum')
    num_vars = ['quantity_tot_an', 'elec_quantity_tot_an', 'tot_mmbtu_tot_an', 'elec_mmbtu_tot_an']
    for var in num_vars:
        ggp[f'{var}_po'] = ggp[var].copy()
        ggp[var] *= ggp['pct_allocation']

    # summary
    print('Allocation by status, 2021 (N, TWh):\n', ggp.loc[ggp.year == 2021].groupby('status').agg({'utility_id':'count', 'net_gen_tot_an':lambda x:(x.sum()/DENOM).round(ROUND)}))
    print('Allocation by energy source, 2021 (N, TWh):\n', ggp.loc[ggp.year == 2021].groupby('energy_source_1').agg({'utility_id':'count', 'net_gen_tot_an':lambda x:(x.sum()/DENOM).round(ROUND)}))

    summ = pd.DataFrame(summ).T
    return ggp, summ



# %%
if __name__ == '__main__':
    # load data
    gdf, pdf, udf, odf, opsdf = load_datasets()
    # clean datasets
    opsgendf, opsppdf, opssumm = clean_eia923(opsdf)
    print('Generation data summary:\n', opssumm)
    gdf, pdf, udf, odf, gensumm = clean_eia860(gdf, pdf, udf, odf)
    print('Generator data summary:\n', gensumm)
    # align datasets
    ggp, alignsumm = align_eia923_to_eia860(gdf, opsgendf, opsppdf)
    print('Aligning generation summary:\n', alignsumm)
    # allocate operations data
    ggp, allsumm = allocate_eia923_operations(ggp)
    print('Allocation summary:\n', allsumm)
    # write to file
    ggp.to_parquet(PATH_PROCESSED + 'df_generation.parquet')

# %%
