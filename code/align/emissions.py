# %%
# SET UP
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# global variables
DENOM, ROUND = 1e6, 2
SHORT_TO_METRIC_TON = 0.907185

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 250)

# %% 
# CLEAN UP CROSSWALK
def clean_xw(xw):
    # drop and rename columns
    xw.columns = xw.columns.str.lower()
    xw.drop(columns=[col for col in xw.columns if col.startswith('mod_')], inplace=True)
    xw.drop(columns=['sequence_number', 'camd_state', 'camd_facility_name',
                    'camd_latitude', 'camd_longitude', 'camd_nameplate_capacity',
                    'eia_state', 'eia_plant_name', 'eia_latitude', 'eia_longitude',
                    'eia_nameplate_capacity', 'plant_id_change_flag',
                    'match_type_boiler'], inplace=True)
    # check: make sure we're not dropping a unique identifying column
    assert len(xw) == len(xw.drop_duplicates())

    # drop CAMD excluded and unmatched
    print('Dropping records labeled "excluded or unmatched"...')
    print('Rows dropped:\t', len(xw.loc[xw.match_type_gen.isin(['Manual CAMD Excluded', 'CAMD Unmatched'])]))
    xw = xw.loc[~xw.match_type_gen.isin(['Manual CAMD Excluded', 'CAMD Unmatched'])]
    # drop boiler id and collapse on remaining unique ids
    print('Dropping EIA boiler IDs and recollapsing...')
    print('Rows dropped:\t', xw.drop(columns='eia_boiler_id').duplicated().sum())
    xw = xw.drop(columns='eia_boiler_id').drop_duplicates()
    # drop nameplate capacity information and recollapse
    print('Dropping CAMD generator IDs and recollapsing...')
    print('Rows dropped:\t', xw.drop(columns=['camd_generator_id']).duplicated().sum())
    xw = xw.drop(columns='camd_generator_id').drop_duplicates()
    # update datatypes
    xw['eia_plant_id'] = pd.to_numeric(xw.eia_plant_id).astype(int)
    # make ids
    xw['camd_uid'] = 'camd_' + xw.camd_plant_id.astype(str) + '_' + xw.camd_unit_id
    xw['eia_gid'] = 'eia_' + xw.eia_plant_id.astype(str) + '_' + xw.eia_generator_id
    # check: rows unique on the identifiers we want
    assert len(xw) == len(xw[['camd_uid', 'eia_gid']].drop_duplicates())
    print('Returning clean and deduped crosswalk...')
    return xw


# %%
# Create subplant associations
def create_subplants(xw):
    # make edgegraph and assert that it's bipartite
    graph = nx.from_pandas_edgelist(
        xw, source='camd_uid', target='eia_gid', edge_attr=True)
    assert nx.algorithms.bipartite.is_bipartite(graph)

    # label subgraphs
    for i, node_set in enumerate(nx.connected_components(graph)):
        subgraph = graph.subgraph(node_set)
        nx.set_edge_attributes(subgraph, name='subplant_id', values=f'sp_{i}')

    xw_sp = nx.to_pandas_edgelist(graph).rename({'source':'camd_uid', 'target':'eia_gid'})
    return xw_sp

# %%
# ALIGN DATASETS
def align_datasets(edf, gdf, xwc, vars_em):
    # 0. make ID tables
    xwid = xwc[['camd_uid', 'eia_gid']].copy()
    edfid = edf[['year', 'facility_id', 'unit_id'] + vars_em].copy()
    edfid['uid'] = 'camd_' + edfid.facility_id.astype(str) + '_' + edfid.unit_id
    gdfid = gdf[['year', 'plant_code', 'generator_id', 'nameplate_capacity_mw', 'status']].copy()
    gdfid['gid'] = 'eia_' + gdfid.plant_code.astype(str) + '_' + gdfid.generator_id
    # CHECK: id tables preserve unique identifiers
    assert edfid.duplicated().sum() == 0
    assert gdfid.duplicated().sum() == 0

    # 1. Merge together
    # 1.a Merge emissions id tables with crosswalk
    edfid_xw = (pd.merge(left=edfid, right=xwid, how='outer',
                        left_on='uid', right_on='camd_uid', indicator='merge_e')
                        .drop(columns='camd_uid'))
    summ_medf = edfid_xw.groupby(['year', 'merge_e'], observed=False).agg({'uid':'nunique', 'co2_mass_tons':'sum'})
    summ_medf.columns = [f'{col}_xw' for col in summ_medf.columns]
    edfid_xw = edfid_xw.loc[edfid_xw.merge_e != 'right_only']
    # 1.b Merge generator id tables with crosswalk
    gdfid_xw = (pd.merge(left=gdfid, right=xwid, how='outer',
                        left_on=['gid'], right_on=['eia_gid'], indicator='merge_g')
                        .drop(columns='eia_gid'))
    summ_mgdf = gdfid_xw.groupby(['year', 'merge_g'], observed=False).agg({'gid':'nunique', 'nameplate_capacity_mw':'sum'})
    summ_mgdf.columns = [f'{col}_xw' for col in summ_mgdf.columns]
    gdfid_xw = gdfid_xw.loc[gdfid_xw.merge_g != 'right_only']
    # 1.c Merge emisisons with generator ids 
    gdfid_xw_edfid = pd.merge(left=gdfid_xw, right=edfid_xw, how='outer',
                            left_on=['year', 'gid', 'camd_uid'], right_on=['year', 'eia_gid', 'uid'], indicator='merge')
    summ_mgdfedf = (gdfid_xw_edfid.groupby(['year', 'merge'], observed=False)
                    .agg({'gid':'nunique', 'uid':'nunique', 
                        'co2_mass_tons':'sum',
                        'nameplate_capacity_mw':'sum'}))
    
    # 2. Summarize merge
    summ_m = pd.concat([summ_medf, summ_mgdf, summ_mgdfedf], axis=1).reset_index()
    summ_m['both'] = summ_m.level_1 == 'both'
    summ_m = summ_m.drop(columns='level_1').groupby(['year', 'both']).sum()
    summ_m_pct = summ_m / summ_m.groupby(['year']).transform('sum')
    summ_m_pct.columns = [f'{col}_pct' for col in summ_m_pct.columns]
    summ_m = pd.concat([summ_m, summ_m_pct], axis=1)
    summ_m = summ_m[summ_m.columns.sort_values()]

    
    return gdfid_xw_edfid, summ_m

# %%
# ALLOCATE EMISSIONS TO GENERATORS
def allocate_emissions_to_generators(g_x_e, vars_em):
    # 1. Calculate the percent of nameplate capacity per unit ID
    g_x_e['pct_allocation'] = g_x_e.nameplate_capacity_mw / g_x_e.groupby(['year', 'uid'])['nameplate_capacity_mw'].transform('sum')
    g_x_e['is_allocated'] = g_x_e.pct_allocation.notna() & g_x_e.co2_mass_tons.notna()

    # 2. Multiply emissions values by % nameplate capacity
    for var in vars_em:
        g_x_e[f'{var}_gen'] = g_x_e[var].copy()
        g_x_e[f'{var}_gen'] *= g_x_e['pct_allocation']

    # 3. Collapse on generator IDs
    agg = {f'{var}_gen':'sum' for var in vars_em}
    agg['is_allocated'] = 'max'
    gdf_e = (g_x_e
            .groupby(['year', 'plant_code', 'generator_id', 'gid', 
                      'nameplate_capacity_mw'], dropna=False)
            .agg(agg)
            .reset_index())
    gdf_e.loc[~gdf_e.is_allocated, [f'{var}_gen' for var in vars_em]] = pd.NA
    
    # 4. Check results
    summ_all = g_x_e[['year', 'camd_uid', 'is_allocated'] + vars_em].drop_duplicates().drop(columns='camd_uid').groupby(['year', 'is_allocated']).sum()/1e6
        
    return gdf_e, summ_all

# %%
# CALCULATE EMISSIONS FROM EIA 923
def calculate_eia_emissions(gdf, gendf, ef):
    efc = ef[['energy_source_code', 'co2_tons_per_mmbtu']].copy()
    gendfc = pd.merge(
        left=gdf[['year', 'plant_code', 'generator_id', 'energy_source_1']],
        right=gendf[['year', 'plant_code', 'generator_id',
                    'elec_mmbtu_tot_an', 'tot_mmbtu_tot_an']],
        how='left', on=['year', 'plant_code', 'generator_id'])
    gendfc_ef = pd.merge(
        left=gendfc, right=efc, how='left',
        left_on='energy_source_1', right_on='energy_source_code')
    gendfc_ef['co2_mass_short_tons_gen_923'] = gendfc_ef.co2_tons_per_mmbtu * gendfc_ef.tot_mmbtu_tot_an
    gendfc_ef['co2_mass_tons_gen_923'] = gendfc_ef.co2_mass_short_tons_gen_923 * SHORT_TO_METRIC_TON
    return gendfc_ef

# %%
def compare_emissions(df, var_em_epa, var_em_eia, path_results, hist_thresh=1e4):
    # table summary
    summ = (df
            .groupby(['year', 'has_emissions_epa', 'has_emissions_eia'], dropna=False)
            [['co2_mass_tons_gen', 'co2_mass_tons_gen_923']].sum()) / 1e6
    
    # scatterplot
    df[[var_em_epa, var_em_eia]] = df[[var_em_epa, var_em_eia]].fillna(0)
    p = sns.FacetGrid(df, col='year', col_wrap=4) 
    p.map(plt.scatter, var_em_eia, var_em_epa, alpha=0.5)
    p.map(plt.axline, xy1=(0,0), slope=1, color='C1')
    plt.savefig(path_results / 'fig_summ_eia_v_epa_scatter.png', dpi=300)

    # histograms
    fig, ax = plt.subplots(nrows=2)
    dif = df[var_em_eia] - df[var_em_epa]
    abs_dif_gt_tresh = np.abs(dif.loc[dif > hist_thresh])
    ax[0].hist(dif, bins=100, label=f'Difference')
    ax[1].hist(abs_dif_gt_tresh, bins=100, color='C1', label=f'Absolute difference > {hist_thresh}')
    ax[0].set_yscale('log')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[0].set_title('Histogram of difference between emissions from EIA923 and EPA')
    ax[1].legend()
    ax[0].legend()
    plt.savefig(path_results / 'fig_summ_eia_v_epa_hist.png', dpi=300)

    return summ


# %% RUN SCRIPT
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('align_emissions')
    
    results_dir = Path(snakemake.params.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    final_csv_dir = Path(snakemake.params.final_csv_dir)
    final_csv_dir.mkdir(parents=True, exist_ok=True)
    
    print('\nReading in data...')
    gdf = pd.read_parquet(snakemake.input.infile_gen)
    gendf = pd.read_parquet(snakemake.input.infile_ops)
    edf = pd.read_parquet(snakemake.input.infile_epa_emissions)
    xw = pd.read_csv(snakemake.input.crosswalk)
    ef = pd.read_csv(snakemake.input.emission_factors, header=1, na_values='*')

    print('\nCleaning crosswalk...')
    xwc = clean_xw(xw.copy())

    print('\nAligning EPA and EIA datasets...')
    vars_em = [
        'sum_of_the_operating_time', 'gross_load_mwh', 'steam_load_1000_lb',
        'so2_mass_tons', 'co2_mass_tons', 'nox_mass_tons',
        'heat_input_mmbtu']
    g_x_e, summ_all = align_datasets(edf, gdf, xwc, vars_em)
    summ_all.to_csv(results_dir / 'df_summ_align.csv')

    # CHECK: How are emissions distributed across statuses?
    print('Check: emissions across generator status:')
    print(g_x_e.loc[g_x_e.co2_mass_tons > 0].groupby(['status'], dropna=False)['co2_mass_tons'].sum())

    print('\nAllocating EPA emissions to EIA generators...')    
    gen_e, summ = allocate_emissions_to_generators(g_x_e, vars_em)
    summ.to_csv(results_dir / 'df_summ_allocate.csv')

    print('\nCalculating emisisons from 923 data...')
    gen_e923 = calculate_eia_emissions(gdf, gendf, ef)
    gen_ee = pd.merge(left=gen_e, how='outer',
                      right=gen_e923[['year', 'plant_code', 'generator_id', 'co2_tons_per_mmbtu', 'co2_mass_tons_gen_923']],
                        on=['year', 'plant_code', 'generator_id'])
    gen_ee['has_emissions_epa'] = gen_ee.co2_mass_tons_gen.notna()
    gen_ee['has_emissions_eia'] = gen_ee.co2_mass_tons_gen_923.notna()
    
    print('Check: how do emissions compare?')
    summ_compare = compare_emissions(gen_ee.copy(), var_em_epa='co2_mass_tons_gen', 
                                     var_em_eia='co2_mass_tons_gen_923',
                                     path_results=results_dir)
    summ_compare.to_csv(results_dir / 'df_summ_eia_v_epa.csv')
    
    print('\nWriting to file...')
    gen_ee.to_parquet(snakemake.output.outfile)
    gen_ee.to_csv(final_csv_dir / snakemake.output.outfile.split('/')[-1].replace('.parquet', '.csv'))
    
    print('\nDone!')


# %%
summ_all[summ_all.columns[~summ_all.columns.str.endswith(('_xw', '_xw_pct'))]]


# %%
