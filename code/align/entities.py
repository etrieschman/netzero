# %%
# SETUP
import pandas as pd
import numpy as np
from pathlib import Path

# global variables
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 150)

# SAMPLE SELECTION
def eia860_sample_selection(gdf, pdf, udf, odf):
    # 0. All
    colname = '0_all'
    pdf[colname] = True
    gdf[colname] = True
    udf[colname] = True
    odf[colname] = True

    # 1. drop non-contig states
    colname = '1_drop_noncontig_state'
    pdf[colname] = ~pdf.state_plant.isin(['AK', 'HI', 'PR'])
    gdf[colname] = gdf.plant_code.isin(pdf.loc[pdf[colname], 'plant_code'].values)
    udf[colname] = udf.utility_id.isin(pdf.loc[pdf[colname], 'utility_id'].values)
    odf[colname] = odf.plant_code.isin(pdf.loc[pdf[colname], 'plant_code'].values)

    # 2. drop CHP plants
    colname = '2_drop_chp'
    pdf[colname] = pdf.sector.isin([1,2,4,6]) | pdf.sector.isna()
    gdf[colname] = gdf.plant_code.isin(pdf.loc[pdf[colname], 'plant_code'].values)
    udf[colname] = udf.utility_id.isin(pdf.loc[pdf[colname], 'utility_id'].values)
    odf[colname] = odf.plant_code.isin(pdf.loc[pdf[colname], 'plant_code'].values)

    # SUMMARIZE
    steps = ['0_all', '1_drop_noncontig_state', '2_drop_chp']
    step_cols = []
    summ = pd.DataFrame()
    for step in steps:
        step_cols += [step]
        for idcol, df in zip(['plant_code', 'gid', 'utility_id', 'oid'], [pdf, gdf, udf, odf]):
            condition = np.array([df[step].values for step in step_cols]).all(0)
            summdf = (df.loc[condition].groupby('year')
                        .agg({idcol:'nunique'})
                        .rename(columns={idcol:(idcol, step)}))

            summ = pd.concat([summdf, summ], axis=1)

    summ = summ[summ.columns.sort_values(ascending=True)]
    summ.columns = pd.MultiIndex.from_tuples(summ.columns)

    # subset and return
    gdf['in_sample'] = np.array([gdf[step].values for step in steps]).all(0)
    pdf['in_sample'] = np.array([pdf[step].values for step in steps]).all(0)
    udf['in_sample'] = np.array([udf[step].values for step in steps]).all(0)
    odf['in_sample'] = np.array([odf[step].values for step in steps]).all(0)
    
    return (gdf, pdf, udf, odf), summ


# SIMPLIFY KEY COLUMNS: ENERGY SOURCE
# categories from EIA: https://www.eia.gov/tools/faqs/faq.php?id=427&t=3
subcat_to_energy_source = {
    'natural_gas':['NG'],
    'coal':['ANT', 'BIT', 'LIG', 'SGC', 'SUB', 'WC', 'RC'],
    'petroleum':['DFO', 'JF', 'KER', 'PC', 'PG', 'RFO', 'SGP', 'WO'],
    'other_gases':['BFG', 'OG'],
    'nuclear':['NUC'],
    'wind':['WND'],
    'hydropower':['WAT'],
    'solar':['SUN'],
    'geothermal':['GEO'],
    'waste_and_tires':['AB', 'MSW', 'OBS', 'WDS', 'WH', 'TDF'],
    'biomass':['OBL', 'SLW', 'BLQ', 'WDL', 'LFG', 'OBG'], 
    'batteries':['MWH'],
    'hydrogen':['H2'],
    'purchased_steam':['PUR'],
    'other':['OTH']
}
cat_to_subcat = {
    'fossil_fuels':['natural_gas', 'coal', 'petroleum', 'other_gases'],
    'nuclear':['nuclear'],
    'renewables':['wind', 'hydropower', 'solar', 'biomass', 'geothermal'],
    'other': ['waste_and_tires', 'biomass', 'batteries', 'hydrogen', 'purchased_steam', 'other']
}

def categorize_fuel(gdf, subcat_to_energy_source, cat_to_subcat):
    vars_source = gdf.columns[gdf.columns.str.contains('source') & ~gdf.columns.str.contains('raw')]
    source_to_subcat = {source:subcat for subcat, sources in subcat_to_energy_source.items() 
                        for source in sources}
    subcat_to_cat = {subcat:cat for cat, subcats in cat_to_subcat.items() 
                        for subcat in subcats}
    for var_source in vars_source:
        gdf.loc[:,f'{var_source}_subcat'] = gdf[var_source].map(source_to_subcat)
        gdf.loc[:,f'{var_source}_cat'] = gdf[f'{var_source}_subcat'].map(subcat_to_cat)
    return gdf

# %%
def fill_column_missing_values(df, ids, col):
    # Sort the DataFrame
    df[f'{col}_raw'] = df[col].copy()
    df.sort_values(by=ids, inplace=True)
    df[col] = df.groupby(ids)[col].ffill().bfill()
    return df


# %%
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('align_entities')
    
    results_dir = Path(snakemake.params.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # LOAD DATA
    gdf = pd.read_parquet(snakemake.input.infile_gen)
    pdf = pd.read_parquet(snakemake.input.infile_plant)
    udf = pd.read_parquet(snakemake.input.infile_util)
    odf = pd.read_parquet(snakemake.input.infile_own)
    # unique ids
    gdf['gid'] = gdf.plant_code.astype(str) + '_' + gdf.generator_id
    odf['oid'] = odf.ownership_id.astype(str) + '_' + odf.plant_code.astype(str) + '_' + odf.generator_id

    # SAMPLE SELECTION
    (gdf_sub, pdf_sub, udf_sub, odf_sub), summ = eia860_sample_selection(gdf, pdf, udf, odf)
    summ.to_csv(results_dir / 'sample_selection.csv')

    # CLEAN ATTRIBUTES
    # GENERATORS
    cols = ['nameplate_capacity_mw', 'dt_operation_start']
    cols += gdf_sub.columns[gdf.columns.str.contains('source')].to_list()
    ids = ['plant_code', 'generator_id', 'year']
    for col in cols:
        gdf_sub = fill_column_missing_values(gdf_sub, ids, col)
    gdf_sub['age'] = gdf_sub.year - gdf_sub.dt_operation_start.dt.year
    # PLANTS
    cols = ['nerc_region', 'latitude', 'longitude']
    ids = ['plant_code', 'year']
    for col in cols:
        pdf_sub = fill_column_missing_values(pdf_sub, ids, col)
    # fuel
    gdf_sub = categorize_fuel(gdf_sub, subcat_to_energy_source, cat_to_subcat)
    
    # WRITE TO FILE
    gdf_sub.to_parquet(snakemake.output.outfile_gen)
    pdf_sub.to_parquet(snakemake.output.outfile_plant)
    udf_sub.to_parquet(snakemake.output.outfile_util)
    odf_sub.to_parquet(snakemake.output.outfile_own)
# %%
