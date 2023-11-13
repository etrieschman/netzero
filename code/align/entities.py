# %%
# SETUP
import pandas as pd
import numpy as np

# global variables
PATH_PROCESSED = '../data/processed/'
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 150)

# %%
# LOAD DATASETS
def load_datasets(path=PATH_PROCESSED):
    gdf = pd.read_parquet(path + 'eia860_generator.parquet')
    pdf = pd.read_parquet(path + 'eia860_plant.parquet')
    udf = pd.read_parquet(path + 'eia860_utility.parquet')
    odf = pd.read_parquet(path + 'eia860_ownership.parquet')

    # unique ids
    gdf['gid'] = gdf.plant_code.astype(str) + '_' + gdf.generator_id
    odf['oid'] = odf.ownership_id.astype(str) + '_' + odf.plant_code.astype(str) + '_' + odf.generator_id

    return gdf, pdf, udf, odf
gdf, pdf, udf, odf = load_datasets()

# %%
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
    gdf_sub = gdf.loc[np.array([gdf[step].values for step in steps]).all(0)]
    pdf_sub = pdf.loc[np.array([pdf[step].values for step in steps]).all(0)]
    udf_sub = udf.loc[np.array([udf[step].values for step in steps]).all(0)]
    odf_sub = odf.loc[np.array([odf[step].values for step in steps]).all(0)]
    
    return (gdf_sub, pdf_sub, udf_sub, odf_sub), summ

(gdf_sub, pdf_sub, udf_sub, odf_sub), summ = eia860_sample_selection(gdf, pdf, udf, odf)
summ

# %%
# GENERATOR RETIREMENT DATE
# 1. Take retirement date from latest available year
gdf_sub['dt_operation_end_from_latest_year'] = np.where(gdf_sub.dt_operation_end.notna(), gdf_sub.year, np.nan)
gdf_sub['dt_operation_end_from_latest_year'] = gdf_sub.groupby('gid')['dt_operation_end_from_latest_year'].transform('max')
gdf_sub['dt_operation_end_from_latest_year'] = np.where(gdf_sub.dt_operation_end_from_latest_year == gdf_sub.year, gdf_sub.dt_operation_end, pd.NaT)
gdf_sub['dt_operation_end_from_latest_year'] = pd.to_datetime(gdf_sub.groupby('gid')['dt_operation_end_from_latest_year'].transform('max'))

# %%
# Does non-retirement status appear after retirement status?
gdf_sub.loc[(gdf_sub.status == 'RE') & 
            (gdf_sub.dt_operation_end_from_latest_year.dt.year > gdf_sub.year)]
# Yes! Look at case study:
# gdf_sub.loc[gdf_sub.gid == '1553_3']
# supported by Wiki: https://en.wikipedia.org/wiki/Gould_Street_Generating_Station


# %%
# TODO: Align EIA entity properties (ones that persist over time)