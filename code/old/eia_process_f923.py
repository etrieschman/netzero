# %%
import pandas as pd
import numpy as np
from tqdm import tqdm

from utils import PATH_EIA, PATH_PROCESSED

YR_END = 2021
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 25)


# %%
# PROCESS DATA METHOD
def process_f923(yr_start):
    years = range(yr_start, YR_END+1)
    ysuff = {y:'Final_Revision' for y in years}
    ysuff[2022] = 'Early_Release'
    odf = pd.DataFrame({})
    for y in tqdm(years):
        # read in data
        df = pd.read_excel(f'{PATH_EIA}f923/{y}/EIA923_Schedules_2_3_4_5_M_12_{y}_{ysuff[y]}.xlsx', 
                        # header=5, sheet_name='Page 1 Generation and Fuel Data', na_values='.')
                        header=5, sheet_name='Page 4 Generator Data', na_values='.')
        df.columns = (df.columns.str.lower()
                    .str.replace(' ', '_')
                    .str.replace('\n', '_')
                    .str.replace('(', '').str.replace(')', ''))

        # transpose data
        # vars_id = ['year', 'plant_id', 'nuclear_unit_id', 'operator_id', 'plant_state',
        #     'eia_sector_number', 'reported_prime_mover',
        #     'reported_fuel_type_code', 'aer_fuel_type_code']
        vars_id = ['year', 'plant_id', 'operator_id', 'nerc_region', 'generator_id', 'reported_prime_mover']
        vars_val_sw = ('quantity', 'elec_quantity', 'tot_mmbtu', 'net_gen')
        vars_val = [col for col in df.columns if col.startswith(vars_val_sw)]
        vars_val = vars_val[:-1] # drop net generatio nyear to date
        dfl = df.melt(id_vars=vars_id, value_vars=vars_val)
        dfl['date'] = pd.to_datetime(dfl.year.astype(str) + 
                                    dfl['variable'].str.split('_').str[-1],
                                    format='%Y%B')
        dfl['variable'] = dfl['variable'].str.split('_').str[:-1].str.join('_')
        dfl = dfl.groupby(vars_id + ['variable'], dropna=False)['value'].sum().reset_index() 

        # concatenate
        odf = pd.concat([dfl, odf], axis=0, ignore_index=True)

        # output
        odf = odf.astype({'generator_id':str})
        odf.to_parquet(PATH_PROCESSED + 'eia_f923_ops.parquet', index=False)
    return odf

# %%
odf = process_f923(2018)


# %%
odf.reported_prime_mover.value_counts()


# %%
(odf
 .groupby(['year', 'variable'])['value'].agg(lambda x: (x>0).count())
 .reset_index()
 .pivot(index='variable', columns='year', values='value')
 .reset_index())
# %%
