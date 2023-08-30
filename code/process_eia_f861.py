# %%
import pandas as pd
import numpy as np
from tqdm import tqdm

from utils import PATH_EIA, PATH_PROCESSED
from utils import readin_eia
YR_START, YR_END = 2018, 2021

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 25)

# %%
# SUMMARIZE DATA AVAILABILITY
frame = readin_eia(YR_START, YR_END, f'{PATH_EIA}f861', 'frame_')
frame = frame.drop(columns='data_year')
vars_id = ['year', 'utility_number', 'utility_name', 'ownership_code',
           'ownership']
vars_form = [col for col in frame.columns if col not in vars_id]
frame['total'] = True
frame['none'] = frame[vars_form].any()
vars_form += ['total', 'none']
frame[vars_form] = frame[vars_form].replace(pd.NA, False).replace('X', True).replace('Y', True)

frame_summ = frame.groupby('year')[vars_form].sum()
frame_summ

# %%
# WRITE TO FILE
frame.to_csv(PATH_PROCESSED + 'eia_f861_uops.csv', index=False)









# %%
# DICTIONARY OF READ-IN HELPERS
rd = {
    'advanced_meters'       :{'sheets':['states', 'territories'], 'header':1},
    'demand_response'       :{'sheets':['demand response_states', 'demand response_territories'], 'header':2},
    'distribution_systems'  :{'sheets':None, 'header':0},
    'dynamic_pricing'       :{'sheets':None, 'header':0},
    'energy_efficiency'     :{'sheets':None, 'header':0},
    'short_form'            :{'sheets':None, 'header':0},
    'mergers'               :{'sheets':None, 'header':0}, 
    'net_metering'          :{'sheets':None, 'header':0},
    'non_net_metering_distributed': {'sheets':None, 'header':0}, 
    'operational_data'      :{'sheets':None, 'header':0},
    'reliablity'            :{'sheets':None, 'header':0},
    'sales_ult_cust'        :{'sheets':None, 'header':0},
    'sales_ult_cust_cs'     :{'sheets':None, 'header':0}, 
    'utility_data'          :{'sheets':None, 'header':0}
}

# EXAMPLE USAGE
years = range(2018, 2023)
ysuff = {y:'.xlsx' for y in years}
y = 2018
path_folder = f'{PATH_EIA}f861'
k = list(rd.keys())[0]
sdf = pd.DataFrame({})
if rd[k]['sheets'] is not None:
    for s in rd[k]['sheets']:
        df = pd.read_excel(f'{path_folder}/{y}/{k}_{y}{ysuff[y]}',
                           sheet_name=s, 
                           header=rd[k]['header'])
        sdf = pd.concat([df, sdf], axis=0, ignore_index=True)
else:
    sdf = pd.read_excel(f'{path_folder}/{y}/{k}_{y}{ysuff[y]}', 
                           header=rd[k]['header'])


# %%
