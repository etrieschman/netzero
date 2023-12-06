# %%
import pandas as pd
import numpy as np
import camelot
from pathlib import Path
# import os
# os.environ['GS_LIB'] = '/usr/local/Cellar/ghostscript/10.02.1/lib/libgs.dylib'
PATH_DRIVE = "C:/Users/etrieschman/My Drive/NetZero/duke energy/irps/"
PATH_RESULTS = '../results/irp/'
Path(PATH_RESULTS).mkdir(parents=True, exist_ok=True)
# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)

# %%
dd = {}
dd[2023] = {
    'load':{
        't':1,
        'i0':1,
        'rename':{0:'year', 7:'total_gwh'},
        'pages':['27']
        },
    'capacity':{
        't':0,
        'i0':3,
        'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
                  7:'reserve_margin_mw'},
        'pages':['73']
    },
    'changes':{
        't':0,
        'i0':7,
        'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
                  8:'dt_start', 9:'dt_retire', 11:'summer_firm_capacity_mw'},
        'pages':['75', '76']
    }}
dd[2022] = {
    'load':{
        't':1,
        'i0':1,
        'rename':{0:'year', 7:'total_gwh'},
        'pages':['25']
        },
    'capacity':{
        't':0,
        'i0':4,
        'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
                  7:'reserve_margin_mw'},
        'pages':['71']
        },
    'changes':{
        't':0,
        'i0':7,
        'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
                  7:'dt_start', 8:'dt_retire', 10:'summer_firm_capacity_mw'},
        'pages':['73']
    }}
dd[2021] = {
    'load':{
        't':1,
        'i0':1,
        'rename':{0:'year', 7:'total_gwh'},
        'pages':['23']
        },
    'capacity':{
        't':0,
        'i0':4,
        'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
                  7:'reserve_margin_mw'},
        'pages':['69']
        },
    'changes':{
        't':0,
        'i0':7,
        'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
                  7:'dt_start', 8:'dt_retire', 10:'summer_firm_capacity_mw'},
        'pages':['71']
    }}
dd[2020] = {
    'load':{
        't':1,
        'i0':1,
        'rename':{0:'year', 7:'total_gwh'},
        'pages':['20']
        },
    'capacity':{
        't':0,
        'i0':3,
        'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
                  7:'reserve_margin_mw'},
        'pages':['65']
        },
    'changes':{
        't':0,
        'i0':7,
        'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
                  7:'dt_start', 8:'dt_retire', 10:'summer_firm_capacity_mw'},
        'pages':['67']
    }}
dd[2019] = {
    'load':{
        't':1,
        'i0':1,
        'rename':{0:'year', 7:'total_gwh'},
        'pages':['22']
        },
    'capacity':{
        't':0,
        'i0':8,
        'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
                  7:'reserve_margin_mw'},
        'pages':['67']
        },
    'changes':{
        't':0,
        'i0':7,
        'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
                  10:'dt_start', 11:'dt_retire', 13:'summer_firm_capacity_mw'},
        'pages':['69']
    }}
# dd[2018] = {
#     'load':{
#         't':0,
#         'i0':8,
#         'rename':{0:'year', 7:'total_gwh'},
#         'pages':['21']
#         },
#     'capacity':{
#         't':0,
#         'i0':6,
#         'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
#                   7:'reserve_margin_mw'},
#         'pages':['66']
#         },
#     'changes':{
#         't':0,
#         'i0':6,
#         'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
#                   10:'dt_start', 11:'dt_retire', 12:'summer_firm_capacity_mw'},
#         'pages':['68']
#     }}
dd[2017] = {
    'load':{
        't':0,
        'i0':9,
        'rename':{0:'year', 7:'total_gwh'},
        'pages':['21']
        },
    'capacity':{
        't':0,
        'i0':8,
        'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
                  7:'reserve_margin_mw'},
        'pages':['66']
        },
    'changes':{
        't':0,
        'i0':6,
        'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
                  9:'dt_start', 10:'dt_retire', 12:'summer_firm_capacity_mw'},
        'pages':['68']
    }}
dd[2016] = {
    'load':{
        't':1,
        'i0':1,
        'rename':{0:'year', 7:'total_gwh'},
        'pages':['18']
        },
    'capacity':{
        't':0,
        'i0':8,
        'rename':{0:'year', 5:'firm_capacity_mw', 6:'summer_peak_load_mw',
                  7:'reserve_margin_mw'},
        'pages':['51']
        },
    'changes':{
        't':0,
        'i0':6,
        'rename':{0:'plant_name', 1:'unit_id', 3:'type', 4:'fuel_prim', 
                  9:'dt_start', 10:'summer_firm_capacity_mw'},
        'pages':['53']
    }}

dfs = {k:pd.DataFrame({}) for k in dd[2023].keys()}
for yr in dd.keys():
    for table in dd[yr].keys():
        print(yr, table)
        d = dd[yr][table].copy()
        for page in d['pages']:
            tables = camelot.read_pdf(PATH_DRIVE + f'Duke Energy Florida {yr}.pdf', 
                                    pages=page,
                                    split_text=False, flavor='stream')
            assert len(tables) > 0
            sdf = tables[d['t']].df[d['i0']:]
            sdf = sdf.rename(columns=d['rename'])[d['rename'].values()]
            sdf.reset_index(drop=True, inplace=True)
            sdf['irp_year'] = yr
            sdf['table'] = table
            dfs[table] = pd.concat([sdf, dfs[table]], axis=0, ignore_index=True)
# %%
# CLEAN LOAD
dfl = dfs['load'].copy()
dfl['year'] = np.where((dfl.year == '') | (dfl.year.str.startswith(('HIS', 'FORE'))), 
                       pd.NA, dfl.year)
dfl['year'] = dfl.year.bfill()
dfl['year'] = dfl.year.str.replace('ll', '11').astype(int)
dfl['total_gwh'] = pd.to_numeric(dfl.total_gwh.str.replace(',', ''), 
                                 errors='coerce')
dfl = dfl.loc[dfl.total_gwh.notna()].reset_index(drop=True)
dfl['year_rel'] = dfl.year - dfl.irp_year
summ_load = dfl.pivot(index='year_rel', columns='irp_year', values='total_gwh')
summ_load.to_csv(PATH_RESULTS + 'summ_load.csv', index=True)
# %%
# CLEAN CAPACITY
dfc = dfs['capacity'].copy()
for col in dfc.columns:
    if col in ['irp_year', 'table']:
        continue
    dfc[col] = pd.to_numeric(dfc[col].str.replace(',', ''), errors='coerce')
dfc['year_rel'] = dfc.year - dfc.irp_year
summ_cap = dfc.pivot(index='year_rel', columns='irp_year', values=['firm_capacity_mw', 'summer_peak_load_mw'])
summ_cap.to_csv(PATH_RESULTS + 'summ_capacity.csv', index=True)
# %%
# CLEAN CHANGES
dfr = dfs['changes'].copy()
dfr = dfr.loc[~dfr.plant_name.str.contains('DEGRADATION')]
dfr = dfr.loc[(dfr.summer_firm_capacity_mw != '') & 
              (dfr.unit_id != 'N/A')]
for col in dfr.columns:
    if col.startswith('dt'):
        dfr[col] = pd.to_datetime(dfr[col], errors='coerce')
    if col.endswith('_mw'):
        dfr[col] = pd.to_numeric(dfr[col].str.replace('(', '-').str.replace(')',''), errors='coerce')
mask_retire = (dfr.dt_start.dt.year >= dfr.irp_year) & (dfr.summer_firm_capacity_mw <= 0)
dfr['dt_retire'] = np.where(mask_retire, dfr.dt_start, dfr.dt_retire)
dfr['fuel_prim'] = dfr.fuel_prim.str.replace('NG\nDFO', 'NG')
dfr['fuel_prim'] = np.where((dfr.type == 'ST') & (dfr.fuel_prim==''), 'BIT', dfr.fuel_prim)
dfr['retirement'] = dfr.summer_firm_capacity_mw < 0
# summarize
summ_changes = dfr.groupby(['irp_year', 'retirement', 'type', 'fuel_prim'])['summer_firm_capacity_mw'].sum()
summ_changes.to_csv(PATH_RESULTS + 'summ_changes.csv', index=True)
# %%
