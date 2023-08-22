# %%
import numpy as np
import pandas as pd
import geopandas as gpd
import plotly.figure_factory as ff
from shapely.geometry import Point
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

PATH_DATA = '../data/'
PATH_EPA = PATH_DATA + 'epa/'

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def readin_epa_data(key, coll_cols=None):
    files = [f for f in os.listdir(PATH_EPA) if key in f.lower()]
    df = pd.DataFrame({})
    for f in tqdm(files):
        df_in = pd.read_csv(PATH_EPA + f, low_memory=False)
        if coll_cols is not None:
            df_in['quarter'] = pd.to_datetime(df_in.Date).dt.to_period('Q')
            df_in = (df_in.groupby(coll_cols['id']+['quarter'])[coll_cols['val']]
                     .sum().reset_index())
        df = pd.concat([df_in, df], ignore_index=True, axis=0)
    return df

# %%
# READ IN FACILITY DATA
facdf = readin_epa_data('facility')

# %%
# READ IN EMISSIONS DATA
col_em = {'id': ['State', 'Facility Name', 'Facility ID', 'Unit ID'],
          'val': ['Gross Load (MWh)', 'Steam Load (1000 lb)', 
                  'SO2 Mass (short tons)', 'CO2 Mass (short tons)', 
                  'NOx Mass (short tons)', 'Heat Input (mmBtu)']}
emdf = readin_epa_data('emissions', col_em)


# %%
# VISUALIZE EMISSIONS DATA
facility_units = emdf[['Facility ID', 'Unit ID']].drop_duplicates().values
for f, u in facility_units[:10]:
    em = (emdf.loc[(emdf['Facility ID'] == f) & (emdf['Unit ID'] == u), ]
          .sort_values('quarter'))
    plt.plot(em.quarter.dt.to_timestamp().values, em['Gross Load (MWh)'].values)
plt.title('Sample unit-level emissions data by quarter')

# %%
# VISUALIZE FACILITY DATA
# temp = (facdf
#         .loc[~facdf.State.isin(['AK', 'HI', 'PR'])]
#         .groupby(['Latitude', 'Longitude', 'Year']).agg({'Facility ID':'count'})
#         .reset_index())
# geometry = [Point(xy) for xy in zip(temp['Longitude'], temp['Latitude'])]
# geotemp = gpd.GeoDataFrame(temp, geometry=geometry)
# for y in geotemp.Year.drop_duplicates().values:
#     geotempyear = geotemp.loc[geotemp.Year == y]
#     geotempyear.plot(column='Facility ID', alpha=0.5, legend=True)
#     plt.title(f'Count of facilities in Year {y}')
#     plt.show()

# %%
# SUBSET DATA
keep_sourcecat = ['Electric Utility', 'Small Power Producer','Cogeneration']
facs = facdf.loc[facdf['Source Category'].isin(keep_sourcecat)]
facs['funitID'] = facs['Facility ID'].astype(str) + '.' + facs['Unit ID'].astype(str)

# %%
# GET COUNTS OF FACILITIES AND UNITS OVER TIME
faclist = (facs.groupby('Year')[['Facility ID', 'funitID']]
           .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
faclist.columns = ['Year', 'FacIDs', 'FunitIDs']

prev_fids, prev_uids = set(), set()
counts = {}
counts['n_foverlap'], counts['n_fnew'], counts['n_fdrop'] = [], [], []
counts['n_uoverlap'], counts['n_unew'], counts['n_udrop'] = [], [], []
for index, row in faclist.iterrows():
    curr_fids, curr_uids = row.FacIDs, row.FunitIDs
    # Calculate overlap and new IDs
    counts['n_foverlap'] += [len(prev_fids.intersection(curr_fids))]
    counts['n_uoverlap'] += [len(prev_uids.intersection(curr_uids))]
    counts['n_fnew'] += [len(curr_fids.difference(prev_fids))]
    counts['n_unew'] += [len(curr_uids.difference(prev_uids))]
    counts['n_fdrop'] += [len(prev_fids.difference(curr_fids))]
    counts['n_udrop'] += [len(prev_uids.difference(curr_uids))]
    # Update prev_ids for the next iteration
    prev_fids, prev_uids = curr_fids, curr_uids

# Add new columns to the DataFrame
for k, v in counts.items():
    faclist[k] = v
faclist = faclist.drop(columns=['FacIDs', 'FunitIDs'])
faclist['n_f'] = faclist.n_foverlap + faclist.n_fnew + faclist.n_fdrop
faclist['n_u'] = faclist.n_uoverlap + faclist.n_unew + faclist.n_udrop
faclist


# %%
facs
# %%
