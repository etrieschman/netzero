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

def readin_data(key, coll_cols=None):
    files = [f for f in os.listdir(PATH_EPA) if key in f.lower()]
    df = pd.DataFrame({})
    for f in tqdm(files):
        df_in = pd.read_csv(PATH_EPA + f, low_memory=False)
        if coll_cols is not None:
            df_in['quarter'] = pd.to_datetime(df_in.Date).dt.to_period('Q')
            df_in = df_in.groupby(coll_cols['id']+['quarter'])[coll_cols['val']].sum().reset_index()
        df = pd.concat([df_in, df], ignore_index=True, axis=0)
    return df

# %%
# READ IN data
# facdf = readin_data('facility')
# emdf = readin_data('emissions')

# %%
# READ IN EMISSIONS
col_em = {'id': ['State', 'Facility Name', 'Facility ID', 'Unit ID'],
          'val': ['Gross Load (MWh)', 'Steam Load (1000 lb)', 
                  'SO2 Mass (short tons)', 'CO2 Mass (short tons)', 
                  'NOx Mass (short tons)', 'Heat Input (mmBtu)']}
# emdf = readin_data('emissions', col_em)


# %%
# VISUALIZE EMISSIONS DATA
# facility_units = emdf[['Facility ID', 'Unit ID']].drop_duplicates().values
# for f, u in facility_units[:10]:
#     em = emdf.loc[(emdf['Facility ID'] == f) & (emdf['Unit ID'] == u), ].sort_values('quarter')
#     plt.plot(em.quarter.dt.to_timestamp().values, em['Gross Load (MWh)'].values)
# plt.title('Sample unit-level emissions data by quarter')

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
