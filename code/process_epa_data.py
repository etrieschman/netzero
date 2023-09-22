# %%
import numpy as np
import pandas as pd
import geopandas as gpd
import plotly.figure_factory as ff
from shapely.geometry import Point
import matplotlib.pyplot as plt
from tqdm import tqdm

from utils import PATH_EPA, PATH_PROCESSED
from utils import readin_epa

YR_START, YR_END = 2013, 2021

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# %%
# PROCESS FACILITY DATA
vars_keep_fac = ['facility_id', 'unit_id', 'facility_name', 'state',
       'year', 'county', 'county_code', 'fips_code', 'latitude', 'longitude',
       'source_category',  'owner_operator', 'unit_type',
       'primary_fuel_type', 'secondary_fuel_type',
       'commercial_operation_date', 'operating_status',
       'max_hourly_hi_rate_mmbtu_hr', 'associated_generators_and_nameplate_capacity_mwe']
facdf = readin_epa(YR_START, YR_END, PATH_EPA+'/facility/', vars_keep=vars_keep_fac)
facdf.to_csv(PATH_PROCESSED + 'epa_facility.csv', index=False)
facdf.head()

# %%
# PROCESS DAILY EMISSIONS DATA
vars_coll_em = {'id': ['facility_id', 'unit_id', 'facility_name', 'state'],
          'val': ['sum_of_the_operating_time', 'gross_load_mwh', 'steam_load_1000_lb', 
                  'so2_mass_short_tons', 'co2_mass_short_tons', 'nox_mass_short_tons', 
                  'heat_input_mmbtu']}
emdf = readin_epa(YR_START, YR_END, PATH_EPA+'/emissions/daily/', vars_coll=vars_coll_em)
emdf.to_csv(PATH_PROCESSED + 'epa_emissions.csv', index=False)
emdf.head()


# %%
# VISUALIZE EMISSIONS DATA
facility_units = emdf[['facility_id', 'unit_id']].drop_duplicates().values
for f, u in facility_units[:10]:
    em = (emdf.loc[(emdf['facility_id'] == f) & (emdf['unit_id'] == u), ]
          .sort_values('quarter'))
    plt.plot(em.quarter.dt.to_timestamp().values, em['gross_load_mwh'].values)
plt.title('Sample unit-level emissions data by quarter')

# %%
# VISUALIZE FACILITY DATA
temp = (facdf
        .loc[~facdf.state.isin(['AK', 'HI', 'PR'])]
        .groupby(['latitude', 'longitude', 'year']).agg({'facility_id':'count'})
        .reset_index())
geometry = [Point(xy) for xy in zip(temp['longitude'], temp['latitude'])]
geotemp = gpd.GeoDataFrame(temp, geometry=geometry)
for y in geotemp.year.drop_duplicates().values:
    geotempyear = geotemp.loc[geotemp.year == y]
    geotempyear.plot(column='facility_id', alpha=0.5, legend=True)
    plt.title(f'Count of facilities in Year {y}')
    plt.show()

# %%
# GET COUNTS OF FACILITIES AND UNITS OVER TIME
facs = facdf
facs['funit_id'] = facs['facility_id'].astype(str) + '.' + facs['unit_id'].astype(str)
faclist = (facs.groupby('year')[['facility_id', 'unit_id']]
           .agg(lambda x: set(x.drop_duplicates().values)).reset_index())
faclist.columns = ['year', 'fac_ids', 'funit_ids']

prev_fids, prev_uids = set(), set()
counts = {}
counts['n_f'], counts['n_fnew'], counts['n_fdrop'] = [], [], []
counts['n_u'], counts['n_unew'], counts['n_udrop'] = [], [], []
for index, row in faclist.iterrows():
    curr_fids, curr_uids = row.fac_ids, row.funit_ids
    # Calculate overlap and new IDs
    counts['n_f'] += [len(curr_fids)]
    counts['n_u'] += [len(curr_uids)]
    counts['n_fnew'] += [len(curr_fids.difference(prev_fids))]
    counts['n_unew'] += [len(curr_uids.difference(prev_uids))]
    counts['n_fdrop'] += [len(prev_fids.difference(curr_fids))]
    counts['n_udrop'] += [len(prev_uids.difference(curr_uids))]
    # Update prev_ids for the next iteration
    prev_fids, prev_uids = curr_fids, curr_uids

# Add new columns to the DataFrame
for k, v in counts.items():
    faclist[k] = v
faclist = faclist.drop(columns=['fac_ids', 'funit_ids'])
faclist
# %%
