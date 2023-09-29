# %%
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

from utils_summ import summarize_id_counts_byyear
from utils import PATH_PROCESSED

# %%
# READ IN DATA
facdf = pd.read_csv(PATH_PROCESSED + 'epa_facility.csv')
emdf = pd.read_csv(PATH_PROCESSED + 'epa_emissions.csv')

# %% 
# SUMMARIZE UNIQUE IDENTIFIERS
# emdf = pd.read_csv(PATH_PROCESSED + 'epa_emissions.csv')
emdf['pid'] = emdf.facility_id
emdf['uid'] = emdf.facility_id.astype(str) + '_' + emdf.unit_id.astype(str)
summarize_id_counts_byyear(emdf, ['pid', 'uid'])

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

