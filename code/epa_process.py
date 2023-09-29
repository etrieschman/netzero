# %%
import numpy as np
import pandas as pd

from utils import PATH_EPA, PATH_PROCESSED, START_YEAR, END_YEAR

from utils_data_epa import readin_epa

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
facdf = readin_epa(START_YEAR, END_YEAR, PATH_EPA+'/facility/', vars_keep=vars_keep_fac)
facdf.to_csv(PATH_PROCESSED + 'epa_facility.csv', index=False)
facdf.head()

# %%
# PROCESS DAILY EMISSIONS DATA
vars_coll_em = {'id': ['facility_id', 'unit_id', 'facility_name', 'state'],
          'val': ['sum_of_the_operating_time', 'gross_load_mwh', 'steam_load_1000_lb', 
                  'so2_mass_short_tons', 'co2_mass_short_tons', 'nox_mass_short_tons', 
                  'heat_input_mmbtu']}
emdf = readin_epa(START_YEAR, END_YEAR, PATH_EPA+'/emissions/daily/', vars_coll=vars_coll_em)
emdf.to_csv(PATH_PROCESSED + 'epa_emissions.csv', index=False)
emdf.head()
