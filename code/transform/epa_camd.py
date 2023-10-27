# %%
import numpy as np
import pandas as pd

from utils_transform import PATH_RAW, PATH_PROCESSED, START_YEAR, END_YEAR
from utils_transform import readin_epa
from utils_summ import summarize_id_counts_byyear

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 25)

# %%
# PROCESS FACILITY DATA
vars_keep_fac = ['facility_id', 'unit_id', 'facility_name', 'state',
       'year', 'county', 'county_code', 'fips_code', 'latitude', 'longitude',
       'source_category',  'owner_operator', 'unit_type',
       'primary_fuel_type', 'secondary_fuel_type',
       'commercial_operation_date', 'operating_status',
       'max_hourly_hi_rate_mmbtu_hr', 'associated_generators_and_nameplate_capacity_mwe']
vars_coll_em = {
       'id': ['facility_id', 'unit_id', 'facility_name', 'state'],
       'val': ['sum_of_the_operating_time', 'gross_load_mwh', 'steam_load_1000_lb', 
                  'so2_mass_short_tons', 'co2_mass_short_tons', 'nox_mass_short_tons', 
                  'heat_input_mmbtu']}

if __name__ == '__main__':
       print('Reading in facility data...')
       facdf = readin_epa(START_YEAR, END_YEAR, PATH_RAW + 'epa/facility/', vars_keep=vars_keep_fac)
       facdf = facdf.astype({'unit_id':str})
       facdf.to_parquet(PATH_PROCESSED + 'epa_facility.parquet', index=False)

       print('Reading in daily emissions data...')
       emdf = readin_epa(START_YEAR, END_YEAR, PATH_RAW + 'epa/emissions/daily/', vars_coll=vars_coll_em)
       emdf = emdf.astype({'unit_id':str})
       emdf.to_parquet(PATH_PROCESSED + 'epa_emissions.parquet', index=False)

       print('Summarizing unique IDs over time...')
       print('Facility data:')
       facdf['uid'] = facdf.facility_id.astype(str) + '_' + facdf.unit_id.astype(str)
       print(summarize_id_counts_byyear(facdf.rename(columns={'facility_id':'fid'}), ['fid', 'uid']))
       print('Emissions data:')
       emdf['uid'] = emdf.facility_id.astype(str) + '_' + emdf.unit_id.astype(str)
       print(summarize_id_counts_byyear(emdf.rename(columns={'facility_id':'fid'}), ['fid', 'uid']))

# %%
