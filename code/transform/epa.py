# %%
import numpy as np
import pandas as pd
from pathlib import Path

from utils_transform import readin_epa
from utils_summ import summarize_id_counts_byyear

SHORT_TO_METRIC_TON = 0.907185

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
       'id': ['year', 'facility_id', 'unit_id', 'facility_name', 'state'],
       'val': ['sum_of_the_operating_time', 'gross_load_mwh', 'steam_load_1000_lb', 
                  'so2_mass_short_tons', 'co2_mass_short_tons', 'nox_mass_short_tons', 
                  'heat_input_mmbtu']}

if __name__ == '__main__':
       if "snakemake" not in globals():
              # readin mock snakemake
              import sys, os
              parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
              sys.path.insert(0, parent_dir)
              from utils import mock_snakemake
              snakemake = mock_snakemake('transform_epa')

       year_start = snakemake.params.year_start
       year_end = snakemake.params.year_end
       indir_facility = Path(snakemake.params.indir_facility)
       indir_emissions = Path(snakemake.params.indir_emissions)
       path_results = Path(snakemake.params.resultsdir)
       path_results.mkdir(parents=True, exist_ok=True)
       outfile_facility = Path(snakemake.output.outfile_facility)
       outfile_facility.parent.mkdir(parents=True, exist_ok=True)
       outfile_emissions = Path(snakemake.output.outfile_emissions)
       outfile_emissions.parent.mkdir(parents=True, exist_ok=True)

       print('Reading in facility data...')
       facdf = readin_epa(year_start, year_end, indir_facility, vars_keep=vars_keep_fac)
       facdf = facdf.astype({'unit_id':str})
       facdf.to_parquet(outfile_facility, index=False)

       print('Reading in daily emissions data...')
       emdf = readin_epa(year_start, year_end, indir_emissions, vars_coll=vars_coll_em)
       emdf = emdf.astype({'unit_id':str})
       emissions = ['so2', 'co2', 'nox']
       for emission in emissions:
              emdf[f'{emission}_mass_tons'] = emdf[f'{emission}_mass_short_tons'] * SHORT_TO_METRIC_TON
       emdf.to_parquet(outfile_emissions, index=False)

       print('Summarizing unique IDs over time...')
       print('Facility data:')
       facdf['uid'] = facdf.facility_id.astype(str) + '_' + facdf.unit_id.astype(str)
       out = summarize_id_counts_byyear(facdf.rename(columns={'facility_id':'fid'}), ['fid', 'uid'])
       out.to_csv(path_results / 'df_summ_final_facility.csv')
       print(out)
       print('Emissions data:')
       emdf['uid'] = emdf.facility_id.astype(str) + '_' + emdf.unit_id.astype(str)
       out = summarize_id_counts_byyear(emdf.rename(columns={'facility_id':'fid'}), ['fid', 'uid'])
       out.to_csv(path_results / 'df_summ_final_emissions.csv')
       print(out)

# %%
