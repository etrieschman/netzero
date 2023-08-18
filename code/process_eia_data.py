# %%
import pandas as pd
import numpy as np

from process_epa_data import readin_data

PATH_DATA = '../data/'
PATH_EIA = PATH_DATA + 'eia/'
PATH_EPA = PATH_DATA + 'epa/'


# %%
# READIN DATA
# TODO: READING IN THE WRONG EIA DATA. NEED PLANT-LEVEL
year = 2021
udf = pd.read_excel(f'{PATH_EIA}{year}/Utility_Data_2021.xlsx', header=1, sheet_name='States')
xwalk = pd.read_csv(f'{PATH_DATA}epa_eia_crosswalk.csv')
fdf = readin_data('facility')

# %%
# MERGE
fdf_ids = fdf[['Facility Name', 'Facility ID', 'Unit ID']].drop_duplicates()
xwalk_ids = xwalk[['CAMD_FACILITY_NAME', 'CAMD_PLANT_ID', 'CAMD_UNIT_ID',
                      'EIA_PLANT_NAME', 'EIA_PLANT_ID', 'EIA_GENERATOR_ID']].drop_duplicates()
udf_ids = udf[['Utility Name', 'Utility Number']].drop_duplicates().astype({'Utility Number':'float64'})
print('EPA units:', len(fdf_ids))
print('Crosswalk units:', len(xwalk_ids))
m = pd.merge(left=fdf_ids, right=xwalk_ids, how='outer', 
             left_on=['Facility ID', 'Unit ID'], 
            right_on=['CAMD_PLANT_ID', 'CAMD_UNIT_ID']).drop_duplicates()
print('Post-merge:', len(m))

print('Post-merge:', len(mm))

# %%
# SUMMARIZE MERGE
eiax_nepa = len(m.loc[m['Facility ID'].isna() & m['EIA_PLANT_ID'].notna()])
epa_nxeia = len(m.loc[m['Facility ID'].notna() & m['EIA_PLANT_ID'].isna()])
epa_dups = sum(m[['Facility Name', 'Facility ID', 'Unit ID']].value_counts() > 1)
print('EIA xwalk, missing in EPA:', eiax_nepa, eiax_nepa / len(m))
print('EPA, missing in EIA xwalk:', epa_nxeia, epa_nxeia / len(m))
print('EPA, duplicates:', epa_dups, epa_dups / len(m))

# %%
fdf
# %%
