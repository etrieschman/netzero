# %%
import pandas as pd
import numpy as np

from process_epa_data import readin_data

PATH_DATA = '../data/'
PATH_EIA = PATH_DATA + 'eia/'
PATH_EPA = PATH_DATA + 'epa/'


# %%
# TEST MERGE
year = 2021
udf = pd.read_excel(f'{PATH_EIA}{year}/Utility_Data_2021.xlsx', header=1, sheet_name='States')
xwalk = pd.read_csv(f'{PATH_DATA}epa_eia_crosswalk.csv')
fdf = readin_data('facility')
# %%
fdf
# %%
