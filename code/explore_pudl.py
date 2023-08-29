# %%
import intake
import pandas as pd
from pudl_catalog.helpers import year_state_filter

# %%
PATH_DATA = '../data/'
REF = 'v2022.11.30'
PUDL_INTAKE_CACHE = f'gs://intake.catalyst.coop/{REF}'
PUDL_INTAKE_PATH = PATH_DATA + 'pudl/'



# %%
pudl_cat = intake.cat.pudl_cat
list(pudl_cat)
# %%
pudl_cat.pudl.fuel_ferc1.read()
# %%
