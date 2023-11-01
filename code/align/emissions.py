# %%
# SETUP
import pandas as pd

# global variables
PATH_PROCESSED = '../data/processed/'
DENOM, ROUND = 1e6, 2

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 150)


# %%
# Make EPA mapping
# Need: EPA unit to EPA plant part mapping. Is this many-to-1?
# If so, then we can aggregate EPA units to EPA plant parts,
# and apportion aggregated emissions from these units to 
# EIA generators in an EIA plant part