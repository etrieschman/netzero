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
# TODO: Align EIA entity properties (ones that persist over time)

# %%
# Remove plants not connected to the grid as SE does
# Data for plants that are not connected to the electrical grid 
# (These include all plants with a plant_id_eia of 88xxxx and the plants 
#  in this table: https://github.com/singularity-energy/open-grid-emissions/blob/main/data/manual/plants_not_connected_to_grid.csv)