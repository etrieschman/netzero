# %%
import numpy as np
import pandas as pd
import os, re

from utils import PATH_PROCESSED

# read in data
generator = pd.read_csv(PATH_PROCESSED + 'eia_f860_generator.csv')
plant = pd.read_csv(PATH_PROCESSED + 'eia_f860_plant.csv')
utility = pd.read_csv(PATH_PROCESSED + 'eia_f860_utility.csv')
owner = pd.read_csv(PATH_PROCESSED + 'eia_f860_ownership.csv')

# %%
print(len(generator), len(owner))
m = pd.merge(left=generator, right=owner, how='outer',
         on=['utility_id', 'plant_code', 'generator_id', 'year'])
# TODO: find owners not in generator file

#%%
