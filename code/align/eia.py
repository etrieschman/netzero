# %%
import pandas as pd

PATH_PROCESSED = '../data/processed/'

# READIN DATA
gdf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_generator.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_plant.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_utility.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'eia_f860_ownership.parquet')
opsdf = pd.read_parquet(PATH_PROCESSED + 'eia_f923_ops.parquet')
epafacdf = pd.read_parquet(PATH_PROCESSED + 'epa_facility.parquet')


# %%
# Align EIA entity properties (ones that persist over time)

# %%
# Make EIA mapping
# Need: EIA plant part to EIA generator mapping. Should be 1-to-many
# 1. split into generator and plant part generation
# 2. aportion operation data to generators and plant parts
#   -- get plant part net generation ex recorded generators
#   -- get % of net generation from recorded generators
# 3. Assign operation data to generators
#     -- Generation: 
#         -- directly where applicable
#         -- remaining generation, apportioned by % of capacity in plant part
#     -- Other operations:
#         -- percent of generation if in generator data
#         -- remaining operations, aportioned by % of capacity in plant part
# Question: How many plant parts don't exist in generator data?
# Question: How much operation gets assigned to discontinued and planning generators?
#         : Should we not assign values to these?

# %%
# Make EPA mapping
# Need: EPA unit to EPA plant part mapping. Is this many-to-1?
# If so, then we can aggregate EPA units to EPA plant parts,
# and apportion aggregated emissions from these units to 
# EIA generators in an EIA plant part


# %%
# Update generator operating dates based on 923 and EPA



