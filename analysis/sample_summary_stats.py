# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm

from utils import plot_mix

# global variables
PATH_DATA = '../data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
PATH_INTERIM = PATH_DATA + 'interim/'
PATH_RESOURCES = PATH_DATA + 'resources/'
PATH_RESULTS = '../results/analysis/sample_summary_stats/'
rpath = Path(PATH_RESULTS)
rpath.mkdir(parents=True, exist_ok=True)

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)

# READ IN DATA
gdf = pd.read_parquet(PATH_PROCESSED + 'df_generators.parquet')
gendf = pd.read_parquet(PATH_PROCESSED + 'df_generation.parquet')
emdf = pd.read_parquet(PATH_PROCESSED + 'df_emissions.parquet')
pdf = pd.read_parquet(PATH_PROCESSED + 'df_plants.parquet')
udf = pd.read_parquet(PATH_PROCESSED + 'df_utilities.parquet')
odf = pd.read_parquet(PATH_PROCESSED + 'df_owners.parquet')

# %%
# MAKE ANALYSIS DATASET
# generation
ggdf = pd.merge(left=gdf, right=gendf, how='outer',
                on=['year', 'utility_id', 'plant_code', 'generator_id'])
assert len(gdf) == len(ggdf)
# emissions
ggedf = pd.merge(left=ggdf, how='left',
                 right=emdf[['year', 'plant_code', 'generator_id', 
                             'co2_mass_tons_gen', 'co2_mass_tons_gen_923', 
                             'has_emissions_epa', 'has_emissions_eia']], 
                 on=['year', 'plant_code', 'generator_id'])
assert len(ggdf) == len(ggedf)
# plants
ggepdf = pd.merge(left=ggedf, how='left',
                 right=pdf[['year', 'utility_id', 'plant_code', 'nerc_region']],
                 on=['year', 'utility_id', 'plant_code'])
assert len(ggedf) == len(ggepdf)

# subset to sample
ggepdf_samp = ggepdf.loc[ggepdf.in_sample]
ggepdf_samp.loc[:,'all'] = 'all'


# %%
# SUMMARIZE EMISSIONS BY YEAR
emissions_summ_nerc = (ggepdf_samp
                  .groupby(['year', 'nerc_region'], dropna=False)
                  [['co2_mass_tons_gen', 'co2_mass_tons_gen_923']]
                  .sum()
                  .reset_index())
emissions_summ = (ggepdf_samp
                  .groupby(['year'], dropna=False)
                  [['co2_mass_tons_gen', 'co2_mass_tons_gen_923']]
                  .sum()
                  .reset_index())

# %%
# PLOT OVERALL
denom = 1e6
plt.plot(emissions_summ.year, emissions_summ.co2_mass_tons_gen/denom,
         linestyle='-', label='source: EPA')
plt.plot(emissions_summ.year, emissions_summ.co2_mass_tons_gen_923/denom,
         linestyle=':', label='source: EIA')
plt.legend()
plt.ylabel('CO2 emissions (M tons)')
plt.title('Emissions reported from EPA [solid] and EIA [dotted]')
plt.savefig(PATH_RESULTS + 'fig_emissions.png', 
            dpi=300, bbox_inches='tight')
plt.show()

# PLOT BY NERC REGION
for i, r in enumerate(emissions_summ_nerc.nerc_region.unique()):
    df = emissions_summ_nerc.loc[emissions_summ_nerc.nerc_region == r]
    plt.plot(df.year, df.co2_mass_tons_gen/denom, 
             linestyle='-', color=f'C{i}', label=r)
    plt.plot(df.year, df.co2_mass_tons_gen_923/denom, 
             linestyle=':', color=f'C{i}')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylabel('CO2 emissions (M tons)')
plt.title('Emissions by NERC region from EPA [solid] and EIA [dotted]')
plt.savefig(PATH_RESULTS + 'fig_emissions_by_nerc.png', 
            dpi=300, bbox_inches='tight')
plt.show()

# %%
# PLOT MIXES
# var, denom, units = 'net_gen_tot_an', 1e6, 'TWh'
var, denom, units = 'nameplate_capacity_mw', 1, 'MW'
col_loc='nerc_region'
show_pct = False
path = PATH_RESULTS + f'fig_{var}_pct{show_pct}.png'
plot_mix(ggepdf_samp, var, col_loc, path, 
         denom=denom, units=units, show_pct=show_pct, 
         ncols=3, scale=[6,5])
# %%
