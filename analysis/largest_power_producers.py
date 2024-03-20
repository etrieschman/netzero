# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm
import sys, os

# readin fuzzy name handling from cdp code
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
target_dir = os.path.join(parent_dir, 'code', 'align')
sys.path.append(target_dir)

from cdp_flag import preprocess_names, fuzzy_search_years

# global variables
PATH_DATA = '../data/'
PATH_FINAL = PATH_DATA + 'final/'
PATH_INTERIM = PATH_DATA + 'interim/'
PATH_RESOURCES = PATH_DATA + 'resources/'
PATH_RESULTS = '../results/analysis/largest_power_producers/'
rpath = Path(PATH_RESULTS)
rpath.mkdir(parents=True, exist_ok=True)

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)

# READ IN DATA
gdf = pd.read_parquet(PATH_FINAL + 'df_generators.parquet')
gendf = pd.read_parquet(PATH_FINAL + 'df_generation.parquet')
emdf = pd.read_parquet(PATH_FINAL + 'df_emissions.parquet')
pdf = pd.read_parquet(PATH_FINAL + 'df_plants.parquet')
udf = pd.read_parquet(PATH_FINAL + 'df_utilities.parquet')
odf = pd.read_parquet(PATH_FINAL + 'df_owners.parquet')

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
# utilities
ggeupdf = pd.merge(left=ggepdf, how='left',
                   right=udf[['year', 'utility_id', 'utility_name']],
                   on=['year', 'utility_id'])
assert len(ggeupdf) == len(ggepdf)

# subset to sample
ggeupdf_samp = ggeupdf.loc[ggeupdf.in_sample]
ggeupdf_samp


# %%
# READIN top100 power companies per MJBradley
toph = pd.read_csv(PATH_RESOURCES + 'top100_pp_mjbradley.csv', encoding='iso-8859-1')
toph.columns = toph.columns.str.lower().str.replace(' ', '_')
toph['name_cln'] = preprocess_names(toph.name)
toph['year'] = 2018

# %%
# FUZZY SEARCH
# format utility names
ggeupdf_samp['utility_name_cln'] = preprocess_names(ggeupdf_samp.utility_name)
ggeupdf_samp_sub = ggeupdf_samp.loc[ggeupdf_samp.year == 2018]

# fuzzy merge
map = fuzzy_search_years(ggeupdf_samp_sub, toph,
                   'utility_name_cln', 'name_cln', num_matches=1)
map = map.loc[map.name_cln_score >= 87]

# %%
# CHECK MERGE
# check toph merge
toph_merge = pd.merge(left=toph, right=map, on='name_cln', indicator=True)
print(toph_merge.groupby('_merge', observed=False)['name_cln'].count())

# check utility_merge
df = pd.merge(left=ggeupdf_samp_sub, right=map, on=['utility_name_cln', 'year'],
              how='left')
df['is_top_hundred'] = df.name_cln.notna()
print(df.is_top_hundred.sum())


# %%
# QUICK SUMMARY STATS MATCHED ONLY ON UTILITY NAME
# summarize split
df.groupby('is_top_hundred').agg({
    'utility_name_cln':'nunique',
    'generator_id':'count',
    'co2_mass_tons_gen_923':'sum',
    'co2_mass_tons_gen':'sum',
    'net_gen_tot_an':'sum'
}).T.style.format('{:,.2f}')

# %%
df.loc[df.is_top_hundred].groupby('name_cln').agg({
    'utility_name_cln':'nunique',
    'generator_id':'nunique',
    'plant_code':'nunique',
    'co2_mass_tons_gen_923':'sum',
    'co2_mass_tons_gen':'sum',
    'net_gen_tot_an':'sum'
}).sort_values('net_gen_tot_an', ascending=False).head(10).style.format('{:,.0f}')

# %%
# RANDOM QUESTION: DOES JENKINS ET AL TRANSFORMATION IMPACT UNCERTAINTY?
# MAKE DATA
import numpy as np
from scipy.stats import t
bootstraps = 5000
n = 75
alpha = 0.1
tval = t.ppf(1 - alpha/2, n-1)
results = {var:{'mean':[], 'var':[], 'std':[], 'unc':[]} for var in ['x', 'y']}

for i in tqdm(range(bootstraps)):
    x = np.random.normal(loc=2.5, scale=0.5, size=n)
    results['x']['mean'] += [x.mean()]
    results['x']['var'] += [x.var() / n]
    results['x']['std'] += [np.sqrt(x.var() / n)]
    results['x']['unc'] += [tval * np.sqrt(x.var() / n) / x.mean()]

    y = np.exp(-2 + 2.5*np.log(x))
    results['y']['mean'] += [y.mean()]
    results['y']['var'] += [y.var() / n]
    results['y']['std'] += [np.sqrt(y.var() / n)]
    results['y']['unc'] += [tval * np.sqrt(y.var() / n) / y.mean()]

# %%
# PLOT
var = 'unc'
plt.hist(results['x'][var], color='C0', alpha=0.5, label='x', density=True)
plt.hist(results['y'][var], color='C1', alpha=0.5, label='y', density=True)
plt.title(var)
plt.legend()
# %%
