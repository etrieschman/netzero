# %%
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

from utils import PATH_DATA, PATH_PROCESSED
from utils_summ import summarize_merge, summarize_plant_inepa, summarize_gen_inepa

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 50)

# %%
# IMPORT
eia_own = pd.read_parquet(PATH_PROCESSED + 'eia_final.parquet')
epa = pd.read_csv(PATH_PROCESSED + 'epa_emissions.csv')
xwalk = pd.read_csv(PATH_DATA + 'epa_eia_crosswalk.csv')
eia_own['key_gen'] = eia_own.plant_code.astype(str) + '_' + eia_own.generator_id
eia_own['key_own'] = eia_own.key_gen + '_' + eia_own.ownership_id.astype(str)
epa['key_unit'] = epa.facility_id.astype(str) + '_' + epa.unit_id.astype(str)
eia = eia_own[[col for col in eia_own.columns if 'own' not in col]].drop_duplicates()

# %%
# CLEAN CROSSWALK
xwalk.columns = xwalk.columns.str.lower()
nrows_xwalk_raw = len(xwalk)
xwalk = xwalk.loc[~xwalk.match_type_gen.isin(
    ['Manual CAMD Excluded', 'CAMD Unmatched'])]
print('CAMD units dropped, either excluded or unmatched:', nrows_xwalk_raw - len(xwalk))
xwalk = xwalk.astype({'eia_plant_id':'Int64', 'camd_plant_id':'Int64'})
# DECISION drop boiler ID and epa generator ID and dedup on remaining crosswalk
keepcols_xwalk = ['camd_plant_id', 'camd_unit_id',
       'eia_plant_id', 'eia_generator_id']
nrows_xwalk_boiler = len(xwalk)
xwalk = xwalk[keepcols_xwalk].drop_duplicates()
print('Crosswalk rows dropped for unique EIA boilers and EIA generators:', nrows_xwalk_boiler - len(xwalk))



# %%
# UNDERSTAND MERGE EIA ON CROSSWALK
eia_on = ['plant_code', 'generator_id']
xwalk_eia_on = ['eia_plant_id', 'eia_generator_id'] 
m_eia = pd.merge(left=eia, right=xwalk, how='outer', 
             left_on=eia_on, 
             right_on=xwalk_eia_on)
summarize_merge(m_eia, eia_on, xwalk_eia_on, r_hasyear=False)

# %%
# UNDERSTAND MERGE EPA ON CROSSWALK
epa_col = epa.groupby(['facility_id', 'unit_id', 'key_unit', 'year'])[['co2_mass_short_tons', 'gross_load_mwh']].sum().reset_index()
epa_on = ['facility_id', 'unit_id']
xwalk_epa_on = ['camd_plant_id', 'camd_unit_id']
m_epa = pd.merge(left=epa_col, right=xwalk, how='outer',
              left_on=epa_on, right_on=xwalk_epa_on)
m_epa['year'] = m_epa.year.astype('Int64')
display(summarize_merge(m_epa, epa_on, xwalk_epa_on, r_hasyear=False))

# %%
# UNDERSTAND OVERALL MERGE
# pull in EIA data
m_epa_eia = pd.merge(left=m_epa, right=eia, how='outer',
                    left_on=['eia_plant_id', 'eia_generator_id', 'year'],
                    right_on=['plant_code', 'generator_id', 'year'])
# drop crosswalk when not in either
m_epa_eia = m_epa_eia.loc[m_epa_eia.key_unit.notna() | m_epa_eia.key_gen.notna()]
display(summarize_merge(m_epa_eia, ['eia_plant_id', 'eia_generator_id'], ['plant_code', 'generator_id'], r_hasyear=True))
# summarize plant mappings
m_epa_eia['n_eia_plants_pepa'] = m_epa_eia.groupby(['facility_id', 'year'])['plant_code'].transform('nunique')
m_epa_eia['n_epa_plants_peia'] = m_epa_eia.groupby(['plant_code', 'year'])['facility_id'].transform('nunique')
display(m_epa_eia[['plant_code', 'facility_id', 'n_eia_plants_pepa']].drop_duplicates()['n_eia_plants_pepa'].value_counts().sort_index())
display(m_epa_eia[['plant_code', 'facility_id', 'n_epa_plants_peia']].drop_duplicates()['n_epa_plants_peia'].value_counts().sort_index())
# summarize generator/unit mappings
m_epa_eia['n_eia_gen_per_epa_unit'] = m_epa_eia.groupby(['key_unit', 'year'])['key_gen'].transform('nunique')
m_epa_eia['n_epa_unit_per_eia_gen'] = m_epa_eia.groupby(['key_gen', 'year'])['key_unit'].transform('nunique')
m_epa_eia.loc[(m_epa_eia.n_eia_gen_per_epa_unit == 1) & (m_epa_eia.n_epa_unit_per_eia_gen == 1), 'eia_epa_merge_type'] = '1_to_1'
m_epa_eia.loc[(m_epa_eia.n_eia_gen_per_epa_unit > 1) & (m_epa_eia.n_epa_unit_per_eia_gen == 1), 'eia_epa_merge_type'] = 'many_to_1'
m_epa_eia.loc[(m_epa_eia.n_eia_gen_per_epa_unit == 1) & (m_epa_eia.n_epa_unit_per_eia_gen > 1), 'eia_epa_merge_type'] = '1_to_many'
m_epa_eia.loc[(m_epa_eia.n_eia_gen_per_epa_unit > 1) & (m_epa_eia.n_epa_unit_per_eia_gen > 1), 'eia_epa_merge_type'] = 'many_to_many'
m_epa_eia.loc[(m_epa_eia.n_epa_unit_per_eia_gen == 0), 'eia_epa_merge_type'] = '1_to_na'
print('Count of unique generators, by merge type')
display(m_epa_eia.groupby(['eia_epa_merge_type', 'year'])['key_gen'].nunique().reset_index().pivot(index='year', columns='eia_epa_merge_type'))

# %%
# COLLAPSE BACK TO EIA-LEVEL DATA AND CREATE ADDITIONAL COLUMNS
m_epa_eia_col = m_epa_eia.copy()
m_epa_eia_col['gen_in_epa'] = m_epa_eia_col.key_unit.notna()
m_epa_eia_col = (m_epa_eia_col.groupby(
    ['year', 'utility_id', 'plant_code', 'generator_id', 'key_gen', 'eia_epa_merge_type', 'state_plant', 'energy_source_simp', 
     'entity_type', 'status', 'status_simp', 'dt_operation_start', 'nameplate_capacity_mw'])['gen_in_epa'].any()
     .reset_index())
m_epa_eia_col['plant_in_epa'] = m_epa_eia_col.groupby(['year', 'plant_code'])['gen_in_epa'].transform('any')
m_epa_eia_col['not_renewable'] = m_epa_eia_col.energy_source_simp != 'renewables'
m_epa_eia_col['plant_has_any_ff'] = m_epa_eia_col.groupby(['plant_code', 'year'])['not_renewable'].transform('any')
m_epa_eia_col['gen_age_yrs'] = (dt.datetime.today() - m_epa_eia_col.dt_operation_start).dt.days / 365.25

# %%
# 2. TODO: drop relevant statuses (e.g., discontinued)
# FINDING: LOOKS LIKE STATUS COLUMN NEEDS TO BE CLEANED BETTER
display(m_epa_eia_col.groupby(['status_simp', 'gen_in_epa', 'eia_epa_merge_type'])['key_gen'].nunique())

# %%
# 3. Summarize at the plant level
# drop plants that have no fossil fuels and all plants not operating or standby
print(m_epa_eia_col.groupby(['plant_has_any_ff', 'plant_in_epa'])['plant_code'].nunique())
eia_col_plant = m_epa_eia_col.loc[m_epa_eia_col.plant_has_any_ff & 
                                  m_epa_eia_col.status_simp.isin(['operating', 'standby'])]
# summarize
plant_summ = pd.DataFrame([])
groups = [None, 'entity_type']
for g in groups:
    plant_summ = pd.concat([summarize_plant_inepa(eia_col_plant, g), plant_summ])
plant_summ

plant_summ.loc[plant_summ.value == 'total'].round(2)


# %%
# 4. SUMMARIZE AT THE GENERATOR LEVEL
# drop nonrenewables
print('Dropping nonrenewable generators:\n', m_epa_eia_col.groupby('not_renewable')['key_gen'].nunique())
eia_col_gen = m_epa_eia_col.loc[m_epa_eia_col.not_renewable & 
                                  m_epa_eia_col.status_simp.isin(['operating', 'standby'])]
# summarize
gen_summ = pd.DataFrame([])
groups = [None, 'entity_type']
for g in groups:
    gen_summ = pd.concat([summarize_gen_inepa(eia_col_gen,g), gen_summ])
gen_summ

gen_summ.loc[gen_summ.value == 'total'].round(1)


# %%
# try plotting this:
import seaborn as sns
df, in_epa_lbl = plant_summ, 'plant_in_epa'
# df, in_epa_lbl = gen_summ, 'gen_in_epa'
df = df.loc[gen_summ.value == 'total',].drop(columns='category')
df.iloc[:, 3:] = df.iloc[:, 3:].apply(lambda x: x.round(3).astype(np.float64))

col_val = df.columns[3:]
# Create subplots
plt.figure(figsize=(10, 10), layout='tight')
for i, col in enumerate(col_val):
    ax = plt.subplot(3, 3, i+1)
    sns.lineplot(data=df, x='year', y=col, hue=in_epa_lbl, ax=ax)


# %%
# FIND HIGHEST EMITTING FACILITY AND PLOT DATA
# find facility
YEAR = 2021
epa_top_facility = epa.loc[epa.year==YEAR].groupby(['facility_id'])['co2_mass_short_tons'].sum().sort_values(ascending=False).index[0]
eia_top_facility = xwalk.loc[xwalk.camd_plant_id == epa_top_facility, 'eia_plant_id'].values[0]
epa_topfac_state = epa.loc[(epa.year==YEAR) & (epa.facility_id == epa_top_facility), 'state'].values[0].lower()
# get EPA emissions data
epa_topfac = pd.read_csv(PATH_DATA + f'raw/epa/emissions/daily/emissions-daily-{YEAR}-{epa_topfac_state}.csv')
epa_topfac = epa_topfac.loc[epa_topfac['Facility ID'] == epa_top_facility]
epa_topfac['datetime'] = pd.to_datetime(epa_topfac.Date)
epa_topfac = epa_topfac.groupby('datetime')['CO2 Mass (short tons)'].sum().reset_index()
# get EIA emissions data
eia_topfac = pd.read_parquet(PATH_PROCESSED + 'eia_emissions.parquet')
eia_topfac = eia_topfac.loc[(eia_topfac.plant_code == eia_top_facility) & (eia_topfac.year == YEAR), 'tons_of_co2_emissions'].sum()
# get EIA plant details
print(f'Highest emitting facility details and emissions (fac={eia_top_facility}, yr={YEAR})')
display(eia.loc[(eia.plant_code == eia_top_facility) & (eia.year == YEAR)])

plt.figure(figsize=(15, 2))
plt.plot(epa_topfac.datetime, epa_topfac['CO2 Mass (short tons)'])
plt.title(f'EPA emissions for {YEAR} top-emitting facility\nEPA total={epa_topfac["CO2 Mass (short tons)"].sum():0,}; EIA total={eia_topfac:0,}')
plt.show()
# %%
