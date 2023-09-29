# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import PATH_DATA, PATH_PROCESSED

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 115)

# %%
# IMPORT
eia = pd.read_parquet(PATH_PROCESSED + 'eia_final.parquet')
epa = pd.read_csv(PATH_PROCESSED + 'epa_emissions.csv')
xwalk = pd.read_csv(PATH_DATA + 'epa_eia_crosswalk.csv')
eia['key_gen'] = eia.plant_code.astype(str) + '_' + eia.generator_id
eia['key_own'] = eia.key_gen + '_' + eia.ownership_id.astype(str)

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
# SUMMARIZE UNIT-TO-GENERATOR MAPPINGS
xwalk_summ = xwalk.copy()
xwalk_summ['n_eia_plants_pepa'] = xwalk.groupby(['camd_plant_id'])['eia_plant_id'].transform('nunique')
xwalk_summ['n_epa_plants_peia'] = xwalk.groupby(['eia_plant_id'])['camd_plant_id'].transform('nunique')
display(xwalk_summ[['camd_plant_id', 'eia_plant_id', 'n_eia_plants_pepa']].drop_duplicates()['n_eia_plants_pepa'].value_counts())
display(xwalk_summ[['camd_plant_id', 'eia_plant_id', 'n_epa_plants_peia']].drop_duplicates()['n_epa_plants_peia'].value_counts())
# NOTE: Anecdotally, we get multiple EIA plants to an EPA plant when an old plant goes offline and is replaced by a new plant or new owner
# display(xwalk.loc[xwalk.groupby(['camd_plant_id'])['eia_plant_id'].transform('nunique') > 1, ['camd_plant_id', 'eia_plant_id']].drop_duplicates())
# eia.loc[eia.plant_code.isin([302, 59002])]
# epa.loc[epa.facility_id == 302]
xwalk_summ['n_eia_gen_per_epa_unit'] = xwalk.groupby(['camd_plant_id', 'camd_unit_id'])['eia_generator_id'].transform('nunique')
xwalk_summ['n_epa_unit_per_eia_gen'] = xwalk.groupby(['eia_plant_id', 'eia_generator_id'])['camd_unit_id'].transform('nunique')
xwalk_gensum = xwalk_summ[['n_eia_gen_per_epa_unit', 'n_epa_unit_per_eia_gen']].value_counts().reset_index()
display(xwalk_gensum.iloc[0:16])
display(xwalk_gensum.iloc[16:])
# DECISION: BASED ON THESE FINDINGS, PROPOSE MAPPING EPA TO EIA AT THE FACILITY LEVEL

# %%
# SUMMARIZE MERGE HELPER FUNCTION
def summarize_merge(df, left_on, right_on, r_hasyear):
    summ = pd.DataFrame([])
    for y in df.year.drop_duplicates().sort_values().dropna().values:
        s = {}
        s['year'] = int(y)
        s['left'] = [len(df.loc[(df.year == y), left_on].drop_duplicates())]
        s['left_nright'] = [len(df.loc[(df.year == y) & df[right_on].isna().all(axis=1), 
                                    left_on].drop_duplicates())]
        s['overlap'] = [len(df.loc[(df.year == y) & df[left_on].notna().all(axis=1) & df[right_on].notna().all(axis=1),
                                    left_on + right_on].drop_duplicates())]
        if r_hasyear:
            y_condition = (df.year == y)
        else:
            y_condition = True
        s['right_nleft'] = [len(df.loc[y_condition & df[left_on].isna().all(axis=1) & df[right_on].notna().all(axis=1), right_on].drop_duplicates())]
        s['right'] = [len(df.loc[y_condition & df[right_on].notna().all(axis=1), right_on].drop_duplicates())]
        summ = pd.concat([summ, pd.DataFrame(s)], ignore_index=True)
    return summ


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
epa_col = epa.groupby(['facility_id', 'unit_id', 'year'])[['co2_mass_short_tons', 'gross_load_mwh']].sum().reset_index()
epa_on = ['facility_id', 'unit_id']
xwalk_epa_on = ['camd_plant_id', 'camd_unit_id']
m_epa = pd.merge(left=epa_col, right=xwalk, how='outer',
              left_on=epa_on, right_on=xwalk_epa_on)
m_epa['year'] = m_epa.year.astype('Int64')
display(summarize_merge(m_epa, epa_on, xwalk_epa_on, r_hasyear=False))

# %%
# FLAG EIA PLANTS THAT HAVE EPA EMISSIONS DATA
# first, add flag to EIA if in EPA data (by year)
epa_flag = m_epa.loc[m_epa.facility_id.notna(), ['facility_id', 'eia_plant_id', 'eia_generator_id', 'year']].drop_duplicates()
epa_flag.rename(columns={'eia_plant_id':'eia_plant_id_in_epa_xwalk'}, inplace=True)
print("Count by year of facilities in EPA that aren't in EIA:")
epa_flag.loc[epa_flag.eia_plant_id_in_epa_xwalk.isna()].groupby('year')['facility_id'].nunique()
# TODO: Look into these facilities that don't appear in EIA. Understand why
epa_flag = epa_flag.drop(columns='facility_id').drop_duplicates()
print('EIA pre-merge shape:', eia.shape)
eia_epaflag = pd.merge(left=eia, right=epa_flag, how='left', left_on=['plant_code', 'generator_id', 'year'], right_on=['eia_plant_id_in_epa_xwalk', 'eia_generator_id', 'year'])
print('EIA post-merge shape:', eia_epaflag.shape)
eia_epaflag['in_epa'] = False
eia_epaflag.loc[eia_epaflag.eia_plant_id_in_epa_xwalk.notna(), 'in_epa'] = True
# next, summarize
# 1. by status
eia_summ_status = eia_epaflag.groupby(['status_simp', 'in_epa']).agg({'key_gen':'nunique', 'nameplate_capacity_mw':'sum'})
eia_summ_status


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