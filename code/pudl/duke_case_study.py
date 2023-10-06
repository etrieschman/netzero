# %%
import sqlite3
import pandas as pd

pd.set_option('display.max_columns', None)

datastore_folder = 'data_pudl'
datastore_path = f'../../{datastore_folder}'

# %%
# GET TABLE NAMES
conn = sqlite3.connect(f'{datastore_path}/sqlite/pudl.sqlite')
cursor = conn.cursor()
cursor.execute("SELECT * FROM sqlite_master WHERE type == 'table'")
table_names = pd.DataFrame(cursor.fetchall(), columns=[column[0] for column in cursor.description])
display(table_names)

# %% 
# IDENTIFY DUKE ENERGY UTILITIES
# xx_entity_eia: clean data across time
# xx_eia: all ids across time
# xx_eia860: data by report_date
table = 'utilities_entity_eia'
cursor.execute(f'select * from {table}')
df = pd.DataFrame(cursor.fetchall(), columns=[column[0] for column in cursor.description])
duke_utilities = df.loc[df.utility_name_eia.str.lower().str.contains('duke')]
print('Unique utilities:', len(duke_utilities.utility_id_eia.drop_duplicates()))
# %%
# IDENTIFY DUKE PLANTS
table = 'plants_eia860'
cursor.execute(f'select * from {table}')
df = pd.DataFrame(cursor.fetchall(), columns=[column[0] for column in cursor.description])
duke_plants = df.loc[df.utility_id_eia.isin(duke_utilities.utility_id_eia.values)].copy()
duke_plants['year'] = pd.to_datetime(duke_plants.report_date).dt.year
print('Unique utilities:', len(duke_plants.utility_id_eia.drop_duplicates()))
print('Unique plants:', len(duke_plants.plant_id_eia.drop_duplicates()))
# %%
# GET DUKE GENERATORS
table = 'generators_eia860'
cursor.execute(f'select * from {table}')
df = pd.DataFrame(cursor.fetchall(), columns=[column[0] for column in cursor.description])
duke_generators = df.loc[df.utility_id_eia.isin(duke_utilities.utility_id_eia.values)].copy()
duke_generators['year'] = pd.to_datetime(duke_generators.report_date).dt.year
print('Unique utilities:', len(duke_generators.utility_id_eia.drop_duplicates()))
print('Unique plants:', len(duke_generators.plant_id_eia.drop_duplicates()))
print('Unique generators:', len(duke_generators[['plant_id_eia', 'generator_id']].drop_duplicates()))

# %%
# PULL GENERATION AND HEAT RATES
table = 'generation_fuel_eia923'
cursor.execute(f'select * from {table}')
df = pd.DataFrame(cursor.fetchall(), columns=[column[0] for column in cursor.description])
duke_pgen = df.loc[df.plant_id_eia.isin(duke_plants.plant_id_eia)]
print(duke_pgen.groupby('report_date')['net_generation_mwh'].sum())

table = 'generation_eia923'
cursor.execute(f'select * from {table}')
df = pd.DataFrame(cursor.fetchall(), columns=[column[0] for column in cursor.description])
duke_ggen = df.loc[df.plant_id_eia.isin(duke_plants.plant_id_eia)]
print(duke_ggen.groupby('report_date')['net_generation_mwh'].sum())




# 'generation_fuel_eia923'

# %%
conn.close()
