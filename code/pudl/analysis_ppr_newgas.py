# %%
import sqlite3
import pandas as pd
pd.set_option('display.max_columns', None)

datastore_folder = 'data_pudl'
datastore_path = f'../../{datastore_folder}'

def get_table(cursor, table_name, where=None):
    if where is None:
        cursor.execute(f'select * from {table_name}')
    else:
        cursor.execute(f'select * from {table_name} where {where}')
    df = pd.DataFrame(cursor.fetchall(), 
                      columns=[column[0] for column in cursor.description])
    return df
# %%
# GET TABLE NAMES
conn = sqlite3.connect(f'{datastore_path}/sqlite/pudl.sqlite')
cursor = conn.cursor()
table_names = get_table(cursor, 'sqlite_master', "type == 'table'").name
# cursor.execute("SELECT name FROM sqlite_master WHERE type == 'table'")
# table_names = pd.DataFrame(cursor.fetchall(), columns=[column[0] for column in cursor.description])
display(table_names)

# %% 
# GET ALL NATURAL GAS GENERATORS
gens = get_table(cursor, 'generators_eia860')
gens = gens.loc[(gens.operational_status == 'existing')] # existing
gens = gens.loc[gens.plant_id_eia.isin(
    gens.loc[gens.fuel_type_code_pudl == 'gas', 'plant_id_eia'])] # at gas plant
gas_plant_ids = gens.plant_id_eia.drop_duplicates().values

# %%
# GET ALL PLANT DATA FOR PLANTS WITH GAS GENERATORS
gas_plant = get_table(cursor, 'plants_eia860')
gas_plant = gas_plant.loc[gas_plant.plant_id_eia.isin(gas_plant_ids)]

# %%
# GET ALL GENERATION FROM THESE PLANTS AND THESE GENERATORS
pwr_gen = get_table(cursor, 'generation_eia923')
pwr_gen = pwr_gen.loc[pwr_gen.plant_id_eia.isin(gas_plant_ids)]
pwr_plant = get_table(cursor, 'generation_fuel_eia923')
pwr_plant = pwr_plant.loc[pwr_plant.plant_id_eia.isin(gas_plant_ids)]

# %%
# MERGE PLANT AND GENERATOR DATA
m = pd.merge(left=gas_plant, right=gens, on=['plant_id_eia', 'utility_id_eia', 'report_date'])
assert len(m) == len(gens)
# %%
# PLANT-LEVEL SUMMARY
# By year
#    - Capacity (% of capacity that is natural gas)
#    - Generation (% of generation from natural gas)
#    - Capacity factor
#    - Emissions
#    - Num generators, average age generators
