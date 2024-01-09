# %%
import pandas as pd
import numpy as np

# import os
# os.environ['GS_LIB'] = '/usr/local/Cellar/ghostscript/10.02.1/lib/libgs.dylib'

temppath = '../data/temp/'

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)


# %%
# CARBON PRICE
# READIN
cp = pd.read_excel(temppath + 'State Policy Drivers_Detail.xlsx', sheet_name='Carbon price', nrows=19)
assert cp.Year.iloc[-1] == 2023
cp.columns = cp.columns.str.lower().str.replace(' ', '_')

# RESHAPE
cpl = cp.melt(id_vars='year')
cpl.loc[:, 'state'] = cpl.variable.str.split('_').str.get(0)
cpl['variable'] = cpl.variable.str.split('_').str[1:].apply('_'.join)
cplw = cpl.pivot(index=['year', 'state'], columns=['variable'], values='value').reset_index()

# REFORMAT
cplw['state'] = cplw.state.str.upper()
cplw['has_carbon_price'] = cplw.has_carbon_price == 'Y'
cplw.loc[~cplw.has_carbon_price, 'price'] = pd.NA

# SAVE
cplw.to_csv(temppath + 'state_policy_drivers_carbonprice.csv', index=False)

# %%
# GHG TARGET
# READIN
ghg = pd.read_excel(temppath + 'State Policy Drivers_Detail.xlsx', sheet_name='GHG Target', nrows=19)
assert ghg.Year.iloc[-1] == 2023
ghg.columns = (ghg.columns.str.lower()
               .str.replace(' ', '_').str.replace('__', '_')
               .str.replace('*', '', regex=False))

# RESHAPE
ghgl = ghg.melt(id_vars='year')
ghgl.loc[:, 'state'] = ghgl.variable.str.split('_').str.get(0)
ghgl['variable'] = ghgl.variable.str.split('_').str[1:].apply('_'.join)
ghgl.loc[ghgl.variable == 'ghg_target', 'variable'] = 'has_ghg_target'
ghgl = ghgl.loc[ghgl.variable != 'target_level']
# target year
ghgl.loc[:, 'variable'] = np.where(ghgl.variable != 'has_ghg_target', 
                                   ghgl.variable.str.rsplit('_').str.get(-1), 'has_ghg_target')
# update netzero and carbon neutral
ghgl.loc[ghgl.value == 'NZ', 'value'] = 1.
ghgl.loc[ghgl.value == 'NZ', 'target_detail'] = 'net_zero'
ghgl.loc[ghgl.value == 'CN', 'value'] = 1.
ghgl.loc[ghgl.value == 'CN', 'target_detail'] = 'carbon_neutral'

ghglw = ghgl.pivot(index=['year', 'state'], columns='variable', values='value').reset_index()

# REFORMAT
ghglw['state'] = ghglw.state.str.upper()
ghglw['has_ghg_target'] = ghglw.has_ghg_target == 'Y'
ghglw.to_csv(temppath + 'state_policy_drivers_ghgtarget.csv', index=False)
