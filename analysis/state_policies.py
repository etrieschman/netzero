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

# %%
# RENEWABLE PORTFOLIO STANDARDS
# I DID THIS ONE MANUALLY
from bs4 import BeautifulSoup
import requests
from tqdm import tqdm
import re
url = 'https://www.ncsl.org/energy/state-renewable-portfolio-standards-and-goals/maptype/'
req = requests.get(url)
soup = BeautifulSoup(req.content, 'html.parser')

# %%
# Initialize an empty list to store your data
data = []

# Loop through each state
for state_tag in tqdm(soup.find_all('h3')):
    state_name = state_tag.get_text()
    if state_name == 'Puerto Rico':
        break
    ul = state_tag.find_next_sibling('ul')
    lis = ul.findAll('li')
    for li in lis:
        # Extract each piece of information from the <li> tag
        text = li.get_text()
        label, value = text.split(':', 1)
        data += [{
            'state': state_name,
            'label': label.strip(),
            'value': value.strip()
        }]

# %%
rps = pd.DataFrame(data)
rps = rps.pivot(index='state', columns='label', values='value').reset_index()
rps.columns = rps.columns.str.lower().str.replace(' ', '_')

# Clean data
# 1. Create categorical columns for applicable sectors
sectors = ['investor-owned utility', 'municipal utilities', 'cooperative utilities', 'local government', 'retail supplier']
for sector in sectors:
    rps[sector] = rps.applicable_sectors.str.lower().str.contains(sector)
rps.drop(columns='applicable_sectors', inplace=True)
rps.to_csv(temppath + 'spd_rps.csv', index=False)

# 2. Parse targets into long format
# cool, but not worth it
# rps['requirements'] = rps['requirement'].apply(lambda x: re.findall(r'(\d+)%\s*(.*?)\s*(by|from)\s*(.*?)\s*(\d+)', x))
# reqs = []
# for __, row in rps.iterrows():
#     for target, year in row['requirements']:
#         reqs.append({'state': row.state, 'target': target, 'year': year})
# reqs = pd.DataFrame(reqs)
# rpsl = pd.merge(left=rps.drop(columns='requirements'), 
#                right=reqs, how='outer', on=['state'])


# %%
rps
# %%
