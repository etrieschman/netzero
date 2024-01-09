# %%
# SETUP
import pandas as pd
import numpy as np
from tqdm import tqdm
import dropbox
import io
import sys, os
import openpyxl

parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, parent_dir)
from credentials import DB_ACCESS_TOKEN

# options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)

# dropbox access
dbx = dropbox.Dropbox(DB_ACCESS_TOKEN)
rps_path = '/Net Zero Cluster/data/State Policy/rps/data_raw/'

# %%
# FROM JOE-ALDY PROVIDED RPS ARCHIVES: 2009-2012
_, j = dbx.files_download(rps_path + 'aldy_state_rps_laws_and_data.xlsx')
j = openpyxl.load_workbook(io.BytesIO(j.content), data_only=False)
sheet_names = j.sheetnames
sheet_names
# %%
# FROM ARCHIVED DATA IN RPS ARCHIVES: 2015-2023
# get relevant archives

entries = dbx.files_list_folder(rps_path)
archives = [f.name for f in entries.entries if f.name.startswith('dsire-')]

# For each archive, download, transform, and stack
# Use the Dropbox API to stream the file's content
cols_keep = ['program_id', 'state_id', 'code', 'name', 'websiteurl', 'summary']
def get_archive(rps_path, archive, cols_keep=cols_keep):
    # read in files
    _, p = dbx.files_download(rps_path + archive + '/program.csv')
    _, t = dbx.files_download(rps_path + archive + '/program_detail.csv')
    _, s = dbx.files_download(rps_path + archive + '/state.csv')
    p = pd.read_csv(io.BytesIO(p.content))
    t = pd.read_csv(io.BytesIO(t.content))
    s = pd.read_csv(io.BytesIO(s.content))

    # clean files
    t['label'] = t.label.str.lower().str.replace(r'\s|//', '_', regex=True)
    t = (t
         .drop(columns=['id', 'display_order', 'template_id'],)
         .loc[~t.label.isin(['alternative_compliance_payment', 'compliance_multipliers', 
              'credit_trading/tracking_system'])])
    p = (p
        .rename(columns={'id':'program_id'})
        .loc[(p.program_type_id == 38), cols_keep])
    p['summary'] = (p.summary
        .str.replace(r'<([^>]*)>', '', regex=True)
        .str.replace(r'\\[a-zA-Z]', '', regex=True)
        .str.replace(r'&#\d+', '', regex=True))
    s.rename(columns={'id':'state_id', 'name':'state_name'}, inplace=True)
    # merge files
    pt = pd.merge(left=p, right=t, how='left', on='program_id')
    pt = pt.pivot(
        index=cols_keep, values='value', columns='label').reset_index()
    pts = pd.merge(left=pt, right=s, how='inner', on='state_id')
    pts = pts.loc[(pts.is_territory == 0) | (pts.abbreviation == 'DC')]
    pts.drop(columns='is_territory', inplace=True)
    return pts

dfs = pd.DataFrame()
for archive in tqdm(archives):
    df = get_archive(rps_path, archive)
    df['date'] = pd.to_datetime(archive[6:])
    dfs = pd.concat([df, dfs], ignore_index=True)
id_cols = ['date', 'program_id', 'state_id', 'abbreviation', 'state_name']
dfs = dfs[id_cols + [col for col in dfs.columns if col not in id_cols]]

# upload to dropbox
try:
    dbx.files_upload(
        dfs.to_csv(index=False).encode(), 
        rps_path + f'df_cln_dsire_{dfs.date.dt.year.min()}_{dfs.date.dt.year.max()}.csv',
        mode=dropbox.files.WriteMode.overwrite)
    print("File uploaded successfully")
except dropbox.exceptions.ApiError as err:
    print(f"API Error: {err}")



# %%
# APPENDIX - RPS USING NCSL DATA
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
rps.to_csv('spd_rps_raw.csv', index=False)

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
