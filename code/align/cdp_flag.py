# %%
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from rapidfuzz import fuzz, process

SCORER = fuzz.WRatio

# # global variables
# PATH_DATA = '../data/'
# PATH_PROCESSED = PATH_DATA + 'processed/'
# PATH_INTERIM = PATH_DATA + 'interim/'
# PATH_RESOURCES = PATH_DATA + 'resources/'
# PATH_RESULTS = '../results/analysis/cdp/'
# DENOM, ROUND = 1e6, 2
# rpath = Path(PATH_RESULTS)
# rpath.mkdir(parents=True, exist_ok=True)

# options
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', 500)

# PREPROCESS NAMES
def preprocess_names(x:pd.Series) -> pd.Series:
    # string updates
    x_pp = (x
        .str.lower()
        .str.replace('\.|\,', '', regex=True)
        .str.replace('-', ' ', regex=True)
        .str.replace('\s\/[a-zA-Z]{2}\/?', '', regex=True)
        .str.replace('\\(.*?\\)', '', regex=True)
        .str.strip()
        .str.replace(r'\s\b(corp|company|corporation|corportation|co)\b(\s|$)', ' co ', regex=True).str.strip()
        .str.replace(r'\s\b(elec|electric)\b(\s|$)', ' electric ', regex=True).str.strip()
        .str.replace(r'\s\b(inc|incorporated)\b(\s|$)', ' inc ', regex=True).str.strip()
        .str.replace(r'\s\b(grp|group)\b(\s|$)', ' group ', regex=True).str.strip()
        .str.replace(r'\s\b(int|international)\b(\s|$)', ' int ', regex=True).str.strip()
        .str.replace(r'\s\b(ltd|limited)\b(\s|$)', ' ltd ', regex=True).str.strip()
        .str.replace(r'\b(generatn|gen)\b', ' generation', regex=True).str.strip()
    )
    # one-offs
    x_pp_oo = (x_pp
               .str.replace('pacific gas & electric', 'pg&e', regex=True)
    )
    return x_pp_oo

# FLAG EIA IDENTIFIERS (ID-YEAR) WITH CDP COMPANIES
def fuzzy_search(df_target, df_search, key_target, key_search, 
                 scorer, num_matches):
    # Step 0: fuzzy search
    searches = df_search[key_search].drop_duplicates()
    searches = searches.tolist()
    targets = df_target[key_target].drop_duplicates()
    target_map = targets.apply(lambda x: process.extract(x, searches, scorer=scorer, limit=num_matches))

    # Step 1: Expand the list of tuples into separate columns
    df_expanded = pd.DataFrame(target_map.tolist(), index=targets)
    df_expanded.columns = [f'pair_{i}' for i in range(num_matches)]
    df_expanded.reset_index(inplace=True)
    
    # Step 2: Transform to Long Form
    df_long = df_expanded.melt(id_vars=key_target, var_name='pair', 
                               value_name='name_score', ignore_index=False)
    
    # Step 3: Split the Tuples into Two Columns
    df_long[[key_search, f'{key_search}_score', 'drop']] = pd.DataFrame(
        df_long['name_score'].tolist(), index=df_long.index)
    
    # Drop columns and rows
    df_long = df_long.drop(['name_score', 'pair', 'drop'], axis=1).reset_index(drop=True)
    # df_long_sub = df_long.loc[df_long[f'{key_search}_score'] >= threshold].reset_index()

    return df_long

def fuzzy_search_years(df_target, df_search, key_target, key_search, 
                       scorer=fuzz.WRatio, num_matches=1):
    dfs = pd.DataFrame()
    for y in tqdm(df_target.year.unique()):
        df = fuzzy_search(df_target.loc[df_target.year == y],
                        df_search.loc[df_search.year == y],
                        key_target, key_search, scorer, num_matches)
        df['year'] = y
        dfs = pd.concat([df, dfs], ignore_index=True)
    return dfs


def fuzzy_merge(df_target, key_target, scorer, num_matches):
    # 1. fuzzy merge parent name
    # 1.a. fuzzy merge
    map_parent = fuzzy_search_years(
        df_target, names_cs, key_target, 'pp_cdp_parname',
        scorer, num_matches)
    # 1.b. pull gvkeys onto parnames
    df_parent = pd.merge(
        left=map_parent, how='left',
        right=names_cs[['year', 'pp_cdp_parname', 'gvkey']].drop_duplicates(),
        on=['year', 'pp_cdp_parname']
    )
    assert len(map_parent) == len(df_parent)
    # 2. fuzzy merge subsidiary names
    # 2.a fuzzy merge
    map_subsid = fuzzy_search_years(
        df_target, names_cs, key_target, 'pp_cdp_subname',
        scorer, num_matches)
    # 2.b parnames and gvkeys onto subsid names
    # NOTE: because subsidiaries are reported for 50% ownership, 
    # we expect some dups of magnitute <=2 since we can't have 
    # more than 2 companies that own 50%
    df_subsid = pd.merge(
        left=map_subsid, how='left',
        right=names_cs[['year', 'pp_cdp_parname', 'pp_cdp_subname', 'gvkey']].drop_duplicates(),
        on=['year', 'pp_cdp_subname'])
    # 2.c collapse parent names into multiple columns
    df_subsid = (
        df_subsid
        .groupby(['year', key_target, 'pp_cdp_subname', 'pp_cdp_subname_score'])
        .agg(n_parents=('pp_cdp_parname', 'nunique'),
            pp_cdp_subsid_parnames=('pp_cdp_parname', lambda x: list(set(x))),
            subsid_gvkeys=('gvkey', lambda x: list(set(x))))
        .reset_index()
    )
    assert len(map_subsid) == len(df_subsid)
    # 3. merge flags together
    df = pd.merge(left=df_parent, right=df_subsid, how='outer',
                  on=['year', key_target])
    assert len(df) == len(df_parent)

    return df



# %%
if __name__ == '__main__':
    if "snakemake" not in globals():
        # readin mock snakemake
        import sys, os
        parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        sys.path.insert(0, parent_dir)
        from utils import mock_snakemake
        snakemake = mock_snakemake('cdp_flag')

    # results path
    rpath = Path(snakemake.params.results_dir)
    rpath.mkdir(parents=True, exist_ok=True)
    
    # analysis options
    scorer = SCORER
    threshold = snakemake.params.threshold

    # READ IN DATA
    gdf = pd.read_parquet(snakemake.input.infile_gen)
    pdf = pd.read_parquet(snakemake.input.infile_plant)
    udf = pd.read_parquet(snakemake.input.infile_util)
    odf = pd.read_parquet(snakemake.input.infile_own)
    # cdp = pd.read_csv(PATH_RESOURCES + 'cdp_elec.csv', encoding = 'ISO-8859-1')
    # write out relevant ids
    # conids = cdp.gvkey.drop_duplicates()
    # conids = conids.astype(str).str.zfill(6)
    # conids.to_csv(PATH_RESULTS + 'in_gvkey_list.csv', index=False)
    # READIN PARENT-SUBSIDIARY MAPPING
    # This was pulled from WRDS Company Subsidiary data using gvkeys in CDP dataset
    csids = pd.read_csv(snakemake.input.infile_csids)
    csids['year'] = csids.fdate.str[:4].astype(int)

    # PREPROCESS NAMES
    # subsidiary (includes year)
    names_cs = csids[['gvkey', 'year', 'coname', 'clean_company']].drop_duplicates().reset_index(drop=True)
    names_cs['pp_cdp_parname'] = preprocess_names(names_cs.coname)
    names_cs['pp_cdp_subname'] = preprocess_names(names_cs.clean_company)
    # EIA datasets (includes year)
    udf['pp_utility_name'] = preprocess_names(udf.utility_name)
    uids = udf[['utility_id', 'utility_name', 'pp_utility_name', 'year']].drop_duplicates().reset_index(drop=True)
    odf['pp_owner_name'] = preprocess_names(odf.owner_name)
    oids = odf[['utility_id', 'plant_code', 'generator_id', 
                'ownership_id', 'owner_name', 'pp_owner_name', 'year']].drop_duplicates().reset_index(drop=True)

    # TODO: [start] delete this after testing
    # temp_min_year = 2021
    # uids = uids.loc[uids.year >= temp_min_year]
    # oids = oids.loc[oids.year >= temp_min_year]
    # names_cs = names_cs.loc[names_cs.year >= temp_min_year]
    # df_target, key_target = uids, 'pp_utility_name'
    # df_search, key_search = names_cs, 'pp_cdp_parname'
    # num_matches = 1
    # TODO: [end] delete this after testing
    
    # CREATE MAPPINGS USING FUZZY MERGE
    df_utility = fuzzy_merge(uids, 'pp_utility_name', scorer, num_matches=1)
    df_owner = fuzzy_merge(oids, 'pp_owner_name', scorer, num_matches=1)
    # write interim data to file
    df_utility.to_csv(snakemake.output.intfile_util)
    df_owner.to_csv(snakemake.output.intfile_own)

    # MAKE CDP FLAGS
    # utility
    udf_cdp = pd.merge(left=uids, right=df_utility, how='outer',
                       on=['year', 'pp_utility_name'])
    assert len(udf_cdp) == len(uids)
    udf_cdp['in_cdp'] = ((udf_cdp.pp_cdp_parname_score >= threshold) | 
                         (udf_cdp.pp_cdp_subname_score >= threshold))
    udf_cdp['in_cdp_as_subsid'] = (udf_cdp.pp_cdp_subname_score >= threshold)
    # owner
    odf_cdp = pd.merge(left=oids, right=df_owner, how='outer',
                       on=['year', 'pp_owner_name'])
    assert len(odf_cdp) == len(oids)
    odf_cdp['in_cdp'] = ((odf_cdp.pp_cdp_parname_score >= threshold) | 
                         (odf_cdp.pp_cdp_subname_score >= threshold))
    odf_cdp['in_cdp_as_subsid'] = (odf_cdp.pp_cdp_subname_score >= threshold)
    # collapse on generator ids
    odf_cdp_col = (odf_cdp
        .groupby(['year', 'utility_id', 'plant_code', 'generator_id'])
        .agg(n_parents_own=('n_parents', 'sum'),
            in_cdp_own=('in_cdp', 'max'),
            in_cdp_own_as_subsid=('in_cdp_as_subsid', 'max'),
            ownersip_ids=('ownership_id', lambda x: list(set(x))),
            pp_owner_names=('pp_owner_name', lambda x: list(set(x))),
            pp_cdp_subnames_own=('pp_cdp_subname', lambda x: list(set(x))))
        .reset_index()
    )
    # generator
    gdf_cdp = pd.merge(
        left=gdf[['year', 'utility_id', 'plant_code', 'generator_id']],
        right=udf_cdp[['year', 'utility_id', 'in_cdp', 'in_cdp_as_subsid']], 
        how='left', on=['year', 'utility_id'])
    gdf_cdp = pd.merge(
        left=gdf_cdp, how='left', 
        right=odf_cdp_col[['year', 'utility_id', 'plant_code', 'generator_id',
                           'in_cdp_own', 'in_cdp_own_as_subsid']],
        on=['year', 'utility_id', 'plant_code', 'generator_id']
    )
    # plant
    pdf_cdp = pd.merge(
        left=pdf[['year', 'utility_id', 'plant_code']],
        right=udf_cdp[['year', 'utility_id', 'in_cdp', 'in_cdp_as_subsid']],
        how='left', on=['year', 'utility_id']
    )

    # WRITE TO FILE
    udf_cdp.to_parquet(snakemake.output.outfile_util)
    odf_cdp.to_parquet(snakemake.output.outfile_own)
    pdf_cdp.to_parquet(snakemake.output.outfile_plant)
    gdf_cdp.to_parquet(snakemake.output.outfile_gen)
# %%
