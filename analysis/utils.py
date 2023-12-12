import pandas as pd


# HELPER FUNCTIONS
def format_sumpct_col(grouper, col, params):
    summ = grouper[[col]].agg('sum') / params['denom']
    cat = grouper.grouper.names[:-1]
    pct = (
        summ[col] / 
        summ.reset_index().groupby(cat, dropna=False)[col].transform('sum').values * 100 )
    out = pd.Series(
        summ[col].round(params['round']).astype(str).str.replace('nan', '-') + ' (' + 
        pct.round(params['round']).astype(str).str.replace('nan', '-') + '%)', name=f'{col}_sum_pct')
    return out.to_frame()

def format_npct_col(grouper, col, params):
    summ = grouper[[col]].agg('nunique') / params['denom']
    cat = grouper.grouper.names[:-1]
    pct = (
        summ[col] / 
        summ.reset_index().groupby(cat, dropna=False)[col].transform('sum').values * 100 )
    out = pd.Series(
        summ[col].astype(int).astype(str).str.replace('nan', '-') + ' (' + 
        pct.round(params['round']).astype(str).str.replace('nan', '-') + '%)', name=f'{col}_n_pct')
    return out.to_frame()

def format_mstd_col(grouper, col, params):
    summ_m = grouper[[col]].mean() / params['denom']
    summ_s = grouper[[col]].std() / params['denom']
    out = pd.Series(
        summ_m[col].round(params['round']).astype(str).str.replace('nan', '-') + '±' + 
        summ_s[col].round(params['round']).astype(str).str.replace('nan', '-'), name=f'{col}_mean_std')
    return out.to_frame()

def format_iqr_col(grouper, col, params):
    summ_l = grouper[[col]].quantile(0.25) / params['denom']
    summ_h = grouper[[col]].quantile(0.75) / params['denom']
    out = pd.Series(
        '[' + summ_l[col].round(params['round']).astype(str).str.replace('nan', '') + '-' + 
        summ_h[col].round(params['round']).astype(str).str.replace('nan', '') + ']', name=f'{col}_iqr')
    return out.to_frame()


# SAMPLE SUMMARY
def sample_summ(df, cat, summdict):
    grouper = df.groupby(cat, dropna=False)
    summ_cat = pd.DataFrame()
    for k in summdict.keys():
        params = summdict[k]['params']
        for aggfn in summdict[k]['aggfns']:
            sdf = aggfn(grouper, k, params)
            summ_cat = pd.concat([sdf, summ_cat], axis=1)
    return summ_cat
