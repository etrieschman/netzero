import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# HELPER FUNCTIONS
def format_sumpct_col(grouper, col, params):
    summ = grouper[[col]].agg('sum') / params['denom']
    cat = grouper.grouper.names[:-1]
    pct = (
        summ[col] / 
        summ.reset_index().groupby(cat, dropna=False)[col].transform('sum').values * 100 )
    out = pd.Series(
        summ[col].round(params['round']).astype(str).str.replace('nan', '-') + ' (' + 
        pct.round(params['round']).astype(str).str.replace('nan', '-') + '%)', 
        name=(col, 'sum_pct'))
    return out.to_frame()

def format_npct_col(grouper, col, params):
    summ = grouper[[col]].agg('nunique') / params['denom']
    cat = grouper.grouper.names[:-1]
    pct = (
        summ[col] / 
        summ.reset_index().groupby(cat, dropna=False)[col].transform('sum').values * 100 )
    out = pd.Series(
        summ[col].astype(int).astype(str).str.replace('nan', '-') + ' (' + 
        pct.round(params['round']).astype(str).str.replace('nan', '-') + '%)', 
        name=(col, 'n_pct'))
    return out.to_frame()

def format_mstd_col(grouper, col, params):
    summ_m = grouper[[col]].mean() / params['denom']
    summ_s = grouper[[col]].std() / params['denom']
    out = pd.Series(
        summ_m[col].round(params['round']).astype(str).str.replace('nan', '-') + '±' + 
        summ_s[col].round(params['round']).astype(str).str.replace('nan', '-'), 
        name=(col, 'mean_std'))
    return out.to_frame()

def format_iqr_col(grouper, col, params):
    summ_l = grouper[[col]].quantile(0.25) / params['denom']
    summ_h = grouper[[col]].quantile(0.75) / params['denom']
    out = pd.Series(
        '[' + summ_l[col].round(params['round']).astype(str).str.replace('nan', '') + '-' + 
        summ_h[col].round(params['round']).astype(str).str.replace('nan', '') + ']', 
        name=(col, 'iqr'))
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
    summ_cat.columns = pd.MultiIndex.from_tuples(summ_cat.columns)
    return summ_cat


# PLOT GENERATION AND CAPACITY MIX
options = {
    'nameplate_capacity_mw':{'units':'MW', 'denom':1},
    'net_gen_tot_an':{'units':'TWh', 'denom':1e6}
}

def plot_mix(df, var, col_loc, path, denom, units, 
             show_pct=True, ncols=3, scale=[5,4]):
    gm = (df
        .groupby(['year', 'energy_source_1_subcat', col_loc])
        .agg({var:'sum'}) / denom).reset_index()
    gm[f'{var}_loc'] = gm.groupby(['year', col_loc])[var].transform('sum')
    gm[f'{var}_pct'] = gm[var] / gm[f'{var}_loc']

    # Pivot the data to get energy type breakdown for each state by year
    var_summ = var if not show_pct else f'{var}_pct'
    var_summ_unit = 'pct' if show_pct else units
    gmt = gm.pivot_table(
        values=var_summ, index=['year', col_loc], 
        columns='energy_source_1_subcat', aggfunc='sum')
    # Fill missing values with 0 (if any)
    gmt = gmt.fillna(0).reset_index()
    
    # Now we plot a stacked bar chart for each state
    locs = gm[col_loc].unique()

    # Set up the matplotlib figure and axes
    nrows = int(np.ceil(len(locs)/ncols))
    fig, axes = plt.subplots(
        ncols=ncols, nrows=nrows,
        figsize=(scale[0]*ncols, scale[1]*nrows), sharex=False, sharey=True)
    # Loop through each location and create a stacked bar chart
    for ax, loc in zip(axes.flatten(), locs):
        loc_data = gmt.loc[gmt[col_loc] == loc]
        loc_data.plot(kind='bar', x='year', legend=False, stacked=True, cmap='tab20', 
                      ax=ax)
        ax.set_title(f'{var} in {loc}', fontsize=20)
        ax.set_xlabel('Year', fontsize=14)
        ax.set_ylabel(f'{var_summ_unit}', fontsize=14)
        ax.tick_params(axis='both', labelsize=12)
    axes[0,-1].legend(title='Energy Type', 
                      loc='upper left', bbox_to_anchor=(1, 1),
                      fontsize=14)

    # Adjust the layout to prevent overlap
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.show()
