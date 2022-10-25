# %%
import scanpy as sc
import pandas as pd
import decoupler as dc
import numpy as np
import logging

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# %%
# Define input and output files
if 'snakemake' in locals():
    data_fp = snakemake.input['data']
    act_fp = snakemake.input['activities']
    output_fp = snakemake.output[0]

    # parameters for decoupler run
    network = snakemake.wildcards.network
    conf = snakemake.params[0]
else:
    network = 'pathways'

    data_fp = '../../../results/structural/plate/merged.h5ad'
    act_fp = '../../../results/functional/plate/activities_{0}.csv'.format(network)

    # parameters for decoupler run
    conf = {'normalisation': 'log1p', 'top_targets': 300, 'method': 'mlm'}

# %%
if network != 'pathways' and conf.get('plot_top') is None:
    logging.warning('Missing \'plot_top\' in config file for {0}. Plotting top 10 most variable across samples'.format(network))
elif network == 'pathways' and conf.get('plot_top') is None:
    logging.warning('Missing \'plot_top\' in config file for {0}. Plotting all 14 Progeny pathways'.format(network))
    conf['plot_top'] = 14

# %%
#Load the data
adata = sc.read_h5ad(data_fp)
print(adata)

#Load previously computed activities
acts = pd.read_csv(act_fp, sep=',', index_col=0)
source_names = acts.columns + '_act'
acts.columns = source_names 
acts.head()

# %%
# Compute max and min values for each source (for plotting color bars)
lims = pd.DataFrame({ 'llim' : [np.min(acts.iloc[:,ii]) for ii in range(acts.shape[1])], 'ulim': [np.max(acts.iloc[:,ii]) for ii in range(acts.shape[1])]}, index = acts.columns)
lims['lim'] = [np.max(abs(acts.iloc[:,ii])) for ii in range(acts.shape[1])]
# print('Max and min values per source')
# print(lims)

# %%
# Add activities to the .obs
temp = adata.obs
temp = pd.merge(
    temp,
    acts,
    how="left",
    left_index=True,
    right_index=True,
    sort=False,
)
adata.obs = temp.loc[adata.obs.index]
# adata

# %%
if conf.get('top_targets') is not None:
    title_suf = ' top ' + str(conf.get('top_targets'))
else:
    title_suf = ''

# %%
def summarize_acts(acts, groupby, obs=None, mode='mean', top=10):
    """
    Summarizes activities obtained per group by their mean or median and removes features that do not change across samples.

    Parameters
    ----------
    acts : AnnData or DataFrame
        Activities obtained after running a method.
    groupby : str
        Column name of obs to use for grouping.
    obs : DataFrame
        None or a data-frame with sample meta-data.
    mode : str
        Wheter to use mean or median to summarize.
    min_std : float
        Minimum std to filter out features. Only features with enough variability will be returned. Decrease it to return more
        features.

    Returns
    -------
    summary : DataFrame
        Dataframe with summaried actvities per group.
    """

    # Extract acts, obs and features
    # if type(acts) is AnnData:
    #     if obs is not None:
    #         raise ValueError('If acts is AnnData, obs needs to be None.')
    #     obs = acts.obs[groupby].values.astype('U')
    #     features = acts.var.index.values.astype('U')
    #     acts = acts.X
    # else:
    #     obs = obs[groupby].values.astype('U')
    #     features = acts.columns.astype('U')
    #     acts = acts.values
    obs = obs[groupby].values.astype('U')
    features = acts.columns.astype('U')
    acts = acts.values

    # Get sizes
    groups = np.unique(obs)
    n_groups = len(groups)
    n_features = acts.shape[1]

    # Init empty mat
    summary = np.zeros((n_groups, n_features), dtype=np.float32)

    for i in range(n_groups):
        msk = obs == groups[i]
        f_msk = np.isfinite(acts[msk])
        for i in range(n_groups):
            msk = obs == groups[i]
            grp_acts = acts[msk]
            for j in range(n_features):
                ftr_acts = grp_acts[:, j]
                f_msk = np.isfinite(ftr_acts)
                if mode == 'mean':
                    summary[i, j] = np.mean(ftr_acts[f_msk])
                elif mode == 'median':
                    summary[i, j] = np.median(ftr_acts[f_msk])
                else:
                    raise ValueError('mode can only be either mean or median.')

    # Filter by min_std
    std = np.std(summary, axis=0, ddof=1)
    if top < std.shape[0]:
        msk = std >= np.sort(std)[-top]
    else:
        msk = std >= -1

    # Transform to df
    summary = pd.DataFrame(summary[:, msk], columns=features[msk], index=groups)

    return summary

# %%
top_sources = summarize_acts(adata.obs.filter(source_names, axis =1), groupby = 'Sample', obs = adata.obs.drop(source_names, axis =1), mode = 'median', top = conf.get('plot_top', 10)).columns

# %%
with PdfPages(output_fp) as output_pdf:
    for source in top_sources:
        fig, axs = plt.subplots(1, 4, figsize=(22.5, 5.5))
        axs = axs.flatten()

        for i, library in enumerate(sorted(adata.obs['plate'].unique())):
            ad = adata[adata.obs.plate == library, :]#.copy()
            sc.pl.spatial(
                ad,
                img_key='lowres',
                library_id='HC-' + library,
                color=source,
                size=1.5,
                legend_loc=None,
                show=False,
                vmin = (lims.loc[source, 'llim']*1.1),
                vmax = (lims.loc[source, 'ulim']*1.1),
                color_map = 'BrBG',
                vcenter = 0,
                ax=axs[i],
            )
            axs[i].set_title(ad.obs['plate'][0])
            axs[i].set_facecolor('#D9D9D9')
            axs[i].set_ylabel('')
            axs[i].set_xlabel('')

        plt.suptitle(source[:-4] + title_suf, fontsize = 15)
        plt.tight_layout()
        output_pdf.savefig(fig)
        plt.close()
        