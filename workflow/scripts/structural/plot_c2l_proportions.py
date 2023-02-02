# %%
import scanpy as sc
import pandas as pd
import decoupler as dc
import numpy as np
import os
import re

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# %%
import seaborn as sns

# %%
# Define input and output files
if 'snakemake' in locals():

    cellprop_cutoff = snakemake.params[0]

    adata_fp = snakemake.input[0]
    ctProp_fp = snakemake.input[1]

    output_fp = snakemake.output[0]

else:

    cellprop_cutoff = 0.03

    adata_fp = 'results/structural/plate/merged.h5ad'
    ctProp_fp = 'results/Misty/cellprop.csv'


# %%
adata = sc.read_h5ad(adata_fp)

# %%
#count number of visium spots per core
spots_per_core = adata.obs.filter(['Sample', 'orig.ident'], axis = 1).groupby(['Sample']).count()
spots_per_core.columns = ['n']

# %%
ctProp = pd.read_csv(ctProp_fp, index_col=0)
ctProp = ctProp.where(ctProp >= cellprop_cutoff)
ctProp = pd.merge(ctProp, adata.obs.filter(['Sample'], axis = 1), left_index=True, right_index=True)

# %%
#compute number of spots w. cell type
counts = ctProp.groupby(['Sample']).count().reset_index()
counts.columns = ['Sample'] + ctProp.columns[0:-1].tolist()

#divide number of spots by total in core
spot_prop = counts.iloc[:,1:].divide(spots_per_core['n'].to_numpy(), axis=0)
spot_prop['Sample'] = counts['Sample']

#go from wide to long, and merge together
counts = pd.melt(counts, id_vars='Sample', var_name='celltype', value_name='count')
spot_prop = pd.melt(spot_prop, id_vars='Sample', var_name='celltype', value_name='prop')
counts = pd.merge(counts, spot_prop, on = ['Sample', 'celltype'])
counts = pd.merge(adata.obs.filter(['Sample', 'Confidence', 'PFI', 'patient'], axis=1).drop_duplicates().reset_index(drop = True), counts, on='Sample')

# %%
medians = ctProp.groupby(['Sample']).median().reset_index()
medians.columns = ['Sample'] + ctProp.columns[0:-1].tolist()
medians = pd.melt(medians, id_vars='Sample', var_name='celltype', value_name='median')
medians = pd.merge(adata.obs.filter(['Sample', 'Confidence', 'PFI', 'patient'], axis=1).drop_duplicates().reset_index(drop = True), medians, on='Sample')

# %%

with PdfPages(output_fp) as output_pdf:

    ##plot on high confidence cores only
    counts1 = counts.loc[counts['Confidence'] == 'High confidence']
    medians2 = medians.loc[medians['Confidence'] == 'High confidence']

    colors = sns.color_palette("tab10").as_hex()[0:9]

    fig, axes = plt.subplots(1, 2 , figsize=(16, 15), dpi = 300)
    axes = axes.flatten()


    for df, name, ax in zip([counts1, medians2], ['prop', 'median'], axes):

        sns.boxplot(data = df, x=name, y= 'celltype', hue = 'PFI', fliersize=0, ax=ax, palette='Pastel1')

        # for mouse, color in zip(sorted(df['patient'].unique()), colors):
        #     mouse_toplot = df[df['patient'] == mouse]
        #     sns.stripplot(x=name, y='celltype', hue='PFI', dodge=True, palette=[color] * 2, marker='o', data=mouse_toplot, ax=ax)
        sns.stripplot(x=name, y='celltype', hue='PFI', dodge=True, color = 'black', marker='o', data=df, ax=ax)

        handles, labels = ax.get_legend_handles_labels()
        handles = handles[1:3] #+ handles[3:-1:3]
        labels = labels[1:3] #+ sorted(counts['patient'].unique())
        ax.legend(handles, labels)

    axes[0].set_title('Proportion of spots in core containing celltype')
    axes[1].set_title('Median proportions in spots')
    axes[0].set_ylabel('')
    axes[0].set_xlabel('% of spots in core')
    axes[1].set_ylabel('')
    axes[0].set_xlim(0, np.max(counts['prop']) + 0.1)
    axes[1].set_xlim(0, np.max(medians['median']) + 0.05)

    plt.suptitle('Cell proportions from cell2location in high confidence cores')
    plt.tight_layout()
    output_pdf.savefig(fig)
    plt.close()


    ## plot on all cores
    colors = sns.color_palette("tab10").as_hex()[0:3]

    fig, axes = plt.subplots(1, 2 , figsize=(16, 15), dpi = 300)
    axes = axes.flatten()


    for df, name, ax in zip([counts, medians], ['prop', 'median'], axes):

        sns.boxplot(data = df, x=name, y= 'celltype', hue = 'PFI', fliersize=0, ax=ax, palette='Pastel1')

        for conf, color in zip(sorted(df['Confidence'].unique()), colors):
            conf_toplot = df[df['Confidence'] == conf]
            sns.stripplot(x=name, y='celltype', hue='PFI', color = 'Confidence', dodge=True, palette=[color] * 2, marker='o', data=conf_toplot, ax=ax)
        # sns.stripplot(x=name, y='celltype', hue='PFI', dodge=True, color = 'black', marker='o', data=df, ax=ax)

        handles, labels = ax.get_legend_handles_labels()
        handles = handles[0:3] + handles[3:-1:3]
        labels = labels[0:3] + sorted(df['Confidence'].unique())
        ax.legend(handles, labels)

    axes[0].set_title('Proportion of spots in core containing celltype')
    axes[1].set_title('Median proportions in spots')
    axes[0].set_ylabel('')
    axes[0].set_xlabel('% of spots in core')
    axes[1].set_ylabel('')
    axes[0].set_xlim(0, np.max(counts['prop']) + 0.1)
    axes[1].set_xlim(0, np.max(medians['median']) + 0.05)

    plt.suptitle('Cell proportions from cell2location')
    plt.tight_layout()
    output_pdf.savefig(fig)
    plt.close()
