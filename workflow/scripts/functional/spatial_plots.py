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
# Define input and output files
if 'snakemake' in locals():
    adata_fp = snakemake.input[0]
    functional_fp = snakemake.input[1]

    output_fp = snakemake.output[0]

else:

    adata_fp = 'results/structural/plate/merged.h5ad'
    functional_fp = 'results/functional/plate/activities_pathways.csv'



# %%
adata = sc.read_h5ad(adata_fp)
del adata.layers['SCT']
adata

# %%
activities = pd.read_csv(functional_fp, index_col=0, sep=',')
activities = activities.loc[adata.obs.index,:]
activities.columns = [re.sub("-", "", func) for func in activities.columns]

adata.obsm['acts'] = activities
acts = dc.get_acts(adata, 'acts')
print(acts)

# %%
lims = pd.DataFrame({ 'llim' : [np.min(acts.X[:,ii]) for ii in range(acts.n_vars)], 'ulim': [np.max(acts.X[:,ii]) for ii in range(acts.n_vars)]}, index = acts.var_names.values)
lims['lim'] = [np.max(abs(acts.X[:,ii])) for ii in range(acts.n_vars)]
print('Max and min values per pathway')
print(lims)

with PdfPages(output_fp) as output_pdf:
    for source in acts.var.index.values:
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
                # vmin = (lims.loc[source, 'llim']*1.1),
                # vmax = (lims.loc[source, 'ulim']*1.1),
                # vcenter = 0,
                ax=axs[i],
            )
            axs[i].set_title(ad.obs['plate'][0])
            axs[i].set_facecolor('#D9D9D9')
            axs[i].set_ylabel('')
            axs[i].set_xlabel('')

        plt.suptitle(source, fontsize = 15)
        plt.tight_layout()
        output_pdf.savefig(fig)
        plt.close()