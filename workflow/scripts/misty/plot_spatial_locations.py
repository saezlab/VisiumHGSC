# %%
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import string

from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

# %%
# Define input and output files
if 'snakemake' in locals():

    adata_fp = snakemake.input[0]

else:


    adata_fp = 'results/structural/plate/merged.h5ad'


# %%
adata = sc.read_h5ad(adata_fp)

#normalise by total number of cells
adata.obsm['q05_cell_abundance_w_sf'] = adata.obsm['q05_cell_abundance_w_sf'].div(adata.obsm['q05_cell_abundance_w_sf'].sum(axis=1), axis = 0)

adata.obs[[name.translate(str.maketrans('', '', string.punctuation)) for name in adata.uns['mod']['factor_names']]] = adata.obsm['q05_cell_abundance_w_sf']

# %%
# select a sample
samples = ['HC-L_OVA20', 'HC-S_OVA37']
library_ids = ['HC-TMA3', 'HC-TMA2']
sadata = adata[adata.obs['Sample'].isin(samples)].copy()

celltypes = ['EOCC5', 'Macrophages', 'Plasmacells']

# %%
#get max and min and grouped by celltypes in adata.obs
maxs = [None] + list(sadata.obs[celltypes].max())
mins = [None] + list(sadata.obs[celltypes].min())

#get 99th and 1st percentile and grouped by celltypes in adata.obs
qmax = [None] + list(sadata.obs[celltypes].quantile(0.99))
qmin = [None] + list(sadata.obs[celltypes].quantile(0.02))
celltypes = [None] + celltypes

# %%
def add_colorbar(mappable):
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

# %%
# make a figure with six subplots
fig, axs = plt.subplots(2, 4, figsize=(15, 7))
axs = axs.flatten()


for i, ax in enumerate(axs):
    
    mini_adata = adata[adata.obs['Sample'] == samples[int(np.floor(i/4))], :].copy()
    
    sc.pl.spatial(mini_adata, color_map=sns.color_palette("rocket", as_cmap=True),
                  library_id=library_ids[int(np.floor(i/4))],
                  color=celltypes[i%4],
                  size=1.3,
                  img_key='hires',
                  colorbar_loc=None,
                  vmin=qmin[i%4], vmax=qmax[i%4],
                  ax=axs[i],
                  show=False,
                 )
    
    if i%4 != 0:
        add_colorbar(ax.collections[0])
        
    axs[i].set_xlabel('')
    axs[i].set_ylabel('')
    axs[i].set_title('')
    
    if i%4 == 0:
        axs[i].set_ylabel(samples[int(np.floor(i/4))], fontsize=12)
        
    if i == 0:
        axs[i].set_title('H&E', fontsize=12)
    elif i < 4:
        axs[i].set_title(celltypes[i], fontsize=12)

plt.tight_layout()
plt.show()

if 'snakemake' in locals():
    with PdfPages(snakemake.output[0]) as pdf:
        pdf.savefig(fig)
    


