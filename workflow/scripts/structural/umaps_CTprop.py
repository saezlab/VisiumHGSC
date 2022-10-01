# %%
import scanpy as sc
import os
import pandas as pd
import decoupler as dc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# %%
if 'snakemake' in locals():
    adata_fp = snakemake.input[0]
    CT_sig_fp = snakemake.input[1]
    c2l = snakemake.output[0]
    dcp = snakemake.output[1]
    spatial1_out = snakemake.output[2]
    spatial2_out = snakemake.output[3]

else:
    adata_fp = '../../../data/merged.h5ad'
    CT_sig_fp = '../../../data/integrated/plate_scores_cell_subtype.csv'

# %%
adata = sc.read_h5ad(adata_fp)

# %%
# adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
prop = adata.obsm['q05_cell_abundance_w_sf']
adata.obs[adata.uns['mod']['factor_names']] = prop.div(prop.sum(axis=1), axis = 0)

# %%
fig, axs = plt.subplots(5, 5, figsize=(23, 15))
axs = axs.flatten()

for ax, CT in zip(axs, sorted(adata.uns['mod']['factor_names'])):

    sc.pl.umap(adata, color = CT, ax = ax, show = False)

for ax in axs[len(adata.uns['mod']['factor_names']):]:
    ax.set_axis_off()

plt.suptitle('Cell 2 location proportions', fontsize = 15)
plt.tight_layout()

if 'snakemake' in locals():
    plt.savefig(c2l, dpi = 300)

# %%
df = pd.read_csv(CT_sig_fp, index_col=0)
df.index = ['TMA' + cell[1] + '_' + cell[0] + '-1' for cell in df.index.str.split('-')]
df = df.loc[adata.obs.index]

# %%
adata.obs['DC_' + df.columns] = df

# %%
fig, axs = plt.subplots(5, 5, figsize=(23, 15))
axs = axs.flatten()

for ax, CT in zip(axs, sorted(df.columns)):

    sc.pl.umap(adata, color = 'DC_' + CT, ax = ax, show = False)
    ax.set_title(CT)

for ax in axs[len(df.columns):]:
    ax.set_axis_off()

plt.suptitle('decoupler activities', fontsize = 15)
plt.tight_layout()
if 'snakemake' in locals():
    plt.savefig(dcp, dpi = 300)

# %%

cell_types = sorted(list(adata.uns['mod']['factor_names']))
plates = adata.obs.filter(['plate'],axis=1).drop_duplicates().sort_values('plate')['plate']
lims = pd.DataFrame({'max':adata.obs.filter(cell_types , axis =1).max()})

with PdfPages(spatial1_out) as output_pdf:
    for plate in plates:

        fig, axs = plt.subplots(5, 5, figsize=(20, 20))
        axs = axs.flatten()

        for i, ct in enumerate(cell_types):

            sc.pl.spatial(
                    adata[adata.obs.plate == plate, :],
                    img_key='lowres',
                    library_id= 'HC-' + plate,
                    color=ct,
                    size=1.5,
                    legend_loc=None,
                    show=False,
                    vmin = 0,
                    vmax = lims.loc[ct,'max'],
                    # vmin = (lims.loc[pathway, 'llim']*1.1),
                    # vmax = (lims.loc[pathway, 'ulim']*1.1),
                    # color_map = 'BrBG',
                    # vcenter = 0,
                    ax=axs[i]
                )
            axs[i].set_title(ct)
            axs[i].set_facecolor('#D9D9D9')
            axs[i].set_ylabel('')
            axs[i].set_xlabel('')

        for ax in axs[len(cell_types):]:
            ax.set_axis_off()

        plt.suptitle(plate, fontsize = 15)
        plt.tight_layout()
        output_pdf.savefig(fig)
        plt.close()

# %%

cell_types = sorted(list(adata.uns['mod']['factor_names']))
plates = adata.obs.filter(['plate'],axis=1).drop_duplicates().sort_values('plate')['plate']
lims = pd.DataFrame({'max':adata.obs.filter(cell_types , axis =1).max()})

with PdfPages(spatial2_out) as output_pdf:
    
    for ct in cell_types:

        fig, axs = plt.subplots(1, 4, figsize=(20, 5))
        axs = axs.flatten()

        for i, plate in enumerate(plates):

            sc.pl.spatial(
                    adata[adata.obs.plate == plate, :],
                    img_key='lowres',
                    library_id= 'HC-' + plate,
                    color=ct,
                    size=1.5,
                    legend_loc=None,
                    show=False,
                    vmin = 0,
                    vmax = lims.loc[ct,'max'],
                    # vmin = (lims.loc[pathway, 'llim']*1.1),
                    # vmax = (lims.loc[pathway, 'ulim']*1.1),
                    # color_map = 'BrBG',
                    # vcenter = 0,
                    ax=axs[i]
                )
            axs[i].set_title(plate)
            axs[i].set_facecolor('#D9D9D9')
            axs[i].set_ylabel('')
            axs[i].set_xlabel('')

        # for ax in axs[len(cell_types):]:
        #     ax.set_axis_off()

        plt.suptitle(ct, fontsize = 15)
        plt.tight_layout()
        output_pdf.savefig(fig)
        plt.close()