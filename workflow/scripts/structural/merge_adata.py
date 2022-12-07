# %%
import scanpy as sc
import sys
import pandas as pd
import decoupler as dc

# %%
if 'snakemake' in locals():
    adata_s_fp = snakemake.input[0]
    adata_fp = snakemake.input[1]

    if len(snakemake.input) > 2:
        adata_fps = snakemake.input[1:]
else:
    adata_fp = '../../../data/sp.h5ad'
    adata_s_fp = '../../../results/integrated/plate.h5ad'


# %%

if 'adata_fps' in locals():
    adata = sc.concat([sc.read_h5ad(file) for file in adata_fps], join="outer")
    temp = sc.read_h5ad(adata_fp)
    adata.uns = {'mod': {'factor_names': temp.uns['mod']['factor_names'].copy()}}
    adata.uns['spatial'] = temp.uns['spatial']
    del temp
else:
    adata = sc.read_h5ad(adata_fp)
print(adata)

# %%
adata_seurat = sc.read_h5ad(adata_s_fp)

#change index naming
adata_seurat.obs.index = ['TMA' + cell[1] + '_' + cell[0] + '-1' for cell in adata_seurat.obs.index.str.split('-')]

# %%
print(adata_seurat)

# %%
#sort renamed indexes in seurat object
idx = adata_seurat.obs.index.sort_values()
adata_seurat = adata_seurat[idx]

#sort first anndata object in the same order as the second
adata = adata[adata_seurat.obs.index].copy()

# %%
# transfer obsm to anndata from Seurat
for id, value in adata.obsm.items():
    if not id in list(adata_seurat.obsm.keys()):
        adata_seurat.obsm[id] = value

# transfer uns to anndata from Seurat
for id, value in adata.uns.items():
    if not id in list(adata_seurat.uns.keys()):
        adata_seurat.uns[id] = value

obs_meta = adata.obs.filter(items=['in_tissue', 'array_row', 'array_col','_indices', '_scvi_batch', '_scvi_labels'])

# %%
#transfer spatial coordinates to anndata from Seurat
adata_seurat.obs[obs_meta.columns] = obs_meta
adata_seurat.obs

adata_seurat.var = adata_seurat.var.rename({'_index':'_indices'},axis = 1).copy()

# %%
print(adata_seurat)

# %%
temp = adata_seurat.raw.to_adata()
temp.var = temp.var.rename({'_index':'_indices'},axis = 1).copy()
adata_seurat.raw = temp

# %%
if 'snakemake' in locals():
    adata_seurat.write_h5ad(snakemake.output[0])
else:
    adata_seurat.write_h5ad('test.h5ad')
