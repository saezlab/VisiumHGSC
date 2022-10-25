import scanpy as sc
from os.path import join, normpath

input_dir = snakemake.input[0]

### Load RNA modalities
print('INFO: Loading Spatial assays')
spatial = sc.read_h5ad(join(input_dir,'Spatial.h5ad'))
# spatial.layers['counts'] = spatial.layers['logcounts'].copy()
del spatial.layers['logcounts']
print(spatial)

# sort genes alphabetically
spatial = spatial[:,spatial.var.index.sort_values()]


sct = sc.read_h5ad(join(input_dir,'SCT.h5ad'))

#change layer name
del sct.layers['logcounts']

# ensure same order of cells, and alphabetically sorted genes
sct = sct[spatial.obs.index, sct.var.index.sort_values()]

# copy counts from spatial assay to raw
sct.raw = spatial[sct.obs.index].copy()

sct.layers['counts'] = spatial[sct.obs.index, sct.var.index].X.copy()

#change keys to lower case
for obsm in list(sct.obsm.keys()):
    sct.obsm['X_' + obsm.lower()] = sct.obsm[obsm].copy()
    del sct.obsm[obsm]

print(sct)

sct.write_h5ad(snakemake.output[0])
print('INFO: Combined adata written to', snakemake.output[0])
