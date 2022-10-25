# %%
import scanpy as sc
import pandas as pd

# %%
if 'snakemake' in locals():
    adata_fp = snakemake.input[0]
    output1_fp = snakemake.output[0]
    output2_fp = snakemake.output[1]
else:
    adata_fp = 'results/structural/plate/merged.h5ad'
    output_fp = 'test.csv'


# %%
adata = sc.read_h5ad(adata_fp)

# %%
coords = adata.obs.filter(['array_row', 'array_col', 'Sample'])
coords[['x','y']] = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'], index = coords.index)

coords.to_csv(output1_fp)

# %%
prop = adata.obsm['q05_cell_abundance_w_sf']

prop = prop.div(prop.sum(axis=1), axis = 0)
prop.columns = adata.uns['mod']['factor_names']
prop.index.name = None

prop.to_csv(output2_fp)