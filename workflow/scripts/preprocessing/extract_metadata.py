# %%
import scanpy as sc
import pandas as pd

# %%
# Define input and output files
if 'snakemake' in locals():
    data_fp = snakemake.input[0]

else:
    data_fp = 'results/integrated/plate.h5ad'

# %%
adata = sc.read_h5ad(data_fp)
metadata = adata.obs.filter(['Sample','Confidence', 'PFI', 'patient'],axis =1).drop_duplicates().reset_index(drop=True)

# %%
if 'snakemake' in locals():
    metadata.to_csv(snakemake.output[0], sep = ',', index = False)


