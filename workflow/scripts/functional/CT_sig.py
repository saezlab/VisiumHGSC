# %%
import scanpy as sc
import pandas as pd
import decoupler as dc

# %%
adata = sc.read_h5ad(snakemake.input[0])
adata.var_names = adata.var.features

# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
adata

# %%
net = pd.read_csv(snakemake.input[1])

# %%
dc.decouple(adata, net = net, source = 'cell_type', weight = 'logFCs', methods = 'consensus', use_raw = False)

# %%
acts = dc.get_acts(adata, obsm_key='consensus_estimate')


# %%
pd.DataFrame(acts.X, index=acts.obs_names, columns=acts.var_names).to_csv(snakemake.output[0], sep=',')

# %%
# sc.pl.umap(acts, color=list(acts.var_names))


