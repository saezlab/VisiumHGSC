# %%
import scanpy as sc
import pandas as pd
import decoupler as dc

# %%

counts_fn = snakemake.input[0]
meta_fn = snakemake.input[1]

level = snakemake.wildcards[0]

# %%
adata = sc.read_text(counts_fn, delimiter='\t').T
meta = pd.read_table(meta_fn, delimiter='\t', index_col=0)
adata.obs = meta
adata.layers['counts'] = adata.X

# %%
keep = adata.obs['treatment_phase'] == 'treatment-naive'
adata = adata[keep,:]

# %%
padata = dc.get_pseudobulk(adata, sample_col='patient_id', groups_col=level, layer='counts', min_prop=0.2, min_smpls=3)
sc.pp.normalize_total(padata, target_sum=1e4)
sc.pp.log1p(padata)
padata

# %%
cell_types = padata.obs[level].unique()
net = []
cutoff = -10
for cell_type in cell_types:
    logFCs, pvals = dc.get_contrast(padata, None, condition_col=level, condition= cell_type, reference='rest', method='t-test')
    deg = dc.format_contrast_results(logFCs, pvals)
    sign = dc.get_top_targets(logFCs, pvals, cell_type + '.vs.rest', sign_thr=0.05, lFCs_thr=0.5)
    sign = sign.reset_index().rename({'index':'target'}, axis = 'columns')
    sign['cell_type'] = cell_type
    keep = sign['logFCs'] > cutoff
    sign = sign[keep]
    net.append(sign)
    print('cell type: ', cell_type, ' done')
net = pd.concat(net, ignore_index=True).drop(['target'], axis = 1).rename({'name':'target'}, axis = 1)

print(net.head())

#remove mitochondrial genes
not_mito = [not gene.startswith('MT-') for gene in net['target']]
net = net[not_mito]

net.to_csv(snakemake.output[0], index=False)

