# %%
import scanpy as sc
import pandas as pd
import decoupler as dc
import logging

# %%
if 'snakemake' in locals():
    #input and output files
    data_fp = snakemake.input['data']
    net_fp = snakemake.input.get('net', '')
    output_fp = snakemake.output[0]

    # parameters for decoupler run
    network = snakemake.wildcards.network
    conf = snakemake.params[0]
    organism = snakemake.params[1]

else:
    #input and output files
    data_fp = '../../../results/integrated/plate.h5ad'
    net_fp = '../../../resources/model_matrix_zscore_clean_full_assocs.csv'

    # parameters for decoupler run
    network = 'pathways'
    conf = {'normalisation': 'log1p', 'top_targets': 300, 'method': 'mlm'}
    organism = 'human'

# %% [markdown]
# Load the data

# %%
adata = sc.read_h5ad(data_fp)
print(adata)

# %% [markdown]
# Check normalisation method that is defined in the config file. It is assumed (for this pipeline) that adata.X contains SCT normalised counts.

# %%
if conf.get('normalisation') == 'log1p':
    print('INFO: using log1p normalised counts')
    adata.layers['SCT'] = adata.X
    adata.X = adata.layers['counts'].copy()
    del adata.layers['counts']

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)


elif conf.get('normalisation') == 'SCT':
    print('INFO: using SCT normalised counts')

elif conf.get('normalisation') is None:
    raise ValueError('The config file is missing a normalisation method for {0}. Set it to either \'log1p\' or \'SCT\'.'.format(network))
else:
    raise ValueError('The normalisation method \'{0}\' is not implemented. Set it to either \'log1p\' or \'SCT\'.'.format(conf.get('normalisation')))

# %% [markdown]
# Setup the regulons passed later on to decoupler, for: progeny pathways, dorothea TFs and Procyto cytokines

# %%
if network == 'pathways':
    model = dc.get_progeny(organism=organism, top= conf.get('top_targets'))

elif network == 'TFs':
    if conf.get('levels') is None:
        logging.warning('No confidence levels were provided for Dorothea regulons. Using high-quality levels ABC.')
    model = dc.get_dorothea(organism=organism, levels=[c for c in conf.get('levels', 'ABC')])

elif network == 'cytokines':
    if net_fp == '':
        raise ValueError('No file was provided for the cytokine regulons!')
    else:
        net = pd.read_csv(net_fp, sep=',', index_col=0)
        if organism == 'mouse': net['target'] = net['target'].str.capitalize()
        
        model = []
        for src in net["source"].unique():
            temp = net[(net["source"] == src)].sort_values(by = "adj.p", axis = 0)

            #select top targets
            if temp.shape[0] < conf.get('top_targets', 100):
                model.append(temp)
            else:
                model.append(temp.iloc[0:conf.get('top_targets', 100),:])

        model = pd.concat(model, axis = 0)
else:
    raise ValueError('The \'network\' wildcard can only take on \'pathways\', \'TFs\' or \'cytokines\' as value, to run either Progeny pathways, Dorothea regulons or cytokine-responsive genes respectively.')

# %% [markdown]
# Run decoupler with the provided method.

# %%
if conf.get('method') is None:
    raise ValueError('The decoupler activity inference method is not defined in the config file!')

dc.decouple(mat=adata, net=model, source='source', target='target', weight='weight', methods = conf.get('method', 'mlm'),  verbose=True, use_raw=False)
print(adata)

# %% [markdown]
# Extract computed activities

# %%
acts = dc.get_acts(adata, obsm_key='{0}_estimate'.format(conf.get('method')))
print(acts)
acts = pd.DataFrame(acts.X, index=acts.obs.index, columns=acts.var.index)
print(acts.head())

# %% [markdown]
# Save to file

# %%
if 'snakemake' in locals():
    acts.to_csv(output_fp)
