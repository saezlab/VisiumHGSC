############################
#   Functional anntation of cells with TFs and pathways
############################

rule get_activities:
    input:
        data = 'results/integrated/{split_type}.h5ad'
    output:
        act = 'results/functional/{split_type}/activities_{network}.csv' #network being either 'pathways' or 'cytokines'
    params:
        normalisation = config['functional'].get("normalisation", 'log1p'), #or SCT
        top_genes = config['functional'].get("top_targets", 300),
        method = config['functional'].get("method", 'mlm')
    script:
        "../scripts/functional/compute_activities.py"