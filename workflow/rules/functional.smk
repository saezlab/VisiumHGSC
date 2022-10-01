
############################
#   Use prior atlas to get cell type markers using decoupler
#   Atlas GEO accession = GSE165897
############################
rule get_refCT_markers:
    input:
        counts = 'data/sc-data/GSE165897_UMIcounts_HGSOC.tsv',
        meta = 'data/sc-data/GSE165897_cellInfo_HGSOC.tsv'
    output:
        'results/sc-data/ref_markers_{CT}.csv'
    script:
        "../scripts/functional/get_refCT_markers.py"

rule CT_scores:
    input:
        integrated = "results/integrated/{split_type}_int.h5ad",
        net = 'results/sc-data/ref_markers_{CT}.csv'
    output:
        'data/integrated/{split_type}_scores_{CT}.csv'
    script:
        "../scripts/functional/CT_sig.py"



############################
#   Functional anntation of cells with TFs and pathways
############################

# rule plot_pathways:
#     input:
#         "data/merged.h5ad"
#     output:
#         'plots/functional/{split_type}_int_pathways.pdf'
#     params:
#         normalisation = config['functional'].get("normalisation", 'log1p'), #or SCT
#         top_genes = config['functional'].get("pathways_top_gene", 300)
#     resources:
#         mem_mb=40000
#     script:
#         '../scripts/functional/plots_pathways.py'

# rule get_PDactivities:
#     input:
#         data = "data/integrated/{split_type}_int.h5ad"
#     output:
#         act = 'data/working/ST/functional/{tissue}_activities_{network}.csv' #network being either 'pathways' or 'TFs'
#     params:
#         normalisation = config['functional'].get("normalisation", 'log1p'), #or SCT
#         top_genes = config['functional'].get("pathways_top_gene", 300),
#         TF_conf = config['functional'].get("TF_confidence", 'ABC'),
#         method = config['functional'].get("pathways_method", 'mlm')
#     script:
#         "../scripts/functional/compute_activities.py"