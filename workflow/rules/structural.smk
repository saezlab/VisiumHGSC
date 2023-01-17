rule cell_density:
    input:
        'data/cell_per_location.csv',
        'results/integrated/sample_metadata.csv'
    output:
        'plots/structural/cell_density.pdf'
    conda:
        '../envs/misty.yaml'
    script:
        '../scripts/structural/plot_cell_density.R'


# merge normalised data from Seurat with
# results from cell2location (provided by collaborators)
rule merge_adata:
    input:
        'results/integrated/{split_type}.h5ad',
        expand('data/cell2location/{sample}_sp.h5ad', sample = config['samples'])
    output:
        merged = 'results/structural/{split_type}/merged.h5ad'
    conda:
        "../envs/scanpy.yaml"
    script:
        '../scripts/structural/merge_adata.py'

#make plots of cell2location output (in proportions)
rule cell_prop:
    input:
        merged = 'results/structural/{split_type}/merged.h5ad',
        dc_CT = 'results/structural/{split_type}/scores_cell_subtype.csv'
    output:
        c2l = 'plots/structural/{split_type}/cell2location_umap.pdf',
        dc = 'plots/structural/{split_type}/dc_umap.pdf',
        spat1 = 'plots/structural/{split_type}/cell2location_spatial1.pdf',
        spat2 = 'plots/structural/{split_type}/cell2location_spatial2.pdf'
    conda:
        "../envs/scanpy.yaml"
    script:
        '../scripts/structural/umaps_CTprop.py'

############################
#   Use prior atlas to get cell type markers using decoupler
#   Atlas GEO accession = GSE165897
############################
rule get_refCT_markers:
    input:
        counts = 'data/sc-data/GSE165897_UMIcounts_HGSOC.tsv',
        meta = 'data/sc-data/GSE165897_cellInfo_HGSOC.tsv'
    output:
        'results/structural/ref_markers_{CT}.csv'
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/structural/get_refCT_markers.py"

# use cell type markers identified in get_refCT_markers
# and infer on our spatial data with decoupler
rule CT_scores:
    input:
        integrated = "results/integrated/{split_type}.h5ad",
        net = 'results/structural/ref_markers_{CT}.csv'
    output:
        'results/structural/{split_type}/scores_{CT}.csv'
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/structural/CT_sig.py"