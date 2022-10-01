rule merge_adata:
    input:
        adata_w_c2l = 'data/sp.h5ad',
        integrated = "data/integrated/plate_int.h5ad"
    output:
        merged = 'data/merged.h5ad'
    script:
        '../scripts/structural/merge_adata.py'

rule cell_prop:
    input:
        merged = 'data/merged.h5ad',
        dc_CT = 'data/integrated/plate_scores_cell_subtype.csv'
    output:
        c2l = 'plots/CTprop/cell2location_umap.pdf',
        dc = 'plots/CTprop/dc_umap.pdf',
        spat1 = 'plots/CTprop/cell2location_spatial1.pdf',
        spat2 = 'plots/CTprop/cell2location_spatial2.pdf'
    script:
        '../scripts/structural/umaps_CTprop.py'