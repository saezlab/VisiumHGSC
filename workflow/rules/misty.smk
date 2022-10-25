rule get_coords_cellprop:
    input:
        data = 'results/structural/plate/merged.h5ad'
    output:
        coords = 'results/Misty/coordinates.csv',
        deconv = 'results/Misty/cellprop.csv'
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/misty/extract_coordinates.py"


rule get_celltype_views:
    input:
        'results/Misty/coordinates.csv',
        'results/Misty/cellprop.csv'
    output:
        view = 'results/Misty/models/{sample}/celltype_view.rds'
    conda:
        "../envs/misty.yaml"
    script:
        "../scripts/misty/make_views.R"

rule run_views:
    input:
        view = 'results/Misty/models/{sample}/{view_type}_view.rds'
    output: 
        directory('results/Misty/models/{sample}/{view_type}_misty_model')
    params:
        seed = config['misty'].get("random_seed", 42),
        bypass_intra = lambda wildcards: config['misty'][wildcards.view_type].get('bypass_intra', False)
    conda:
        "../envs/misty.yaml"
    threads: 6
    script:
        "../scripts/misty/run.R"