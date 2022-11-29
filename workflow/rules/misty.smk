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
        view = 'results/Misty/celltype/views/{sample}_view.rds',
        paraview = 'results/Misty/celltype/views/{sample}_paraview.csv'
    conda:
        "../envs/misty.yaml"
    script:
        "../scripts/misty/make_views.R"

rule get_combine_paraviews:
    input:
        expand('results/Misty/celltype/views/{sample}_paraview.csv', sample = config['samples'])
    output:
        'results/Misty/{view_type}/paraviews.csv'
    shell:
        "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} >> {output}"

rule run_views:
    input:
        view = 'results/Misty/{view_type}/views/{sample}_view.rds'
    output: 
        directory('results/Misty/{view_type}/models/{sample}')
    params:
        seed = config['misty'].get("random_seed", 42),
        bypass_intra = lambda wildcards: config['misty'][wildcards.view_type].get('bypass_intra', False)
    conda:
        "../envs/misty.yaml"
    threads: 6
    script:
        "../scripts/misty/run.R"