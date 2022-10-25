############################
#   Functional anntation of cells with TFs and pathways
############################

def activities_inputs(wildcards):
    files = {'data': 'results/structural/{wildcards.split_type}/merged.h5ad'.format(wildcards=wildcards)}
    if wildcards.network == 'cytokines':
        files['net'] = 'resources/model_matrix_zscore_clean_full_assocs.csv'
    return files


rule get_activities:
    input:
        unpack(activities_inputs)
    output:
        act = 'results/functional/{split_type}/activities_{network}.csv' #network being either 'pathways', TFs or 'cytokines'
    params:
        lambda w: config['functional'][w.network], # set parameters in config file
        organism = config.get('organism', 'human')
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/functional/compute_activities.py"

rule plot_activities:
    input:
        data = 'results/structural/{split_type}/merged.h5ad',
        activities = 'results/functional/{split_type}/activities_{network}.csv'
    output:
        'plots/functional/{split_type}/spatial_activities_{network}.pdf'
    params:
        lambda w: config['functional'][w.network]
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/functional/plot_spatial_activities.py"