############################
#   Do different types of normalisation
#   Integrate w. harmony
############################

checkpoint split_seurat:
    input:
        raw = "data/merged_samples.RDS"
    output:
        directory("results/split/{split_type}")
    log:
        "logs/split/{split_type}.log"
    conda:
        '../envs/preprocessingR.yaml'
    script:
        "../scripts/preprocessing/split_seurat.R"

rule normalise_seurat:
    input:
        raw = "results/split/{split_type}/{i}_raw.RDS"
    output:
        normalised =  "results/normalised/{split_type}/{i}_norm.RDS"
    log:
        "logs/normalised/{split_type}_{i}.log"
    conda:
        '../envs/preprocessingR.yaml'
    script:
        "../scripts/preprocessing/normalise_seurat.R"

def aggregate_norm_data(wildcards):
    checkpoint_output = checkpoints.split_seurat.get(**wildcards).output[0]
    return expand("results/normalised/{split_type}/{i}_norm.RDS", split_type=wildcards.split_type, i=glob_wildcards(os.path.join(checkpoint_output, "{i}_raw.RDS")).i)

rule integrate_seurat:
    input: 
        samples = aggregate_norm_data,
        raw = "data/merged_samples.RDS"
    output:
        integrated = "results/integrated/{split_type}_int.RDS"
    log:
        "logs/integrated/{split_type}_int.log"
    conda:
        '../envs/preprocessingR.yaml'
    script:
        "../scripts/preprocessing/integrate_seurat.R"

rule integration_plots:
    input:
        integrated = "results/integrated/{split_type}_int.RDS"
    output:
        umaps = "plots/integration/{split_type}_umaps.pdf",
        clusters = "plots/integration/{split_type}_clusters.pdf"
    log:
        "logs/plots/{split_type}_int.log"
    conda:
        '../envs/preprocessingR.yaml'
    script:
        "../scripts/preprocessing/integration_plots.R"

############################
#   Convert from Seurat to Anndata
############################

rule seurat_to_folder:
    input:
        integrated = "results/integrated/{split_type}_int.RDS"
    output:
        h5ad = temp(directory('results/convert/{split_type}'))
    conda:
        "../envs/preprocessingR.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"

rule make_anndata:
    input:
        rules.seurat_to_folder.output
    output:
        ad = 'results/integrated/{split_type}.h5ad'
    conda:
        "../envs/preprocessingR.yml"
    script:
        "../scripts/preprocessing/ST_to_adata.py"