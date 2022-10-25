# Define input and output files
if 'snakemake' in locals():
    data_fp = snakemake.input['data']
    act_fp = snakemake.input['activities']
    output_fp = snakemake.output[0]

    # parameters for decoupler run
    network = snakemake.wildcards.network
    conf = snakemake.params[0]

print(locals())