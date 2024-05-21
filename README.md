# Spatial analysis of ovarian cancer tissue micro-arrays

## Description
This repository contains the code necessary to reproduce the spatial analysis using the Misty package from the paper "Opening the black box: spatial transcriptomics and the relevance of AI-detected prognostic regions in high grade serous carcinoma". It is wrapped in a snakemake pipeline that takes care of installing the required dependencies and running the analysis.

## Installation
You need a working `snakemake` installation to run the pipeline. The pipeline has only been tested using snakemake version up to 7.24.0.

You can reproduce the analysis locally by running the following commands:

```bash
snakemake --use-conda -c 1
```

This will create a conda environment with the required dependencies and run the analysis. We also provide an example snakemake profile in `config/slurm` that can be used to run the pipeline on a SLURM cluster.

```bash
snakemake --profile config/slurm
```

By default, the pipeline produces three pdf files that contain the plots shown in Figure 7 and in supplementary figure 6 of the paper. The plots are saved in the `plots` directory.
- Figures 7A-B correspond to pages 10 and 11 of `plots/Misty/celltype_misty.pdf`
- Figure 7C is page 2 of `plots/Misty/celltype_interaction_correlations.pdf`
- Figure 7D is page 1 of `plots/spatial_locations.pdf`
- Figure S6A and S6C correspond to pages 12 and 13 of `plots/Misty/celltype_misty.pdf`
- Figure S6B is page 3 of `plots/Misty/celltype_interaction_correlations.pdf`

## References
Anna Ray Laury, Shuyu Zheng, et al., **Opening the black box: spatial transcriptomics and the relevance of AI-detected prognostic regions in high grade serous carcinoma**, Modern Pathology, 2024, 100508, ISSN 0893-3952, https://doi.org/10.1016/j.modpat.2024.100508.
