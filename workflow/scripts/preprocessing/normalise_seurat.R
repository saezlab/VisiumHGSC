library(tidyverse)
library(Seurat)

if(exists("snakemake")){
  log <- file(snakemake@log[[1]], open="wt", encoding = "UTF-8")
  sink(log, type = c("output", "message"))
}

# define input output files -----------------------------------------------

if(exists("snakemake")){
  input_fp = normalizePath(snakemake@input$raw)
  output_fp <- snakemake@output[[1]]
}else{
  input_fp = normalizePath("../../../data/normalised/plate/TMA1_norm.RDS")
  output_fp <- "test_norm.h5Seurat"
}


# load an normalise -------------------------------------------------------

data <- readRDS(input_fp) %>% SCTransform(assay = "Spatial")
cat("INFO: normalisation done")
# data <- RunPCA(data, assay = "SCT")
# cat("INFO: PCA done on SCT")

if(exists("snakemake")){
  
  saveRDS(data, file = output_fp)
  
  sink(file = NULL)
}