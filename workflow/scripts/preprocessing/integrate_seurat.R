library(tidyverse)
library(Seurat)
library(harmony)

if(exists("snakemake")){
  log <- file(snakemake@log[[1]], open="wt", encoding = "UTF-8")
  sink(log, type = c("output", "message"))
}

# define input output files -----------------------------------------------

if(exists("snakemake")){
  raw_fp = normalizePath(snakemake@input$raw)
  input_fps = normalizePath(unlist(snakemake@input$samples))
  split_type <- snakemake@wildcards[[1]]
  output_fp <- snakemake@output[[1]]
}else{
  raw_fp = normalizePath("../../../data/merged_samples.RDS")
  input_fps = normalizePath(c("../../../data/normalised/plate/TMA1_norm.RDS", "../../../data/normalised/plate/TMA2_norm.RDS", "../../../data/normalised/plate/TMA3_norm.RDS", "../../../data/normalised/plate/TMA4_norm.RDS"))
  split_type <- "plate"
  output_fp <- "test_int.RDS"
}

# load and merge data -----------------------------------------------------

data <- lapply(input_fps, readRDS)
features <- SelectIntegrationFeatures(object.list = data, nfeatures = 3000)

if(split_type == "raw"){
  data <- data[[1]]
}else{
  data <- merge(data[[1]], data[2:length(data)], merge.data = TRUE)  
}


# data@images <- data@images[c("TMA1", "TMA2", "TMA3", "TMA4")]



data <- RunPCA(data, assay = "SCT", features = features)

data <- RunHarmony(data, c("patient", "plate"), assay.use = "SCT", reduction.save = "harmony", max.iter.harmony = 20)

data <- FindNeighbors(data, reduction = "harmony")
print('Neighbors added')
data <- FindClusters(data, algorithm = 4)
print('Clusters added')

data <- RunUMAP(data, reduction = "harmony", dims = 1:30)
print('Umaps added')

if(exists("snakemake")){
  raw <- readRDS(raw_fp)
  data@images <- raw@images
  print('Images added')
  
  saveRDS(data, file = output_fp)
  
  sink(file = NULL)
}