library(tidyverse)
library(Seurat)

if(exists("snakemake")){
  log <- file(snakemake@log[[1]], open="wt", encoding = "UTF-8")
  sink(log, type = c("output"))
}


# define input output files -----------------------------------------------

if(exists("snakemake")){
  raw_fp = normalizePath(snakemake@input$raw)
  split_type <- snakemake@wildcards[[1]]
  output_dir <- normalizePath(snakemake@output[[1]])
}else{
  raw_fp = normalizePath("../../../data/merged_samples.RDS")
  split_type <- "plate"
  output_dir <- "../../../data/split/plate"
}


# load data ---------------------------------------------------------------

data <- readRDS(raw_fp)
if(split_type == "raw"){
  data@meta.data$to.split <- "raw"
}else{
  data@meta.data$to.split <- data@meta.data[[split_type]]
}

#clean object
data@active.assay <- "Spatial"
data@reductions <- list()
data@commands <- list()
data <- DietSeurat(data, counts = TRUE, assays = "Spatial")
#reset ident
data <- SetIdent(data, value = "patient")

#clean metadata
data@meta.data <- data@meta.data %>% select(-contains("SCT")) %>% rename(seurat_clusters_old = seurat_clusters)

#get total number of plates
plates <- unique(data@meta.data$plate)

#define which way to split the seurat object
groups <- unique(data@meta.data$to.split)

cat("INFO: Will create raw Seurat objects for the following levels:\n", groups, "\n")


# do the split ------------------------------------------------------------

datas <- lapply(groups, function(group){
  temp_data <- subset(data, subset = to.split == group)

  # #remove unneeded images
  # images.remove <- setdiff(plates, unique(temp_data@meta.data$plate))
  # if(length(images.remove) != 0) temp_data@images[images.remove] <- NULL
  
  #remove unneeded split column
  temp_data@meta.data <- temp_data@meta.data %>% select(-to.split)
  
  return(temp_data)
})
names(datas) <- groups


# output file handling ----------------------------------------------------

if(exists("snakemake")){
  
  if(dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
  dir.create(output_dir)
  output_dir <- normalizePath(output_dir)
  
  for (group in groups) {
    output_fp <- file.path(output_dir, paste(group, "raw.RDS", sep = "_"))
    saveRDS(datas[[group]], output_fp)
    cat("INFO: saved Seurat object to:", output_fp, "\n")
  }
  
  if(Sys.info()['sysname'] == "Darwin" & Sys.info()['effective_user'] == "demian"){
    line <- paste("say Robin, the", split_type, "has been split", sep = " ")
    system(line)
  }
  
  sink(file = NULL)
}
