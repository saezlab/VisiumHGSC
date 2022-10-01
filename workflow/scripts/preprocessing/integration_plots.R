library(tidyverse)
library(Seurat)

if(exists("snakemake")){
  log <- file(snakemake@log[[1]], open="wt", encoding = "UTF-8")
  sink(log, type = c("output", "message"))
}

# define input output files -----------------------------------------------

if(exists("snakemake")){
  input_fp = normalizePath(snakemake@input$integrated)
  split_type <- snakemake@wildcards[[1]]
  output_umaps <- snakemake@output$umaps
  output_clusters <- snakemake@output$clusters
}else{
  raw_fp = normalizePath("../../../data/merged_samples.RDS")
  input_fp = normalizePath("../../../data/integrated/plate_int.RDS")
  split_type <- "plate"
  output_umaps <- "test_umaps.pdf"
  output_clusters <- "test_clusters.pdf"
}

data <- readRDS(input_fp)

# qc umaps ----------------------------------------------------------------

pdf(file=output_umaps, width=20, height=10)
to.plot <- c("plate","patient", "PFI", "annotation_manual", "PFI_confidence", "seurat_clusters")

DimPlot(data, ncol = 3, reduction = "umap", group.by = to.plot) + patchwork::plot_annotation(title = paste("Split and normalised by", split_type), theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))

dev.off()

pdf(file=output_clusters, width=15, height=5)
plots <- lapply(names(data@images), function(image){
  p <- SpatialDimPlot(data, images = image) + theme(plot.title = element_text(hjust = 0.5)) + labs(title = image)
  
  if(image != tail(names(data@images), n=1)){
    p <- p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  }
  
  return(p)
  
})
patchwork::wrap_plots(plots, ncol = length(plots)) + patchwork::plot_annotation(title = paste("Seurat clusters (normalised by ", split_type,")",sep=""), theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))
dev.off()

# DimPlot(data, group.by = "nCount_Spatial") #+ theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + labs(title = "TMA1")
