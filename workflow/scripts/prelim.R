library(tidyverse)
library(Seurat)
library(harmony)

fp <- normalizePath("../../../data/merged_samples.RDS")
# fp <- normalizePath("../../../data/integrated/plate_int.RDS")

data <- readRDS(fp)

# images ------------------------------------------------------------------

SpatialFeaturePlot(data, features = c("nCount_Spatial", "nFeature_Spatial"))

VlnPlot(data, features = c("nCount_SCT", "nFeature_SCT"), group.by = "PFI")








# temp --------------------------------------------------------------------

a <- data@meta.data %>% distinct(Sample, plate, Confidence, PFI) 

a %>% count(plate, PFI) %>% ggplot(aes(x=PFI, y=n, fill=PFI)) +
  geom_bar(stat="identity") + facet_wrap(~plate)



data <- RunUMAP(data, reduction = "harmony", dims = 1:30)

data <- SetIdent(data, value = "Sample")
plots <- lapply(names(data@images), function(image){
  p <- SpatialDimPlot(data, images = image) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + labs(title = image)
  
  return(p)
  
})
patchwork::wrap_plots(plots, ncol = length(plots)) + patchwork::plot_annotation(title = paste("Samples"), theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))





# get cluster -------------------------------------------------------------

plots <- lapply(names(data@images), function(image){
  p <- SpatialDimPlot(data, images = image) + theme(plot.title = element_text(hjust = 0.5)) + labs(title = image)
  
  if(image != tail(names(data@images), n=1)){
    p <- p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  }
  
  return(p)
  
})
patchwork::wrap_plots(plots, ncol = length(plots)) + patchwork::plot_annotation(title = "Pre-existing clusters", theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))

DimPlot(data, reduction = "umap", group.by = c("patient", "Sample", "seurat_clusters", "cell_type"))

data <- SetIdent(data, value = "cell_type")
SpatialDimPlot(data)
