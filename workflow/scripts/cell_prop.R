library(tidyverse)
library(Seurat)
library(harmony)

if(exists("snakemake")){
  log <- file(snakemake@log[[1]], open="wt", encoding = "UTF-8")
  sink(log, type = c("output", "message"))
}

# define input output files -----------------------------------------------

if(exists("snakemake")){
  input_fp = normalizePath(unlist(snakemake@input[[1]]))
  split_type <- snakemake@wildcards[[1]]
  output_fp <- snakemake@output[[1]]
}else{
  input_fp = normalizePath("../../../data/integrated/plate_int.RDS")
  cell_score_fp = normalizePath("../../../data/integrated/plate_scores_cell_subtype.csv")
  groupby <- c("")
  output_fp <- "test_int.RDS"
}

# load and merge data -----------------------------------------------------

data <- readRDS(input_fp)




# cell proportions from annot and clusters --------------------------------

manual <- data@meta.data %>% filter(annotation_manual != '', Confidence == 'High confidence') %>% add_count(Sample, name = 'n_tot') %>% group_by(Sample, Confidence, PFI, n_tot, annotation_manual) %>% tally(name = 'n_spots') %>% mutate(prop = n_spots/n_tot)

seurat <- data@meta.data %>% filter(Confidence == 'High confidence') %>% add_count(Sample, name = 'n_tot') %>% group_by(Sample, Confidence, PFI, n_tot, seurat_clusters) %>% tally(name = 'n_spots') %>% mutate(prop = n_spots/n_tot)

ggplot(manual, aes(x=PFI, y=prop, fill=PFI)) + geom_boxplot() + 
  labs(title = 'Prop. of manually annotated cell types in HighConf cores') +
  geom_jitter(shape=16, position=position_jitter(width = 0.2, height = 0)) + facet_wrap(~annotation_manual, scales = 'free_y')

ggplot(seurat, aes(x=PFI, y=prop, fill=PFI)) + geom_boxplot() +
  labs(title = 'Prop. of seurat clusters in HighConf cores') +
  geom_jitter(shape=16, position=position_jitter(width = 0.2, height = 0)) + facet_wrap(~seurat_clusters)#, scales = 'free_y')





# cell sig scores ---------------------------------------------------------

cell_sig = read.csv(cell_score_fp)
types <- colnames(cell_sig)[2:ncol(cell_sig)]

data@meta.data <- inner_join(data@meta.data %>% tibble::rownames_to_column(var='cell_id'), cell_sig, by = c('cell_id' = 'X')) %>% tibble::column_to_rownames(var = 'cell_id')

scores <- data@meta.data %>% filter(Confidence == 'High confidence') %>% select(Sample, PFI, all_of(types))

scores <- pivot_longer(scores, cols = 3:ncol(scores), names_to ="cell_subtype", values_to = "score")

# ggplot(scores %>% filter(score > 0), aes(x=PFI, y=score, fill=PFI)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
#   labs(title = 'Decoupler cell type scores in HighConf cores') +
#   #geom_jitter(shape=16, position=position_jitter(width = 0.2, height = 0)) +
#   facet_wrap(~cell_subtype, nrow = 3)


ggplot(scores %>% filter(score > 0) %>% group_by(Sample, PFI, cell_subtype) %>% summarise(meanscore = mean(score)), aes(x=PFI, y=meanscore, fill=PFI)) + geom_boxplot() +
  labs(title = 'Mean Decoupler cell type scores in HighConf cores') +
  geom_jitter(shape=16, position=position_jitter(width = 0.2, height = 0)) + facet_wrap(~cell_subtype)#, scales = 'free_y')






# spatial representations -------------------------------------------------


cell_type_plots <- lapply(types, function(celltp){

  plots <- lapply(names(data@images), function(image){
    p <- SpatialPlot(data, features = celltp, images = image) + theme(plot.title = element_text(hjust = 0.5), legend.position = "right",) + labs(title = image)
    p <-p + scale_fill_gradient2(low = scales::muted("blue"),
                                 mid = "white",
                                 high = scales::muted("red"),
                                 midpoint = 0,
                                 space = "Lab",
                                 na.value = "grey50",
                                 guide = "colourbar",
                                 aesthetics = "fill",
                                 limits = c(min(data@meta.data[celltp]), max(data@meta.data[celltp]))
                                 )
    
    coord <- GetTissueCoordinates(data, image)
    maxX <- max(GetTissueCoordinates(data)$imagerow) + 20
    
    meta <- data@meta.data %>% filter(plate == image) %>% select(Sample, Confidence, PFI) %>% filter(Confidence == 'High confidence')
    coord <- inner_join(coord %>% tibble::rownames_to_column(var='cell_id'), meta %>% tibble::rownames_to_column(var='cell_id'), by = 'cell_id')
    centroids <- coord %>% select(Sample, PFI, imagerow, imagecol) %>%  group_by(Sample, PFI) %>% 
      summarise(meanrow = mean(imagerow), meancol = mean(imagecol), .groups = 'drop_last') %>% 
      mutate(meanrow = abs(maxX-meanrow))
    
    p <- p + geom_label(mapping = aes(y = meanrow, x= meancol, label = Sample, color = PFI), data = centroids, size = 3, inherit.aes = FALSE) + guides(color = FALSE)
    
    if(image != tail(names(data@images), n=1)){
      p <- p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    }
    
    return(p)
    
  })
  p <- patchwork::wrap_plots(plots, ncol = length(plots)) + patchwork::plot_annotation(title = celltp, theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))
  return(p)
  
})
names(cell_type_plots) <- types

lapply(names(cell_type_plots), function(celltp){
  
  filename <- paste('../../../plots/cell_type_scores/' , gsub('\\.','', celltp), '.pdf',sep = '')
  
  # pdf(file=filename, width=15, height=5)
  print(cell_type_plots[celltp])
  # dev.off()
})













