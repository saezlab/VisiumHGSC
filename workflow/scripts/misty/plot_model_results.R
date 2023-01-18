# load libraries ----------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(mistyR)
# library(factoextra)

# set input output files --------------------------------------------------

cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  
  view <- snakemake@wildcards$view_type
  
  plot_params <- snakemake@params[[1]]
  sig.cutoff <- snakemake@params$sig_cutoff
  
  helpers_fp <- snakemake@input[[1]]
  metadata_fp <- snakemake@input[[2]]
  result_folders <- unlist(snakemake@input[3:length(snakemake@input)])
  
}else{
  
  # view <- 'celltype'
  view <- 'celltype'
  
  plot_params <- list(trim = 1, cutoff = 1)
  sig.cutoff <- snakemake@params$sig_cutoff

  helpers_fp <- 'workflow/scripts/helpers/misty.R'
  result_folders <- paste('results/Misty', view, 'models', sep = .Platform$file.sep) %>% list.files(full.names = TRUE)
  
  #files for testing in Rstudio
  metadata_fp <- 'results/integrated/sample_metadata.csv'
  
}

source(helpers_fp)

# load data ---------------------------------------------------------------

result_folders <-  data.frame(path = result_folders) %>% mutate(sample = basename(path))

metadata <- read.csv(metadata_fp) %>% rename(sample =Sample)
metadata$Confidence <- factor(metadata$Confidence, levels = c('High confidence', 'Low confidence', 'Benign'))
metadata$PFI <- factor(metadata$PFI, levels = c('Short', 'Long', 'Benign'))

metadata <- left_join(metadata, result_folders, by = 'sample')

results <- metadata %>% pull(path) %>% collect_results() %>% reformat_samples()


if(view == 'celltype'){
  intra_name <- 'intra'
  cleaning <- FALSE
}

if(view == 'pathwaysCT'){
  intra_name <- 'intra_act'
  cleaning <- TRUE
}

if(view == 'CTpathways'){
  intra_name <- 'intra_act'
  cleaning <- TRUE
}


# pca plots ---------------------------------------------------------------

if(exists("snakemake")) pdf(snakemake@output[[1]])

imp.signature <- extract_signature(results, type = "importance")
imp.signature.pca <- prcomp(imp.signature %>% select(-sample))
# imp.signature.umap <- uwot::umap(imp.signature %>% select(-sample))

ggplot(
  left_join(bind_cols(sample = metadata %>%  pull(sample), imp.signature.pca$x), metadata, by = 'sample'),
  aes(x = PC1, y = PC2)
) +
  geom_point(aes(color = as.factor(PFI)), size = 1) +
  labs(color = "PFI") +
  theme_classic()


ggplot(
  left_join(bind_cols(sample = metadata %>%  pull(sample), imp.signature.pca$x), metadata, by = 'sample'),
  aes(x = PC1, y = PC2)
) +
  geom_point(aes(color = as.factor(Confidence)), size = 1) +
  labs(color = "Confidence") +
  theme_classic()

ggplot(
  left_join(bind_cols(sample = metadata %>%  pull(sample), imp.signature.pca$x), metadata, by = 'sample'),
  aes(x = PC1, y = PC2)
) +
  geom_point(aes(color = as.factor(patient)), size = 1) +
  labs(color = "patient") +
  theme_classic()

# fviz_pca_var(imp.signature.pca,
#              col.var = "cos2", select.var = list(cos2 = 15), repel = TRUE,
#              gradient.cols = c("#666666", "#377EB8", "#E41A1C"), col.circle = NA
# ) + theme_classic()

# general plots -----------------------------------------------------------

if(view != 'functional'){
  results %>% plot_improvement_stats()
  
  results %>% plot_improvement_stats("intra.R2")
  
  results %>% plot_view_contributions(trim = 1)
  
  results %>% plot_interaction_heatmap(intra_name, cutoff = 0.5, clean = cleaning)
  
  results %>% plot_interaction_heatmap('para', cutoff = 0.5, trim = 0.5)
}



# BG vs. HC ---------------------------------------------------------------

contrast.interactions <- get_differential_interactions(metadata, 'Confidence', levels(metadata %>% pull(matches('Confidence')))[1:2])

views <- contrast.interactions %>% dplyr::pull(.data$view) %>% unique() %>% sort()

views %>% purrr::walk(function(current.view){
  
  long.data <- contrast.interactions %>% dplyr::filter(.data$view == current.view) %>%
    dplyr::mutate(sig = -log10(.data$p.adj) * sign(.data$t.value), is.sig = .data$p.adj < 0.05)
  
  inter.plot <- long.data %>% ggplot2::ggplot(ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$sig)) + ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggplot2::coord_equal() + 
    ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "#8DA0CB", midpoint = 0) + ggplot2::labs(fill='-log10(p)') +
    geom_tile(data = long.data %>% dplyr::filter(.data$is.sig), aes(fill = .data$sig, color = is.sig), size = 0.5) + 
    scale_color_manual(guide = FALSE, values = c(`TRUE` = "black")) +
    ggplot2::ggtitle(paste('Condition specific interactions in ', current.view, '\nHC vs BG', sep = ''))
  
  print(inter.plot)
  
})


# long vs short -----------------------------------------------------------

contrast.interactions <- get_differential_interactions(metadata %>% filter(Confidence == 'High confidence'), 'PFI', levels(metadata %>% pull(matches('PFI')))[1:2])

views <- contrast.interactions %>% dplyr::pull(.data$view) %>% unique() %>% sort()

views %>% purrr::walk(function(current.view){
  
  long.data <- contrast.interactions %>% dplyr::filter(.data$view == current.view) %>%
    dplyr::mutate(sig = -log10(.data$p.adj) * sign(.data$t.value), is.sig = .data$p.adj < 0.05)
  
  inter.plot <- long.data %>% arrange(-p.adj) %>% ggplot2::ggplot(ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$sig)) + ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggplot2::coord_equal() + 
    ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "#8DA0CB", midpoint = 0) + ggplot2::labs(fill='-log10(p)') +
    geom_tile(aes(x = .data$Predictor, y = .data$Target, fill = .data$sig, color = is.sig), size = 0.5) + 
    scale_color_manual(guide = FALSE, values = c(`TRUE` = "black", `FALSE` = "#00000000")) +
    ggplot2::ggtitle(paste('Condition specific interactions in ', current.view, '\nShort vs Long (in HC)', sep = ''))
  
  print(inter.plot)
  
})

# long vs short -----------------------------------------------------------

contrast.interactions <- get_differential_interactions(metadata %>% filter(Confidence == 'Low confidence'), 'PFI', levels(metadata %>% pull(matches('PFI')))[1:2])

views <- contrast.interactions %>% dplyr::pull(.data$view) %>% unique() %>% sort()

views %>% purrr::walk(function(current.view){
  
  long.data <- contrast.interactions %>% dplyr::filter(.data$view == current.view) %>%
    dplyr::mutate(sig = -log10(.data$p.adj) * sign(.data$t.value), is.sig = .data$p.adj < 0.05)
  
  inter.plot <- long.data %>% arrange(-p.adj) %>% ggplot2::ggplot(ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$sig)) + ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggplot2::coord_equal() + 
    ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "#8DA0CB", midpoint = 0) + ggplot2::labs(fill='-log10(p)') +
    geom_tile(aes(x = .data$Predictor, y = .data$Target, fill = .data$sig, color = is.sig), size = 0.5) + 
    scale_color_manual(guide = FALSE, values = c(`TRUE` = "black", `FALSE` = "#00000000")) +
    ggplot2::ggtitle(paste('Condition specific interactions in ', current.view, '\nShort vs Long (in BG)', sep = ''))
  
  print(inter.plot)
  
})

if(exists("snakemake")) dev.off()

