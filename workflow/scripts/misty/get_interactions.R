# load libraries ----------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(mistyR)


# set input output files --------------------------------------------------

cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  
  view <- snakemake@wildcards$view_type
  
  plot_params <- snakemake@params[[1]]
  
  metadata_fp <- snakemake@input[[1]]
  result_folders <- unlist(snakemake@input[2:length(snakemake@input)])
  
}else{
  
  view <- 'celltype'
  
  plot_params <- list(trim = 1, cutoff = 1)
  
  
  result_folders <- paste('results/Misty', view, 'models', sep = .Platform$file.sep) %>% list.files(full.names = TRUE)
  
  #files for testing in Rstudio
  metadata_fp <- 'results/integrated/sample_metadata.csv'
  
  source('workflow/scripts/helpers/misty.R') 
}

# load data ---------------------------------------------------------------

if(view == 'celltype'){
  intra_name <- 'intra'
  cleaning <- FALSE
  
}

result_folders <-  data.frame(path = result_folders) %>% mutate(sample = basename(path))

metadata <- read.csv(metadata_fp) %>% rename(sample =Sample)
metadata$Confidence <- factor(metadata$Confidence, levels = c('High confidence', 'Low confidence', 'Benign'))
metadata$PFI <- factor(metadata$PFI, levels = c('Short', 'Long', 'Benign'))

metadata <- left_join(metadata, result_folders, by = 'sample')


# get interactions --------------------------------------------------------

#BG vs. HC
ROI.interactions <- get_differential_interactions(metadata, 'Confidence', levels(metadata %>% pull(matches('Confidence')))[1:2])

#long vs short PFI in HC
responder.interactions <- get_differential_interactions(metadata %>% filter(Confidence == 'High confidence'), 'PFI', levels(metadata %>% pull(matches('PFI')))[1:2])




