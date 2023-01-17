# load libraries ----------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)

# set input output files --------------------------------------------------

cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  
  density_fp <- snakemake@input[[1]]
  metadata <- snakemake@input[[2]]
  
}else{
  
  density_fp <- 'data/cell_per_location.csv'
  metadata_fp <- 'results/integrated/sample_metadata.csv'
  
}


# load data ---------------------------------------------------------------

density <- read.csv(density_fp) %>% filter(!is.na(cell_per_location))
metadata <- read.csv(metadata_fp)

density <- density %>%  left_join(metadata, by ='Sample')


# plot --------------------------------------------------------------------

if(exists("snakemake")) pdf(snakemake@output[[1]], width = 7, height = 7)

density %>% ggplot(aes(x = PFI, y = cell_per_location, fill = PFI)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0) + facet_wrap(~Confidence, scale = 'free_x') +
  labs(title = 'Cell density across conditions', y = 'Cells per spot') + theme_bw() +
  theme(axis.title.x=element_blank()) 
  
if(exists("snakemake")) dev.off()