library(mistyR)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)


# files -------------------------------------------------------------------

if(exists("snakemake")){
  #files and parameters input by snakemake
  
  view_type <- snakemake@wildcards$view_type
  sample <- snakemake@wildcards$sample
  
  correlation.type <- snakemake@params$corr
  
  interactions_fp <- snakemake@input$interactions
  view_fp <- snakemake@input$view

  output_fp <- snakemake@output[[1]]
  
  
}else{
  #files for testing in Rstudio
  
  view_type <- 'celltype'
  sample <- 'BG-L_OVA15-16'
  
  correlation.type <- 'pearson'
  
  interactions_fp <- str_glue('results/Misty/{view_type}/diffInteractions.csv', view_type=view_type, .sep = "")
  view_fp <- str_glue('results/Misty/{view_type}/views/{sample}_view.rds', sample=sample, view_type=view_type, .sep = "")
  

}


# load data ---------------------------------------------------------------

interactions <- read.csv(interactions_fp) %>% select(view, Target, Predictor) %>% dplyr::distinct(.keep_all = TRUE)

misty.views <- readRDS(view_fp) %>% lapply(function(view){
    if(class(view) != "character") colnames(view$data) <- gsub('\\.|_','',colnames(view$data))
    return(view)
})



# do correlation spot by spot ---------------------------------------------

targets <- intersect(interactions$Target, colnames(misty.views$intraview$data))


output <- targets %>% map_dfr(function(target){
  
  current.inter <- interactions %>% filter(.data$Target == target)
  
  # if(view_type == 'celltype' | view_type == 'CTpathways'){
  #   temp.views <- misty.views %>% filter_views(NA, view = 'mask', .data[[target]]) %>% remove_views('mask') 
  # }else{
    temp.views <- misty.views
  # }
  
  inter_cor <- current.inter %>% pull(.data$view) %>% unique() %>% map_dfr(function(current.view){
    
    #select columns of target-predictors
    if(current.view == 'intra'){
      df <- temp.views$intraview$data
    }else if (current.view == 'intra_act'){
      df <- cbind(temp.views$intraview$data %>% select(all_of(target)),
                  temp.views$intra_act$data %>% select(any_of(current.inter %>%
                                                                filter(.data$view == current.view) %>% pull(.data$Predictor) %>% unique())))
    }else if (current.view == 'para'){
      df <- cbind(temp.views$intraview$data %>% select(all_of(target)),
                  temp.views$paraview$data %>% select(any_of(current.inter %>%
                                                                filter(.data$view == current.view) %>% pull(.data$Predictor) %>% unique())))
    }
  
    correlations <- df %>%
      cor(method = correlation.type) %>%
      as.data.frame %>%
      rownames_to_column(var = 'var1') %>%
      gather(var2, value, -var1) %>% rename(corr = value)
    
    inner_join(current.inter %>% filter(.data$view == current.view), correlations, by = c('Predictor' = 'var1', 'Target' = 'var2'))
    
  })
  
  
  
  
}) %>% mutate(sample = sample)

if(exists("snakemake")){
  write.csv(output, file = output_fp, sep = ",")
}
