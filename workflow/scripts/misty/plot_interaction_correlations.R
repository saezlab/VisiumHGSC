library(mistyR)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)


# files -------------------------------------------------------------------

if(exists("snakemake")){
    #files and parameters input by snakemake
    
    view_type <- snakemake@wildcards$view_type
    
    metadata_fp <- snakemake@input$metadata
    interactions_fp <- snakemake@input$interactions
    correlations_fp <- snakemake@input$correlations
    
    output_fp <- snakemake@output[[1]]
    
    
}else{
    #files for testing in Rstudio
    
    view_type <- 'celltype'
    
    metadata_fp <- 'results/integrated/sample_metadata.csv'
    interactions_fp <-  str_glue('results/Misty/{view_type}/diffInteractions.csv', view_type=view_type, .sep = "")
    correlations_fp <-  str_glue('results/Misty/{view_type}/model_correlations.csv', view_type=view_type, .sep = "")
    
    
}


# load data ---------------------------------------------------------------

metadata <- read.csv(metadata_fp)
interactions <- read.csv(interactions_fp)

title_suffix <- ''

if((interactions %>% filter(p.adj <= 0.05) %>% nrow()) < 1){
  interactions <- interactions %>% arrange(p.adj) %>% head(4)
  title_suffix <- '\n(only n.s.)'
}else{
  interactions <- interactions %>% filter(p.adj <= 0.05)
}


 # select only significant interactions
correlations <- read.csv(correlations_fp)


# making plots ------------------------------------------------------------

if(exists("snakemake")) pdf(snakemake@output[[1]])

#iterate over contrasts
interactions %>% pull(contrast) %>% unique() %>% walk(function(contrast){
    
    tit <- 'Correlations for differential interactions'
    
    #data selection based on contrast of interest
    if(contrast == 'HCvsBG'){
        meta_df <- metadata %>% filter(Confidence != 'Benign')
        grouping_var <- 'Confidence'
        groups <- c('High confidence', 'Low confidence')
    }else if (contrast == 'ShortvsLong'){
        meta_df <- metadata %>% filter(Confidence == 'High confidence')
        grouping_var <- 'PFI'
        groups <- c('Short', 'Long')
        tit <- paste(tit, '\nin high confidence cores')
    }
    
    #combine all three dataframes
    to.plot <- left_join(interactions %>% select(view, Predictor, Target), correlations %>% filter(sample %in% meta_df$Sample)) %>% 
        left_join(meta_df, by =c('sample' = 'Sample')) %>%  unite(col = interaction, view, Predictor, sep = ": ") %>% 
        unite(col = interaction, interaction, Target, sep = " -> ")
    print(groups)
    print(to.plot)
    to.plot[grouping_var] <- factor(to.plot %>% pull(grouping_var), levels = groups)
    
    if(to.plot %>% nrow() > 0){
        #make boxplot
        plot <- ggplot(data = to.plot, mapping = aes(y = corr, x = !!as.symbol(grouping_var), fill = !!as.symbol(grouping_var))) + geom_boxplot(outlier.shape = NA) + 
            geom_jitter(width = 0.3, height = 0) + facet_wrap(~interaction, scale="free_x") + theme_bw() +
            theme(axis.title.x=element_blank()) +
            labs(y = 'Correlation', title = paste(tit, title_suffix, sep = ''))
        print(plot)
    }
})

if(exists("snakemake")) dev.off()










