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
interactions <- read.csv(interactions_fp) %>% filter(p.value <= 0.05) # select only significant interactions
correlations <- read.csv(correlations_fp)


# making plots ------------------------------------------------------------

if(exists("snakemake")) pdf(snakemake@output[[1]])

#iterate over contrasts
interactions %>% pull(contrast) %>% unique() %>% walk(function(contrast){
    
    tit <- 'Correlations for differential interactions'
    
    #data selection based on contrast of interest
    if(contrast == 'HCvsBG'){
        meta_df <- metadata
        grouping_var <- 'Confidence'
        groups <-levels(metadata %>% pull(matches('Confidence')))[1:2]
    }else if (contrast == 'LvsS'){
        meta_df <- metadata %>% filter(Confidence == 'High confidence')
        grouping_var <- 'PFI'
        groups <- c('Short', 'Long')
        tit <- paste(tit, '\nin high confidence cores')
    }
    
    #combine all three dataframes
    to.plot <- left_join(interactions %>% select(view, Predictor, Target), correlations %>% filter(sample %in% meta_df$Sample)) %>% 
        left_join(meta_df, by =c('sample' = 'Sample')) %>%  unite(col = interaction, view, Predictor, sep = ": ") %>% 
        unite(col = interaction, interaction, Target, sep = " -> ")
    
    to.plot[grouping_var] <- factor(to.plot %>% pull(grouping_var), levels = groups)
    
    #make boxplot
    plot <- ggplot(data = to.plot, mapping = aes(y = corr, x = interaction, fill = !!as.symbol(grouping_var))) + geom_boxplot(outlier.shape = NA) + 
        facet_wrap(~interaction, scale="free") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        labs(y = 'Correlation', title = tit)
    print(plot)
})

if(exists("snakemake")) dev.off()









