# library(mistyR)
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
    
    sig.cutoff <- snakemake@params$sig_cutoff
    
    metadata_fp <- snakemake@input$metadata
    interactions_fp <- snakemake@input$interactions
    correlations_fp <- snakemake@input$correlations
    
    output_fp <- snakemake@output[[1]]
    
    
}else{
    #files for testing in Rstudio
    
    view_type <- 'celltype'
    
    sig.cutoff <- 0.05
    
    metadata_fp <- 'results/integrated/sample_metadata.csv'
    interactions_fp <-  str_glue('results/Misty/{view_type}/diffInteractions.csv', view_type=view_type, .sep = "")
    correlations_fp <-  str_glue('results/Misty/{view_type}/model_correlations.csv', view_type=view_type, .sep = "")
    
    
}


# load data ---------------------------------------------------------------

metadata <- read.csv(metadata_fp)
interactions <- read.csv(interactions_fp)


 # select only significant interactions
correlations <- read.csv(correlations_fp)


# making plots ------------------------------------------------------------

if(exists("snakemake")) pdf(snakemake@output[[1]])

#iterate over contrasts
interactions %>% pull(contrast) %>% unique() %>% walk(function(c){
  
    contrast.inter <- interactions %>% filter(contrast == c)
    
    if((contrast.inter %>% filter(p.adj <= sig.cutoff) %>% nrow()) < 1){
      contrast.inter <- contrast.inter %>% arrange(p.adj) %>% head(4)
      title_suffix <- '\n(only n.s.)'
    }else{
      contrast.inter <- contrast.inter %>% filter(p.adj <= sig.cutoff)
      title_suffix <- ''
    }
    
    tit <- 'Correlations for differential interactions'
    
    #data selection based on contrast of interest
    if(c == 'HCvsBG'){
        meta_df <- metadata %>% filter(Confidence != 'Benign')
        grouping_var <- 'Confidence'
        groups <- c('High confidence', 'Low confidence')
        tit <- paste(tit, '\nin high vs. low confidence cores')
    }else if (c == 'ShortvsLongHC'){
        meta_df <- metadata %>% filter(Confidence == 'High confidence')
        grouping_var <- 'PFI'
        groups <- c('Short', 'Long')
        tit <- paste(tit, '\nin high confidence cores')
    }else if (c == 'ShortvsLongBG'){
      meta_df <- metadata %>% filter(Confidence == 'Low confidence')
      grouping_var <- 'PFI'
      groups <- c('Short', 'Long')
      tit <- paste(tit, '\nin low confidence cores')
    }
    
    #combine all three dataframes
    to.plot <- left_join(contrast.inter %>% select(view, Predictor, Target, p.adj), correlations %>% filter(sample %in% meta_df$Sample)) %>% 
        left_join(meta_df, by =c('sample' = 'Sample')) %>%  unite(col = interaction, view, Predictor, sep = ": ") %>% 
        unite(col = interaction, interaction, Target, sep = " -> ") %>% unite(col = interaction, interaction, p.adj, sep = '\np.adj = ', remove = FALSE)
    to.plot[grouping_var] <- factor(to.plot %>% pull(grouping_var), levels = groups)
    to.plot$interaction <- factor(to.plot$interaction, levels = to.plot %>% arrange(p.adj) %>% pull(interaction) %>% unique())
    
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










