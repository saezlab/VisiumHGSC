library(tidyr)
library(dplyr)
library(ggplot2)
library(mistyR)

# helpers -----------------------------------------------------------------

# reformat_samples <- function(misty.results, view){
#     
#     results <- lapply(misty.results, function(x){
#         if('sample' %in% colnames(x)){
#             x$sample <- x$sample %>% basename() 
#         }
#         return(x)
#     })
#     
#     return(results)
# }
# 
# extract_contrast_interactions <- function (misty.results.from, misty.results.to, 
#                                            views = NULL, cutoff.from = 1, cutoff.to = 1, 
#                                            trim = -Inf, trim.measure = c(
#                                                "gain.R2", "multi.R2","intra.R2",
#                                                "gain.RMSE", "multi.RMSE", "intra.RMSE"
#                                            )){
#     
#     trim.measure.type <- match.arg(trim.measure)
#     
#     #check that both results collections are properly formatted
#     assertthat::assert_that(("importances.aggregated" %in% names(misty.results.from)), 
#                             msg = "The first provided result list is malformed. Consider using collect_results().")
#     assertthat::assert_that(("improvements.stats" %in% names(misty.results.from)), 
#                             msg = "The provided result list is malformed. Consider using collect_results().")
#     assertthat::assert_that(("importances.aggregated" %in% names(misty.results.to)), 
#                             msg = "The second provided result list is malformed. Consider using collect_results().")
#     assertthat::assert_that(("improvements.stats" %in% names(misty.results.to)), 
#                             msg = "The provided result list is malformed. Consider using collect_results().")
#     
#     
#     if (is.null(views)) {
#         #check that both result collections contain the same views
#         assertthat::assert_that(rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::pull(.data$view) %>% unique(),
#                                                         misty.results.to$importances.aggregated %>% dplyr::pull(.data$view) %>% unique())),
#                                 msg = "The requested views do not exist in both result lists.")
#         
#         views <- misty.results.from$importances.aggregated %>% dplyr::pull(.data$view) %>% unique()
#     } else {
#         #check that all views are present in both result collections
#         assertthat::assert_that(all(views %in% (misty.results.from$importances.aggregated %>% dplyr::pull(.data$view))) &
#                                     all(views %in% (misty.results.to$importances.aggregated %>% dplyr::pull(.data$view))),
#                                 msg = "The requested views do not exist in both result lists.")
#         
#     }
#     
#     #check that both collections contain exactly the same target-predictor interactions in each view
#     assertthat::assert_that(all(views %>% purrr::map_lgl(function(current.view) {
#         
#         rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
#                                     dplyr::pull(.data$Predictor) %>% unique(),
#                                 misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
#                                     dplyr::pull(.data$Predictor) %>% unique())) &
#             
#             rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
#                                         dplyr::pull(.data$Target) %>% unique(),
#                                     misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
#                                         dplyr::pull(.data$Target) %>% unique()))
#         
#     })), msg = "Incompatible predictors and targets.")
#     
#     
#     inv <- sign((stringr::str_detect(trim.measure.type, "gain") |
#                      stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5)
#     
#     
#     targets <- misty.results.from$improvements.stats %>% dplyr::filter(.data$measure == trim.measure.type, inv * .data$mean >= inv * trim) %>% dplyr::pull(.data$target)
#     
#     interactions <- views %>% purrr::map_dfr(function(current.view) {
#         
#         from.view.wide <- misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view, .data$Target %in% targets) %>%
#             tidyr::pivot_wider(names_from = "Target", values_from = "Importance", -c(.data$view, .data$nsamples))
#         
#         to.view.wide <- misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view, .data$Target %in% targets) %>%
#             tidyr::pivot_wider(names_from = "Target", values_from = "Importance", -c(.data$view, .data$nsamples))
#         
#         
#         
#         
#         mask <- ((from.view.wide %>% dplyr::select(-.data$Predictor)) < cutoff.from) & ((to.view.wide %>% dplyr::select(-.data$Predictor)) >= cutoff.to)
#         
#         
#         if(sum(mask, na.rm = TRUE) > 0){ 
#             # assertthat::assert_that(sum(mask, na.rm = TRUE) > 0, msg = paste0("All values are cut off while contrasting."))
#             
#             
#             masked <- ((to.view.wide %>% tibble::column_to_rownames("Predictor")) * mask)
#             
#             
#             masked %>% dplyr::slice(which(masked %>% rowSums(na.rm = TRUE) > 0)) %>%
#                 dplyr::select(which(masked %>% colSums(na.rm = TRUE) > 0)) %>% tibble::rownames_to_column("Predictor") %>%
#                 tidyr::pivot_longer(names_to = "Target", values_to = "Importance", -.data$Predictor) %>% dplyr::filter(Importance > 0) %>%
#                 dplyr::mutate(view = current.view)
#         }else{
#             mat = matrix(ncol = 0, nrow = 0)
#             df=data.frame(mat)
#         }
#         
#     })
#     
#     return(interactions)
# }
# 
# get_differential_interactions <- function(metadata, grouping_var, groups, cutoff.from = 1, cutoff.to = 1, output.importances = FALSE){
#     
#     grouped.results <- lapply(groups, function(group){
#         group.results <- metadata %>% filter(!!as.symbol(grouping_var) == group) %>% pull(path) %>%
#             collect_results() %>% reformat_samples()
#     })
#     names(grouped.results) <- groups
#     
#     if(output.importances) imp <- list()
#     
#     
#     contrast.interactions <- names(grouped.results) %>%  purrr::map_dfr(function(to.group){
#         
#         from.group <- names(grouped.results)[which(!grepl(to.group, names(grouped.results)))]
#         
#         # get interactions that are only in to.group, but not from.group
#         # using default cutoff of 1 (on Importance), as 'being present'
#         interactions <- extract_contrast_interactions(grouped.results[[from.group]], grouped.results[[to.group]], cutoff.from = cutoff.from, cutoff.to = cutoff.to) %>% 
#             tidyr::unite(col = 'Interaction', .data$view, .data$Predictor, .data$Target, sep = '_')
#         
#         # extract the importances per sample for these interactions and add metadata
#         importances <- bind_rows(list(grouped.results[[1]]$importances, grouped.results[[2]]$importances)) %>% tidyr::unite(col = 'Interaction', .data$view, .data$Predictor, .data$Target, remove = FALSE) %>% 
#             dplyr::filter(.data$Interaction %in% (interactions %>% dplyr::pull(.data$Interaction))) %>% dplyr::select(-.data$Predictor, -.data$Target) %>% dplyr::left_join(metadata, by = 'sample')
#         
#         if(output.importances) imp <<- append(imp, list(importances))
#         
#         # t-test over conditions and do BH p.value adjustment
#         stats <- importances %>% dplyr::group_by(.data$Interaction) %>% rstatix::t_test(data =., as.formula(paste('Importance', '~', grouping_var))) %>% 
#             dplyr::left_join(importances %>% dplyr::select(.data$Interaction, .data$view) %>% dplyr::distinct(), by = 'Interaction') %>% dplyr::group_by(.data$view) %>% 
#             rstatix::adjust_pvalue(method = "BH") %>% dplyr::select(.data$view, .data$Interaction, .data$statistic, .data$p, .data$p.adj) %>% dplyr::rename(t.value = .data$statistic, p.value = .data$p) %>% dplyr::ungroup()
#         
#         stats %>% plyr::adply(.margins = 1, function(x){
#             x$Interaction <- base::gsub(paste('^', x$view, '_' , sep=''), '', x$Interaction)
#             return(x) 
#         }) %>% tidyr::separate(col=.data$Interaction, into = c('Predictor', 'Target'), sep = '_') %>% dplyr::mutate(only.in = to.group)
#         
#     })
#     
#     if(output.importances){
#         imp <- dplyr::bind_rows(imp)
#         return(list(interactions = contrast.interactions, importances = imp))
#     }else{
#         return(contrast.interactions)    
#     }
#     
# }
# 
# 

# set input and output paths ----------------------------------------------

cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  
  tissue <- snakemake@wildcards$tissue
  model <- snakemake@wildcards$view_type
  contrast <- snakemake@wildcards$contrast
  
  helpers_fp <- snakemake@input[[1]]
  metadata_fp <- snakemake@input[[2]]
  result_folders <- unlist(snakemake@input[3:length(snakemake@input)])
  
  importances_fp <- snakemake@output[[1]]
  interactions_fp <- snakemake@output[[2]]
  
}else{
  tissue <- 'brain'
  model <- 'celltype'
  
  helpers_fp <- 'workflow/helpers/misty.R'
  
  samples <- list.files(paste('data/original/ST/visium_data', tissue, sep = '_')) %>% sort()
  
  result_folders <- paste0(paste('results/ST/Misty', tissue, sep = .Platform$file.sep), '/', samples, '/', model, '_misty_model')
  
  #files for testing in Rstudio
  metadata_fp <- paste('data/original/ST/metadata_visium_', tissue,'.csv', sep = '')
  
}

source(helpers_fp)

# load data ---------------------------------------------------------------

result_folders <-  data.frame(path = result_folders) %>% mutate(sample = basename(path))

metadata <- read.csv(metadata_fp) %>% rename(sample =Sample)
metadata$Confidence <- factor(metadata$Confidence, levels = c('High confidence', 'Low confidence', 'Benign'))
metadata$PFI <- factor(metadata$PFI, levels = c('Short', 'Long', 'Benign'))

metadata <- left_join(metadata, result_folders, by = 'sample')

results <- metadata %>% pull(path) %>% collect_results() %>% reformat_samples()


# define contrasts --------------------------------------------------------

if(contrast == 'HCvsBG'){
    grouping_var <- 'Confidence'
    groups <-levels(metadata %>% pull(matches('Confidence')))[1:2]
}else if (contrast == 'ShortvsLong'){
    metadata <- metadata %>% filter(Confidence == 'High confidence')
    grouping_var <- 'PFI'
    groups <- levels(metadata %>% pull(matches('PFI')))[1:2]
}
  
# grouped results ---------------------------------------------------------

inter <- get_differential_interactions(metadata, grouping_var, groups, correction = 'Target', cutoff.from = 1, cutoff.to = 1, output.importances = TRUE)
    
if(exists("snakemake")){
  
  write.csv(inter$importances %>% mutate(model = model, contrast = contrast), file = importances_fp, sep = ",")
  write.csv(inter$interactions %>% mutate(model = model, contrast = contrast), file = interactions_fp, sep = ",")
  
}


