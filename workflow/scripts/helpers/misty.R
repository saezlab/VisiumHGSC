library(tidyr)
library(dplyr)
library(ggplot2)
library(mistyR)

options(dplyr.legacy_locale = TRUE)

reformat_samples <- function(misty.results, view){
  
  results <- lapply(misty.results, function(x){
    if('sample' %in% colnames(x)){
      x$sample <- x$sample %>% basename() 
    }
    return(x)
  })
  
  return(results)
}

extract_contrast_interactions <- function (misty.results.from, misty.results.to, 
                                           views = NULL, cutoff.from = 1, cutoff.to = 1, 
                                           trim = -Inf, trim.measure = c(
                                             "gain.R2", "multi.R2","intra.R2",
                                             "gain.RMSE", "multi.RMSE", "intra.RMSE"
                                           )){
  
  trim.measure.type <- match.arg(trim.measure)
  
  #check that both results collections are properly formatted
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results.from)), 
                          msg = "The first provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that(("improvements.stats" %in% names(misty.results.from)), 
                          msg = "The provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results.to)), 
                          msg = "The second provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that(("improvements.stats" %in% names(misty.results.to)), 
                          msg = "The provided result list is malformed. Consider using collect_results().")
  
  
  if (is.null(views)) {
    #check that both result collections contain the same views
    assertthat::assert_that(rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::pull(.data$view) %>% unique(),
                                                    misty.results.to$importances.aggregated %>% dplyr::pull(.data$view) %>% unique())),
                            msg = "The requested views do not exist in both result lists.")
    
    views <- misty.results.from$importances.aggregated %>% dplyr::pull(.data$view) %>% unique()
  } else {
    #check that all views are present in both result collections
    assertthat::assert_that(all(views %in% (misty.results.from$importances.aggregated %>% dplyr::pull(.data$view))) &
                              all(views %in% (misty.results.to$importances.aggregated %>% dplyr::pull(.data$view))),
                            msg = "The requested views do not exist in both result lists.")
    
  }
  
  #check that both collections contain exactly the same target-predictor interactions in each view
  assertthat::assert_that(all(views %>% purrr::map_lgl(function(current.view) {
    
    rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                              dplyr::pull(.data$Predictor) %>% unique(),
                            misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                              dplyr::pull(.data$Predictor) %>% unique())) &
      
      rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                                dplyr::pull(.data$Target) %>% unique(),
                              misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                                dplyr::pull(.data$Target) %>% unique()))
    
  })), msg = "Incompatible predictors and targets.")
  
  
  inv <- sign((stringr::str_detect(trim.measure.type, "gain") |
                 stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5)
  
  
  targets <- misty.results.from$improvements.stats %>% dplyr::filter(.data$measure == trim.measure.type, inv * .data$mean >= inv * trim) %>% dplyr::pull(.data$target)
  
  interactions <- views %>% purrr::map_dfr(function(current.view) {
    
    from.view.wide <- misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view, .data$Target %in% targets) %>%
      tidyr::pivot_wider(names_from = "Target", values_from = "Importance", -c(.data$view, .data$nsamples))
    
    to.view.wide <- misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view, .data$Target %in% targets) %>%
      tidyr::pivot_wider(names_from = "Target", values_from = "Importance", -c(.data$view, .data$nsamples))
    
    
    
    
    mask <- ((from.view.wide %>% dplyr::select(-.data$Predictor)) < cutoff.from) & ((to.view.wide %>% dplyr::select(-.data$Predictor)) >= cutoff.to)
    
    
    if(sum(mask, na.rm = TRUE) > 0){ 
      # assertthat::assert_that(sum(mask, na.rm = TRUE) > 0, msg = paste0("All values are cut off while contrasting."))
      
      
      masked <- ((to.view.wide %>% tibble::column_to_rownames("Predictor")) * mask)
      
      
      masked %>% dplyr::slice(which(masked %>% rowSums(na.rm = TRUE) > 0)) %>%
        dplyr::select(which(masked %>% colSums(na.rm = TRUE) > 0)) %>% tibble::rownames_to_column("Predictor") %>%
        tidyr::pivot_longer(names_to = "Target", values_to = "Importance", -.data$Predictor) %>% dplyr::filter(Importance > 0) %>%
        dplyr::mutate(view = current.view)
    }else{
      mat = matrix(ncol = 0, nrow = 0)
      df=data.frame(mat)
    }
    
  })
  
  return(interactions)
}

get_differential_interactions <- function(metadata, grouping_var, groups, cutoff.from = 1, cutoff.to = 1, correction = 'Target', output.importances = FALSE){
  
  grouped.results <- lapply(groups, function(group){
    group.results <- metadata %>% filter(!!as.symbol(grouping_var) == group) %>% pull(path) %>%
      collect_results() %>% reformat_samples()
  })
  names(grouped.results) <- groups
  
  if(output.importances) imp <- list()
  
  
  contrast.interactions <- names(grouped.results) %>%  purrr::map_dfr(function(to.group){
    
    from.group <- names(grouped.results)[which(!grepl(to.group, names(grouped.results)))]
    
    # get interactions that are only in to.group, but not from.group
    # using default cutoff of 1 (on Importance), as 'being present'
    interactions <- extract_contrast_interactions(grouped.results[[from.group]], grouped.results[[to.group]], cutoff.from = cutoff.from, cutoff.to = cutoff.to) %>% 
      tidyr::unite(col = 'Interaction', .data$view, .data$Predictor, .data$Target, sep = '_')
    
    # extract the importances per sample for these interactions and add metadata
    importances <- bind_rows(list(grouped.results[[1]]$importances, grouped.results[[2]]$importances)) %>% tidyr::unite(col = 'Interaction', .data$view, .data$Predictor, .data$Target, remove = FALSE) %>% 
      dplyr::filter(.data$Interaction %in% (interactions %>% dplyr::pull(.data$Interaction))) %>% dplyr::select(-.data$Predictor, -.data$Target) %>% dplyr::left_join(metadata, by = 'sample')
    
    if(output.importances) imp <<- append(imp, list(importances))
    
    # t-test over conditions and do BH p.value adjustment
    stats <- importances %>% dplyr::group_by(.data$Interaction) %>% rstatix::t_test(data =., as.formula(paste('Importance', '~', grouping_var))) %>% 
      dplyr::left_join(importances %>% dplyr::select(.data$Interaction, .data$view) %>% dplyr::distinct(), by = 'Interaction') %>%
      dplyr::select(.data$view, .data$Interaction, .data$statistic, .data$p) %>% dplyr::rename(t.value = .data$statistic, p.value = .data$p) %>% dplyr::ungroup()
    
    stats %>% plyr::adply(.margins = 1, function(x){
      x$Interaction <- base::gsub(paste('^', x$view, '_' , sep=''), '', x$Interaction)
      return(x) 
    }) %>% tidyr::separate(col=.data$Interaction, into = c('Predictor', 'Target'), sep = '_') %>% dplyr::mutate(only.in = to.group)
    
  })
  print(colnames(contrast.interactions))
  contrast.interactions <- contrast.interactions %>% dplyr::group_by(!!as.symbol(correction)) %>% rstatix::adjust_pvalue(method = "BH", p.col = 'p.value', output.col = 'p.adj')
  
  
  if(output.importances){
      imp <- dplyr::bind_rows(imp)
      return(list(interactions = contrast.interactions, importances = imp))
  }else{
      return(contrast.interactions)
  }
  
}

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

#' Convenient functions
new_scale_fill <- function() {
  new_scale("fill")
}

new_scale_color <- function() {
  new_scale("colour")
}

new_scale_colour <- function() {
  new_scale("colour")
}

#' Special behaviour of the "+" for adding a `new_aes` object
#' It changes the name of the aesthethic for the previous layers, appending
#' "_new" to them. 
ggplot_add.new_aes <- function(object, plot, object_name) {
  plot$layers <- lapply(plot$layers, bump_aes, new_aes = object)
  plot$scales$scales <- lapply(plot$scales$scales, bump_aes, new_aes = object)
  plot$labels <- bump_aes(plot$labels, new_aes = object)
  plot
}


bump_aes <- function(layer, new_aes) {
  UseMethod("bump_aes")
}

bump_aes.Scale <- function(layer, new_aes) {
  old_aes <- layer$aesthetics[remove_new(layer$aesthetics) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  layer$aesthetics[layer$aesthetics %in% old_aes] <- new_aes
  
  if (is.character(layer$guide)) {
    layer$guide <- match.fun(paste("guide_", layer$guide, sep = ""))()
  }
  layer$guide$available_aes[layer$guide$available_aes %in% old_aes] <- new_aes
  layer
}

bump_aes.Layer <- function(layer, new_aes) {
  original_aes <- new_aes
  
  old_aes <- names(layer$mapping)[remove_new(names(layer$mapping)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  old_geom <- layer$geom
  
  old_setup <- old_geom$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup(data, params)
  }
  
  new_geom <- ggplot2::ggproto(paste0("New", class(old_geom)[1]), old_geom,
                               handle_na = new_setup)
  
  new_geom$default_aes <- change_name(new_geom$default_aes, old_aes, new_aes)
  new_geom$non_missing_aes <- change_name(new_geom$non_missing_aes, old_aes, new_aes)
  new_geom$required_aes <- change_name(new_geom$required_aes, old_aes, new_aes)
  new_geom$optional_aes <- change_name(new_geom$optional_aes, old_aes, new_aes)
  
  layer$geom <- new_geom
  
  old_stat <- layer$stat
  
  old_setup2 <- old_stat$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup2(data, params)
  }
  
  new_stat <- ggplot2::ggproto(paste0("New", class(old_stat)[1]), old_stat,
                               handle_na = new_setup)
  
  new_stat$default_aes <- change_name(new_stat$default_aes, old_aes, new_aes)
  new_stat$non_missing_aes <- change_name(new_stat$non_missing_aes, old_aes, new_aes)
  new_stat$required_aes <- change_name(new_stat$required_aes, old_aes, new_aes)
  new_stat$optional_aes <- change_name(new_stat$optional_aes, old_aes, new_aes)
  
  layer$stat <- new_stat
  
  layer$mapping <- change_name(layer$mapping, old_aes, new_aes)
  layer
}

bump_aes.list <- function(layer, new_aes) {
  old_aes <-  names(layer)[remove_new(names(layer)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  names(layer)[names(layer) %in% old_aes] <- new_aes
  layer
}

change_name <- function(list, old, new) {
  UseMethod("change_name")
}

change_name.character <- function(list, old, new) {
  list[list %in% old] <- new
  list
}

change_name.default <- function(list, old, new) {
  nam <- names(list)
  nam[nam %in% old] <- new
  names(list) <- nam
  list
}

change_name.NULL <- function(list, old, new) {
  NULL
}

remove_new <- function(aes) {
  stringi::stri_replace_all(aes, "", regex = "(_new)*")
}

