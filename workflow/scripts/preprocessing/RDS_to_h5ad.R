library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

if(exists("snakemake")){
  input_fp = normalizePath(snakemake@input$integrated)
  output_fp <- snakemake@output[[1]]
}else{
  input_fp = normalizePath("results/integrated/plate_int.RDS")
  output_fp <- "data/working/MO/MO_brain_annotated"
}

data <- readRDS(input_fp)
  
if(dir.exists(output_fp)){
  unlink(output_fp, recursive=TRUE)
}
dir.create(output_fp)

assays <- base::intersect(names(data@assays), c('Spatial', 'SCT'))

lapply(assays, function (x){
  cat('INFO: converting assay:', as.character(x), '\n')
  temp <- data
  temp@active.assay <- x
  keep <- grepl(paste('^',paste(c(paste(paste0(x, c('_pca', '_harmony', '_umap'))), 'pca', 'harmony', 'umap'), collapse = '$|^'), '$', sep = ''), names(temp@reductions))
  temp@reductions[!keep] <- NULL
  
  cat('INFO: with reductions', names(temp@reductions), '\n')

  temp <- as.SingleCellExperiment(temp)
  
  if(exists('snakemake')){
    writeH5AD(temp, file = paste(output_fp, .Platform$file.sep, gsub('_', '-', x), '.h5ad', sep = '' )) 
    
  }
  
})
  
