library(tidyverse)
library(Seurat)
library(harmony)

data.dir <- normalizePath("data/sc-data")
counts.fp <- normalizePath(paste(data.dir, 'GSE165897_UMIcounts_HGSOC.tsv', sep = .Platform$file.sep))
meta.fp <- normalizePath(paste(data.dir, 'GSE165897_cellInfo_HGSOC.tsv', sep = .Platform$file.sep))

counts <- read.csv(counts.fp, sep = '\t', header = TRUE)
meta <- read.csv(meta.fp, sep = '\t', header = TRUE) %>% mutate(cell = gsub('-', '.', cell))

counts <- counts %>% tibble::column_to_rownames(var = 'X')
meta <- meta %>% tibble::column_to_rownames(var = 'cell')

ref.data <- CreateSeuratObject(counts = counts, project = 'Ovisium_ref')

ref.data <- AddMetaData(ref.data, metadata = meta)

rm(counts, meta)

# eda ---------------------------------------------------------------------

ref.data <- subset(ref.data, subset = treatment_phase == 'treatment-naive')


VlnPlot(ref.data, features = c('nFeature_RNA'), group.by ='patient_id')



# clustering --------------------------------------------------------------

ref.data <- NormalizeData(ref.data)

ref.data <- FindVariableFeatures(ref.data, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(ref.data)
ref.data <- ScaleData(ref.data, features = all.genes)

ref.data <- RunPCA(ref.data, features = VariableFeatures(object = ref.data))

ref.data <- RunHarmony(ref.data, c("patient_id"), assay.use = "RNA", reduction.save = "harmony", max.iter.harmony = 20)

ref.data <- FindNeighbors(ref.data, reduction = "harmony")

ref.data <- RunUMAP(ref.data, reduction = "harmony", dims = 1:30)


DimPlot(ref.data, reduction = "umap", group.by = 'cell_subtype')

system('say Robin, the program has finished running, get your ass back here and work!')
