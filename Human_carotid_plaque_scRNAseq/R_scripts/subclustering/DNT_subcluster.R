library(Seurat)
library(CVRCFunc)

double_neg <- readRDS('T_cells/DNT.rds')

DefaultAssay(double_neg) <- "integrated"
double_neg <- RunPCA(double_neg, verbose = FALSE, npcs = 30)
double_neg <- RunUMAP(double_neg, dims = 1:30, verbose = TRUE)
double_neg <- FindNeighbors(double_neg, verbose = TRUE, dims = 1:30)
double_neg <- FindClusters(double_neg, verbose = TRUE)
double_neg <- NormalizeData(double_neg, assay = 'RNA')

double_neg <- ScaleData(double_neg, assay = 'RNA')
FindMarkersBulk(double_neg, clus_ident = 'seurat_clusters', sample_ident = 'new.ident')
