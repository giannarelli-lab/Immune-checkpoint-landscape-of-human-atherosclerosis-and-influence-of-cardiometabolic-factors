library(Seurat)
library(CVRCFunc)

double_pos <- readRDS('T_cells/DPT.rds')

DefaultAssay(double_pos) <- "integrated"
double_pos <- RunPCA(double_pos, verbose = FALSE, npcs = 30)
double_pos <- RunUMAP(double_pos, dims = 1:30, verbose = TRUE)
double_pos <- FindNeighbors(double_pos, verbose = TRUE, dims = 1:30)
double_pos <- FindClusters(double_pos, verbose = TRUE)
double_pos <- NormalizeData(double_pos, assay = 'RNA')

double_pos <- ScaleData(double_pos, assay = 'RNA')
FindMarkersBulk(double_pos, clus_ident = 'seurat_clusters', sample_ident = 'new.ident')
