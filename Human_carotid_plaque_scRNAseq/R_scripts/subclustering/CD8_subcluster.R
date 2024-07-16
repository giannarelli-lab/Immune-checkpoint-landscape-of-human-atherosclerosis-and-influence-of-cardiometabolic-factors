library(Seurat)
library(CVRCFunc)

CD8_Tcells <- readRDS('T_cells/CD8.rds')

Idents(CD8_Tcells) <- 'new.ident'
CD8_Tcells <- subset(CD8_Tcells, idents = 'GSM7018579_Sample1', invert = T)

DefaultAssay(CD8_Tcells) <- "integrated"
CD8_Tcells <- RunPCA(CD8_Tcells, verbose = FALSE, npcs = 30)
CD8_Tcells <- RunUMAP(CD8_Tcells, dims = 1:30, verbose = TRUE)
CD8_Tcells <- FindNeighbors(CD8_Tcells, verbose = TRUE, dims = 1:30)
CD8_Tcells <- FindClusters(CD8_Tcells, verbose = TRUE, resolution = .75)

CD8_Tcells <- ScaleData(CD8_Tcells, assay = 'RNA')
FindMarkersBulk(CD8_Tcells, clus_ident = 'seurat_clusters', sample_ident = 'new.ident')
