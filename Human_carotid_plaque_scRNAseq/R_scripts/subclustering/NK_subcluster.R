library(Seurat)
library(CVRCFunc)

load('NK-cell_subset.rda')
DefaultAssay(NK_cells) <- "integrated"
NK_cells <- RunPCA(NK_cells, verbose = FALSE, npcs = 30)
NK_cells <- RunUMAP(NK_cells, dims = 1:30, verbose = TRUE)
NK_cells <- FindNeighbors(NK_cells, verbose = TRUE, dims = 1:30)
NK_cells <- FindClusters(NK_cells, verbose = TRUE, resolution = 0.3)
