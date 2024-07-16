library(Seurat)
library(CVRCFunc)

load('B-cell_subset.rda')
DefaultAssay(B_cells) <- "integrated"
B_cells <- ScaleData(B_cells, assay = 'integrated')
B_cells <- RunPCA(B_cells, verbose = FALSE, npcs = 30)
B_cells <- RunUMAP(B_cells, dims = 1:30, verbose = TRUE)
B_cells <- FindNeighbors(B_cells, verbose = TRUE, dims = 1:30)
B_cells <- FindClusters(B_cells, verbose = TRUE, resolution = .2)

B_cells <- ScaleData(B_cells, assay = 'RNA')
FindMarkersBulk(B_cells, clus_ident = 'seurat_clusters', sample_ident = 'new.ident')
