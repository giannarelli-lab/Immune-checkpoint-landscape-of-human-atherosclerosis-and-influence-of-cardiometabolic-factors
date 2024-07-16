library(Seurat)
library(CVRCFunc)

load('Myeloid-cell_subset.rda')

#Subcluster
#remove sample CV7209 very few cells
Idents(Myeloid_cells) <- 'new.ident'
Myeloid_cells <- subset(Myeloid_cells, idents = 'GSM7018579_Sample1', invert = T)
DefaultAssay(Myeloid_cells) <- "integrated"
Myeloid_cells <- RunPCA(Myeloid_cells, verbose = FALSE, npcs = 30)
Myeloid_cells <- RunUMAP(Myeloid_cells, dims = 1:30, verbose = TRUE)
Myeloid_cells <- FindNeighbors(Myeloid_cells, verbose = TRUE, dims = 1:30)
Myeloid_cells <- FindClusters(Myeloid_cells, verbose = TRUE)

Myeloid_cells <- ScaleData(Myeloid_cells, assay = 'RNA')
FindMarkersBulk(Myeloid_cells, clus_ident = 'seurat_clusters', sample_ident = 'new.ident')
