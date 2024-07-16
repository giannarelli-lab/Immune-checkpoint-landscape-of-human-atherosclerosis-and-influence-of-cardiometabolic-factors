library(Seurat)
library(CVRCFunc)

CD4 <- readRDS(file = "T_cells/CD4.rds")

# Remove samples less with 100 cells and cluster 14 (only present in 1 sample and markers don't make sense with T-cells).
CD4_stats <- as.data.frame(table(CD4$new.ident))
Idents(CD4) <- 'new.ident'
CD4 <- subset(CD4, idents = c('GSM7018579_Sample1','GSM7018581_Sample3','GSM7018583_Sample4'), invert = T)
Idents(CD4) <- 'seurat_clusters'
CD4 <- subset(CD4, idents = c('14'), invert = T)

DefaultAssay(CD4) <- "integrated"
CD4 <- RunPCA(CD4, verbose = FALSE, npcs = 30)
CD4 <- RunUMAP(CD4, dims = 1:30, verbose = TRUE)
CD4 <- FindNeighbors(CD4, verbose = TRUE, dims = 1:30)
CD4 <- FindClusters(CD4, verbose = TRUE, resolution = .83)

CD4 <- ScaleData(CD4, assay = 'RNA')
FindMarkersBulk(seurat = CD4, clus_ident = 'annotation_fine', sample_ident = 'new.ident', alpha = 0.05)
