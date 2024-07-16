library(Seurat)
library(CVRCFunc)
library(SeuratDisk)


pbmc_obj_rm_doub_integrated <- readRDS('pbmc_obj_rm_doub_integrated.rds')
DefaultAssay(pbmc_obj_rm_doub_integrated) <- "RNA"
pbmc_obj_rm_doub_integrated <- NormalizeData(pbmc_obj_rm_doub_integrated, assay = 'RNA')

celltypist <- read.csv('celltypist_labels.csv')
celltypist$majority_voting <- gsub(celltypist$majority_voting, pattern = '/', replacement = '_')
pbmc_obj_rm_doub_integrated$celltypist_majvote <- celltypist$majority_voting
pbmc_obj_rm_doub_integrated$celltypist <- celltypist$predicted_labels
saveRDS(pbmc_obj_rm_doub_integrated, file = 'pbmc_obj_rm_doub_integrated.rds')

# Carotid subclustering/annotation strategy
## T_cells
Idents(pbmc_obj_rm_doub_integrated) <- 'seurat_clusters'
T_cells <- subset(pbmc_obj_rm_doub_integrated, idents = c('18', '5', '16', '9', '7', '10', '1', '12', '21', '2', '3', '6', '13', '27', '19'))
DimPlot(T_cells, reduction = 'umap', label = T)

CD4_cells <- WhichCells(T_cells, expression = CD4 > 0 & CD8A == 0 & CD8B == 0)
CD4 <- subset(T_cells, cells=CD4_cells)
saveRDS(CD4, file = 'CD4.rds')

CD8_cells <- WhichCells(T_cells, expression = CD4 == 0 & CD8A + CD8B > 0)
CD8 <- subset(T_cells, cells=CD8_cells)
saveRDS(CD8, file = 'CD8.rds')

DNT <- WhichCells(T_cells, expression = CD4 == 0 & CD8A + CD8B == 0)
DNT <- subset(T_cells, cells=DNT)
saveRDS(DNT, file = 'DNT.rds')

DPT <- WhichCells(T_cells, expression = CD4 > 0 & CD8A + CD8B > 0)
DPT <- subset(T_cells, cells=DPT)
saveRDS(DPT, file = 'DPT.rds')

## NK
Idents(pbmc_obj_rm_doub_integrated) <- 'seurat_clusters'
NK_cells <- subset(pbmc_obj_rm_doub_integrated, idents = c('20', '4'))
DimPlot(NK_cells, reduction = 'umap', label = T)
saveRDS(NK_cells, file = 'NK.rds')

## B
Idents(pbmc_obj_rm_doub_integrated) <- 'seurat_clusters'
B_cells <- subset(pbmc_obj_rm_doub_integrated, idents = c('8', '15', '11', '17', '25', '29'))
DimPlot(B_cells, reduction = 'umap', label = T)
saveRDS(B_cells, file = 'B.rds')

#include other plasma cell cluster (06092023)
Idents(pbmc_obj_rm_doub_integrated) <- 'seurat_clusters'
B_cells <- subset(pbmc_obj_rm_doub_integrated, idents = c('8', '15', '11', '17', '25', '29','30'))
DimPlot(B_cells, reduction = 'umap', label = T)
saveRDS(B_cells, file = 'B_updated.rds')

## Myeloid
Idents(pbmc_obj_rm_doub_integrated) <- 'seurat_clusters'
myeloid_cells <- subset(pbmc_obj_rm_doub_integrated, idents = c('0', '22', '24', '14', '26'))
DimPlot(myeloid_cells, reduction = 'umap', label = T)
saveRDS(myeloid_cells, file = 'Myeloid.rds')

#Ran this code after initial rounds of annotation
#re extract CD4 based on manual annotation. Some B cells clusters were re-asasigned
Idents(pbmc_obj_rm_doub_integrated) <- 'annotation_major'
CD4 <- subset(pbmc_obj_rm_doub_integrated, idents = 'CD4 T cell')
saveRDS(CD4, file = 'CD4_updated.rds')
#re extract DNT based on manual annotation. Some B cells clusters were re-asasigned
Idents(pbmc_obj_rm_doub_integrated) <- 'annotation_major'
DNT <- subset(pbmc_obj_rm_doub_integrated, idents = 'DN T cell')
saveRDS(DNT, file = 'DNT_updated.rds')
