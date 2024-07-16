library(Seurat)
library(ggplot2)
library(presto)
library(viridis)
library(pheatmap)
library(Hmisc)
library(CVRCFunc)
library(stringr)
library(RColorBrewer)
library(gridExtra)

seurat <- readRDS('/gpfs/data/moorelab/Gabriel/scRNAseq/PlaqueCD45/IntegratedMousescRNAseq/data.integrated.ref.rds')
DefaultAssay(seurat) <- "integrated"
seurat <- RunPCA(seurat, verbose = FALSE, npcs = 30)
seurat <- RunUMAP(seurat, dims = 1:30, verbose = TRUE)
seurat <- FindNeighbors(seurat, verbose = TRUE, dims = 1:30)
seurat <- FindClusters(seurat, verbose = TRUE, resolution = .75)

#have to run on HPC
anchors <- c('Mki67', 'Pcna', 'Mcm3', 'Top2a', 'Ccnb2')
anchor_ind <- which(row.names(seurat@assays$RNA@data) %in% anchors)
names(anchor_ind) <- row.names(seurat@assays$RNA@data)[anchor_ind]

cor_mat <- matrix(dimnames = list(row.names(seurat@assays$RNA@data), names(anchor_ind)), nrow = nrow(seurat@assays$RNA@data), ncol = 5)
for(i in 1:nrow(cor_mat)){
  for(j in 1:ncol(cor_mat)){
    cor_mat[i,j] <- cor(seurat@assays$RNA@data[i,], seurat@assays$RNA@data[anchor_ind[j],]) 
  }
}
saveRDS(cor_mat, file = 'cor_mat_mouse.rds')


cor_mat <- readRDS('cor_mat.rds')
putative_cc_genes <- cor_mat[which(apply(cor_mat, MARGIN = 1, function(x) max(x)) > 0.15),]
dim(putative_cc_genes)
head(cor_mat)
hist(putative_cc_genes, breaks = 20, col = 'steel blue', xlab = 'Pearson correlation coefficient')
pheatmap(putative_cc_genes)
write.csv(row.names(putative_cc_genes), file = 'putative_cc_genes_mouse.csv')

