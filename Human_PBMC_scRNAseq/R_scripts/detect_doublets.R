### Doublet Finder on PBMC chord data
#### Ravneet Kaur 
#### 05/10/2021

library(Seurat)
library(DoubletFinder)


load("pbmc_obj_list.rda")

for(i in 1:length(pbmc_obj_list)){
  
  # sample name 
  pbmc.list <- pbmc_obj_list[i]
  filename <- names(pbmc.list)
  
  print(filename)
  
  pbmc.seurat <- pbmc.list[[1]]
  doublet.seurat <- pbmc.list[[1]]
  
  ## Pre-process Seurat object (standard) 
  doublet.seurat <- NormalizeData(doublet.seurat)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat<- SCTransform(doublet.seurat, assay = 'RNA', new.assay.name = 'SCT',variable.features.n = 3000)
  doublet.seurat <- RunPCA(doublet.seurat)
  
  ## pK Identification
  sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:30, sct =TRUE)
  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]
  pK <- as.numeric(levels(pK))[pK]
  
  EDFR=0.075 
  homotypic.prop <- modelHomotypic(doublet.seurat@meta.data$predicted.celltype.l2)
  nExp_poi <- round(EDFR*nrow(doublet.seurat@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)
  
  class <- paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')
  doublet.seurat@meta.data[['DoubletClass']] <- doublet.seurat@meta.data[[class]]
  
  ## Put doublet annotations into "original.seurat"
  original.seurat@meta.data[['DoubletClass']] <- doublet.seurat@meta.data[['DoubletClass']]
  
  ## REMOVE DOUBLETS in original.seurat, I didn't do this because I wanted to visualize the doublets in a UMAP.
  pbmc_obj_list[[i]] <- original.seurat
  
}

save(pbmc_obj_list,file="pbmc_obj_list_rm_doublet.rda")



