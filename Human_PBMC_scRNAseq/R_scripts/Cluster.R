#### Ravneet Kaur 
library(Seurat)

load("~/PBMCs_CHORD/sc_analysis/pbmc_obj_doub_integrated.rda")


# Remove doublets
Idents(pbmc_obj_doub_integrated)<-"DoubletClass"
pbmc_obj_rm_doub_integrated <-subset(pbmc_obj_doub_integrated,idents=c("Singlet"))
pbmc_obj_rm_doub_integrated


# Scale data
pbmc_obj_rm_doub_integrated <- ScaleData(pbmc_obj_rm_doub_integrated,assay="RNA")

# Select the RNA counts slot to be the default assay
DefaultAssay(pbmc_obj_rm_doub_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
pbmc_obj_rm_doub_integrated<- NormalizeData(pbmc_obj_rm_doub_integrated, verbose = FALSE)
DefaultAssay(pbmc_obj_rm_doub_integrated) <- "integrated"
pbmc_obj_rm_doub_integrated <- RunPCA(pbmc_obj_rm_doub_integrated, npcs = 30)
pbmc_obj_rm_doub_integrated <- RunUMAP(pbmc_obj_rm_doub_integrated, reduction = "pca", dims = 1:30) 
pbmc_obj_rm_doub_integrated <- FindNeighbors(pbmc_obj_rm_doub_integrated, reduction = "pca", dims = 1:30)
pbmc_obj_rm_doub_integrated <- FindClusters(pbmc_obj_rm_doub_integrated, resolution = 0.7)


saveRDS(pbmc_obj_rm_doub_integrated,file="pbmc_obj_rm_doub_integrated.rds")


Collapse










