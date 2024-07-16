library(Seurat)
library(SeuratDisk)

load("final_object_list.rda")
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

# Apply additional filtering based on downstream analysis
final_obj_list$GSM8393427_CEA_scRNAseq_3 <-  subset(final_obj_list$GSM8393427_CEA_scRNAseq_3, subset = nFeature_RNA > 200 & nCount_RNA < 16000 & nCount_RNA > 1500 & percent.mt < 10)
final_obj_list$GSM8393426_CEA_scRNAseq_2 <-  subset(final_obj_list$GSM8393426_CEA_scRNAseq_2, subset = nFeature_RNA > 200 & nCount_RNA < 16000 & nCount_RNA > 1500 & percent.mt < 10)


#Annotate
anchors <- list()
for (i in 1:length(final_obj_list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = final_obj_list[[i]],
    k.filter = 200,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
}


for (i in 1:length(final_obj_list)) {
  final_obj_list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = final_obj_list[[i]],
    reference = reference, 
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}


#Integrate
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
final_obj_list <- lapply(X = final_obj_list, FUN = function(x) {
  x <- SCTransform(x, assay = 'RNA', new.assay.name = 'SCT',variable.features.n = 3000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE,assay='SCT')
  x <- SCTransform(x, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("S.Score", "G2M.Score"),variable.features.n = 3000)
})
DefaultAssay(final_obj_list[[1]]) <- "SCT"
DefaultAssay(final_obj_list[[2]]) <- "SCT"
DefaultAssay(final_obj_list[[3]]) <- "SCT"
DefaultAssay(final_obj_list[[4]]) <- "SCT"
DefaultAssay(final_obj_list[[5]]) <- "SCT"
DefaultAssay(final_obj_list[[6]]) <- "SCT"
DefaultAssay(final_obj_list[[7]]) <- "SCT"
DefaultAssay(final_obj_list[[8]]) <- "SCT"
DefaultAssay(final_obj_list[[9]]) <- "SCT"
DefaultAssay(final_obj_list[[10]]) <- "SCT"
DefaultAssay(final_obj_list[[11]]) <- "SCT"
DefaultAssay(final_obj_list[[12]]) <- "SCT"
DefaultAssay(final_obj_list[[13]]) <- "SCT"
DefaultAssay(final_obj_list[[14]]) <- "SCT"
DefaultAssay(final_obj_list[[15]]) <- "SCT"
DefaultAssay(final_obj_list[[16]]) <- "SCT"
DefaultAssay(final_obj_list[[17]]) <- "SCT"
DefaultAssay(final_obj_list[[18]]) <- "SCT"
DefaultAssay(final_obj_list[[19]]) <- "SCT"
DefaultAssay(final_obj_list[[20]]) <- "SCT"
DefaultAssay(final_obj_list[[21]]) <- "SCT"
DefaultAssay(final_obj_list[[22]]) <- "SCT"
DefaultAssay(final_obj_list[[23]]) <- "SCT"
features <- SelectIntegrationFeatures(final_obj_list, nfeatures = 3000)
final_obj_list <- PrepSCTIntegration(final_obj_list, anchor.features = features)
final_obj_list <- lapply(X = final_obj_list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors  <- FindIntegrationAnchors(final_obj_list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reference = c('GSM8393427_CEA_scRNAseq_3','GSM8393428_CEA_scRNAseq_4','GSM8393432_CEA_scRNAseq_8','GSM8393436_CEA_scRNAseq_12'), reduction = "rpca")
final_obj_integration_rev <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:30)

#Fix some annotations
final_obj_integration_rev$new.ident <- final_obj_integration_rev$orig.ident

#Combine GSM7018581_Sample3 and GSM7018581_Sample3a samples # From the same patient
final_obj_integration_rev@meta.data$new.ident[which(final_obj_integration_rev@meta.data$new.ident == 'GSM7018581_Sample3')] = "GSM7018581_Sample3"
final_obj_integration_rev@meta.data$new.ident[which(final_obj_integration_rev@meta.data$new.ident == 'GSM7018581_Sample3a')] = "GSM7018581_Sample3"

#Add metadata
metadata <- read.csv('metadata_updated_02092023.csv')
metadata <- metadata[which(metadata$sample %in% unique(final_obj_integration_rev$orig.ident)),]
conditions <- final_obj_integration_rev$orig.ident
#Conditions (diabetes)
for(i in 1:nrow(metadata)){
  conditions[which(conditions == metadata$sample[i])] <- metadata$conditions[i] 
}
final_obj_integration_rev$conditions <- conditions

save(final_obj_integration_rev,file="final_obj_integration_cellcycle_reg.rda")
