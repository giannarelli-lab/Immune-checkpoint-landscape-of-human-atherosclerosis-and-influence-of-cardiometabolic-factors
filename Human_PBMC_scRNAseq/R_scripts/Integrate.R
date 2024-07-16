#### Ravneet Kaur 
## Sep 15 2022

library("Seurat")

load("~/PBMCs_CHORD/sc_analysis/pbmc_obj_comb_genes_rna.rda")

pbmc_obj_list <- SplitObject(pbmc_obj_comb_genes, split.by = "orig.ident")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmc_obj_list <- lapply(X = pbmc_obj_list, FUN = function(x) {
  x <- SCTransform(x, assay = 'RNA', new.assay.name = 'SCT',variable.features.n = 3000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE,assay='SCT')
  x <- SCTransform(x, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("S.Score", "G2M.Score"),variable.features.n = 3000)
})

DefaultAssay(pbmc_obj_list[[1]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[2]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[3]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[4]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[5]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[6]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[7]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[8]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[9]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[10]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[11]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[12]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[13]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[14]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[15]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[16]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[17]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[18]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[19]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[20]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[21]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[22]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[23]]) <- "SCT"
DefaultAssay(pbmc_obj_list[[24]]) <- "SCT"


features <- SelectIntegrationFeatures(pbmc_obj_list, nfeatures = 3000)
pbmc_obj_list <- PrepSCTIntegration(pbmc_obj_list, anchor.features = features)

pbmc_obj_list <- lapply(X = pbmc_obj_list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors  <- FindIntegrationAnchors(pbmc_obj_list, normalization.method = "SCT", anchor.features = features, dims = 1:30, reference = c(1,3,7,15,17), reduction = "rpca")
pbmc_obj_doub_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:30)

save(pbmc_obj_doub_integrated,file="~/PBMCs_CHORD/sc_analysis/pbmc_obj_doub_integrated.rda")

# Export for celltypist
SaveH5Seurat(pbmc_obj_rm_doub_integrated, filename = "pbmc_obj_rm_doub_integrated.h5Seurat")
Convert("pbmc_obj_rm_doub_integrated.h5Seurat", dest = "pbmc_obj_rm_doub_RNA.h5ad")
Convert("pbmc_obj_rm_doub_integrated.h5Seurat", dest = "pbmc_obj_rm_doub_integrated.h5ad", assay = 'integrated')

