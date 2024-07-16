library(Seurat)

load('final_obj_integration_cellcycle_reg.rda')
DefaultAssay(final_obj_integration_rev) <- "integrated"
final_obj_integration_rev <- RunPCA(final_obj_integration_rev, npcs = 30)
final_obj_integration_rev <- RunUMAP(final_obj_integration_rev, reduction = "pca", dims=1:30)
final_obj_integration_rev <- FindNeighbors(final_obj_integration_rev, verbose = FALSE, dims = 1:30)
final_obj_integration_rev <- FindClusters(final_obj_integration_rev, verbose = FALSE, resolution = .5)
save(final_obj_integration_rev, file = "final_obj_integration_cellcycle_reg_clus.rda")