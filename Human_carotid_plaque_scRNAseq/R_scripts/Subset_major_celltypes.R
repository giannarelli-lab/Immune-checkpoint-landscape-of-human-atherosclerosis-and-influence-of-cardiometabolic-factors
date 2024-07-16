library(Seurat)

load("final_obj_integration_cellcycle_reg_clus.rda")
Idents(final_obj_integration_rev) <- "seurat_clusters"

Myeloid_cells <- subset(final_obj_integration_rev, idents = c(8,10,15,17))
save(Myeloid_cells, file = "Myeloid-cells/Myeloid-cell_subset.rda")

B_cells <- subset(final_obj_integration_rev, idents = c(7,14))
save(B_cells, file = "B-cells/B-cell_subset.rda")

NK_cells <- subset(final_obj_integration_rev, idents = c(3,5,19))
save(NK_cells, file = "NK-cells/NK-cell_subset.rda")

T_cells <- subset(final_obj_integration_rev, idents = c(1,9,0,18,12,11,16,2,6,4))
save(T_cells, file = "T_cells/T-cell_subset.rda")

CD4_cells <- WhichCells(T_cells, expression = CD4 > 0 & CD8A == 0 & CD8B == 0)
CD4 <- subset(T_cells, cells=CD4_cells)
saveRDS(CD4, file = 'T_cells/CD4.rds')

CD8_cells <- WhichCells(T_cells, expression = CD4 == 0 & CD8A + CD8B > 0)
CD8 <- subset(T_cells, cells=CD8_cells)
saveRDS(CD8, file = 'T_cells/CD8.rds')

DNT <- WhichCells(T_cells, expression = CD4 == 0 & CD8A + CD8B == 0)
DNT <- subset(T_cells, cells=DNT)
saveRDS(DNT, file = 'T_cells/DNT.rds')

DPT <- WhichCells(T_cells, expression = CD4 > 0 & CD8A + CD8B > 0)
DPT <- subset(T_cells, cells=DPT)
saveRDS(DPT, file = 'T_cells/DPT.rds')