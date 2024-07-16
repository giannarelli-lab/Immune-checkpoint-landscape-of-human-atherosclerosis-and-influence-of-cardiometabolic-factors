#### Ravneet Kaur 
library(Seurat)

load("/gpfs/data/giannarellilab/PBMCs_CHORD/sc_analysis/pbmc_obj_list_rm_doublet.rda")

pbmc_obj_comb <- merge(pbmc_obj_list[[1]], y = pbmc_obj_list[2:length(pbmc_obj_list)])

save(pbmc_obj_comb,file="/gpfs/data/giannarellilab/PBMCs_CHORD/sc_analysis/pbmc_obj_comb_rna.rda")

genes<-read.csv("/gpfs/data/giannarellilab/PBMCs_CHORD/sc_analysis/rm_sample/cell_cycle_module_genes.csv")

counts <- GetAssayData(pbmc_obj_comb, assay = "RNA")
g<-genes$cell_cycle_module

counts <- counts[-(which(rownames(counts) %in% c(g))),]

#And then subset your data based upon the genes (row names)
pbmc_obj_comb_genes <- subset(pbmc_obj_comb, features = rownames(counts))

save(pbmc_obj_comb_genes,file="/gpfs/data/giannarellilab/PBMCs_CHORD/sc_analysis/pbmc_obj_comb_genes_rna.rda")