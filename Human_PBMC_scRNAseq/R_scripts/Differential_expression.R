library(Seurat)
library(CVRCFunc)

DefaultAssay(pbmc_obj_rm_doub_integrated) <- "RNA"
pbmc_obj_rm_doub_integrated <- NormalizeData(pbmc_obj_rm_doub_integrated, assay = 'RNA')

#Pseudobulk DE annotation
FindMarkersBulk(pbmc_obj_rm_doub_integrated, clus_ident = 'annotation_fine', sample_ident = 'orig.ident', assay = 'RNA', out_dir = 'FindMarkersBulk_out_annotation_fine')


###Type 2 vs no diabetes in NT
### Wilcox
Idents(pbmc_obj_rm_doub_integrated) <- 'Condition'
pbmc_obj_rm_doub_integrated_bl <- subset(pbmc_obj_rm_doub_integrated, idents = 'NT')
Idents(pbmc_obj_rm_doub_integrated_bl) <- 'annotation_fine'
sig_up <- vector()
sig_down <- vector()
for(i in unique(pbmc_obj_rm_doub_integrated_bl$annotation_fine)){
  tr <- tryCatch(markers <- Seurat::FindMarkers(pbmc_obj_rm_doub_integrated_bl, ident.1 = 'DM2', ident.2 ='CONTROL', group.by = 'Diagnosis', subset.ident = i, assay = 'RNA', slot = 'data'),
                 error = function(e) e) 
  if(inherits(tr, "error")){
    markers <- data.frame(row.names = ICI_targets_all, 
                          p_val = rep(1, length(ICI_targets_all)),
                          avg_log2FC = rep(0, length(ICI_targets_all)),
                          pct.1 = rep(0, length(ICI_targets_all)),
                          pct.2 = rep(1, length(ICI_targets_all)),
                          p_val_adj = rep(1, length(ICI_targets_all)))
  }
  write.csv(markers, file = paste('celltypist_wilcox_diabetes_bl/',i,'.csv',sep = ''), quote = F)
  sig_up <- c(sig_up, length(which(markers$p_val_adj < 0.1 & markers$avg_log2FC > 0)))
  sig_down <- c(sig_down, length(which(markers$p_val_adj < 0.1 & markers$avg_log2FC < 0)))
}
summary_results <- data.frame(cluster = unique(pbmc_obj_rm_doub_integrated_bl$annotation_fine), sig_up = sig_up, sig_down = sig_down, sig_all = sig_up+sig_down)
write.csv(file = paste('celltypist_wilcox_diabetes_bl/summary_results.csv'), summary_results, quote = F, row.names = F)
rm(pbmc_obj_rm_doub_integrated_bl)

