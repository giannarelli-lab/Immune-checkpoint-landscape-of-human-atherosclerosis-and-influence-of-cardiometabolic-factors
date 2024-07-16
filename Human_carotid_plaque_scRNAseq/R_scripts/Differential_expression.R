library(Seurat)
library(CVRCFunc)


#refined file curated by Gabriel 05232023
ICI_targets <-read.csv('ICI_target_pairs_refined.csv', header = T)
ICI_targets_all <- unique(c(ICI_targets$Receptors, ICI_targets$Ligand))

Whole_dataset_annotated_human.rds <- readRDS(file = 'Whole_dataset_annotated_human.rds')
Whole_dataset_annotated_human.rds <- ScaleData(Whole_dataset_annotated_human.rds, assay = 'RNA')

#Pseudobulk find markers
FindMarkersBulk(Whole_dataset_annotated_human.rds, clus_ident = 'annotation_fine', sample_ident = 'new.ident', out_dir = 'annotation_fine_findmarkersout')

#Wilcoxon differential expression type 2 vs no diabetes
Idents(Whole_dataset_annotated_human.rds) <- 'annotation_fine'
#diabetes
sig_up <- vector()
sig_down <- vector()
for(i in unique(Whole_dataset_annotated_human.rds$annotation_fine)){
  tr <- tryCatch(markers <- Seurat::FindMarkers(Whole_dataset_annotated_human.rds, ident.1 = 'Type 2', ident.2 ='No Diabetes', group.by = 'conditions', subset.ident = i, assay = 'RNA', slot = 'data'),
                 error = function(e) e) 
  if(inherits(tr, "error")){
    markers <- data.frame(row.names = ICI_targets_all, 
                          p_val = rep(1, length(ICI_targets_all)),
                          avg_log2FC = rep(0, length(ICI_targets_all)),
                          pct.1 = rep(0, length(ICI_targets_all)),
                          pct.2 = rep(1, length(ICI_targets_all)),
                          p_val_adj = rep(1, length(ICI_targets_all)))
  }
  write.csv(markers, file = paste('annotation_fine_wilcox_diabetes/',i,'.csv',sep = ''), quote = F)
  sig_up <- c(sig_up, length(which(markers$p_val_adj < 0.1 & markers$avg_log2FC > 0)))
  sig_down <- c(sig_down, length(which(markers$p_val_adj < 0.1 & markers$avg_log2FC < 0)))
}
summary_results <- data.frame(cluster = unique(Whole_dataset_annotated_human.rds$annotation_fine), sig_up = sig_up, sig_down = sig_down, sig_all = sig_up+sig_down)
write.csv(file = paste('annotation_fine_wilcox_diabetes/summary_results.csv'), summary_results, quote = F, row.names = F)