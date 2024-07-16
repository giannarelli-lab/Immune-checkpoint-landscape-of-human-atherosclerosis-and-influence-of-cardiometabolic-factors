# Create a list of seurat objects for each sample
library(Seurat)

df_samples <- read.csv("samples.csv", row.names = 1)
samples <- df_samples$sample
diabetic_raw <- list()
diabetic_Seurat <- list()
for(i in 1:length(samples)){
  diabetic_raw[[i]] <- Read10X(data.dir = paste0("/gpfs/data/giannarellilab/scRNAseq_FastQ_Data/", samples[[i]],"/filtered_gene_bc_matrices/"))
  colnames(diabetic_raw[[i]]) <- paste0(samples[i], "_",colnames(diabetic_raw[[i]]))
  diabetic_Seurat[[i]] <- CreateSeuratObject(diabetic_raw[[i]], min.cells = 0, min.features = 0, names.delim = "_")
  diabetic_Seurat[[i]][["percent.mt"]] <- PercentageFeatureSet(diabetic_Seurat[[i]], pattern = "^MT-")
}
saveRDS(diabetic_Seurat, 'diabetic_Seurat_04202021.rds')


