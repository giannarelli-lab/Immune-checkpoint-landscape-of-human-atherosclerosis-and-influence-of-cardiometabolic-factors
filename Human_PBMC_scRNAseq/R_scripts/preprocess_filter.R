#### Ravneet Kaur 
library(Seurat)
library(dplyr)
library(ggplot2)

all_samples <- read.csv("samples.csv")
samples <- all_samples$sample
Diagnosis <- all_samples$Diagnosis
Condition <- all_samples$Condition
pbmc_raw <- list()
pbmc_seurat <- list()
for(i in 1:length(samples)){
  pbmc_raw[[i]] <- Read10X(data.dir = paste0("cellranger/",
                                             samples[[i]],"/outs/filtered_feature_bc_matrix/"))
  colnames(pbmc_raw[[i]]) <- paste0(samples[i],
                                    "_",colnames(pbmc_raw[[i]]))
  pbmc_seurat[[i]] <- CreateSeuratObject(pbmc_raw[[i]],
                                         min.cells = 3,
                                         min.features = 200,
                                         names.delim = "_")
  pbmc_seurat[[i]]@meta.data$Condition <- Condition[i]
  pbmc_seurat[[i]]@meta.data$Diagnosis <- Diagnosis[i]
}

head(pbmc_seurat[[1]])

pbmc <- Reduce(function(x, y) merge(x, y, do.normalize = F), pbmc_seurat)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

head(pbmc)

### Rename nCount_RNA and nFeature_RNA slots

pbmc$num_UMIs = pbmc$nCount_RNA
pbmc$num_genes = pbmc$nFeature_RNA

message("Number of cells:", {ncol(pbmc)})

message("Number of genes:", {nrow(pbmc)})

message("unfiltered mean num genes:", {round(mean(pbmc$num_genes), 3)})
message("unfiltered median num genes:", {median(pbmc$num_genes)})
message("unfiltered mean num UMIs:", {round(mean(pbmc$num_UMIs), 3)})
message("unfiltered median num UMIs:", {median(pbmc$num_UMIs)})
message("unfiltered mean percent mito:", {round(mean(pbmc$percent.mt), 3)})
message("unfiltered median percent mito:", {round(median(pbmc$percent.mt),3)})

# Call function for volin plots
g1 <- lapply(c("num_genes", "num_UMIs", "percent.mt"), function(features){
  VlnPlot(object = pbmc, features = features, ncol = 3, pt.size = 0)+
    theme(axis.text.x = element_text(size=15),legend.position="none")
})

options(repr.plot.width=50,repr.plot.height=20)
print(plot_grid(g1[[1]]+ggtitle("num_genes_RNA before filteration")))#+
# scale_y_log10(limits = c(50,9000))))

options(repr.plot.width=30,repr.plot.height=10)
print(plot_grid(g1[[2]]+ggtitle("num_UMIs_RNA before filteration")))#+
#scale_y_log10(limits = c(500,100000))))

options(repr.plot.width=30,repr.plot.height=10)

print(plot_grid(g1[[3]]+ggtitle("percent mito % before filteration")))#+
# ylim(c(0,40))))

cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc, feature1 = "num_UMIs", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')

cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc, feature1 = "num_genes", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc, feature1 = "num_UMIs", feature2 = "num_genes",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

pbmc_seurat <- lapply(X = pbmc_seurat, FUN = function(x) {
  x<-PercentageFeatureSet(object = x, pattern = "^MT-", col.name = "percent.mt") 
})

pbmc_seurat_filt<-pbmc_seurat

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_1
pbmc_seurat_filt[[1]]
head(pbmc_seurat_filt[[1]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[1]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[1]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[1]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[1]] <- subset(pbmc_seurat_filt[[1]], subset =  nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

#QC plots after fileration
options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[1]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[1]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[1]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")
print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

message("Number of cells:", {ncol(pbmc_seurat_filt[[1]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[1]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_1
pbmc_seurat_filt[[2]]
head(pbmc_seurat_filt[[2]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[2]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[2]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[2]] <- subset(pbmc_seurat_filt[[2]], subset =  nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)


#QC plots after fileration
options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[2]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[2]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")
print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

message("Number of cells:", {ncol(pbmc_seurat_filt[[2]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[2]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_2
pbmc_seurat_filt[[3]]
head(pbmc_seurat_filt[[3]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[3]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[3]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[3]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[3]] <- subset(pbmc_seurat_filt[[3]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

#QC plots after fileration
options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[3]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[3]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[3]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")
print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

message("Number of cells:", {ncol(pbmc_seurat_filt[[3]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[3]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_2
pbmc_seurat_filt[[4]]
head(pbmc_seurat_filt[[4]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[4]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[4]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[4]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[4]] <- subset(pbmc_seurat_filt[[4]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

#QC plots after fileration
options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[4]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[4]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[4]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")
print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

message("Number of cells:", {ncol(pbmc_seurat_filt[[4]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[4]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_3
pbmc_seurat_filt[[5]]
head(pbmc_seurat_filt[[5]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[5]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[5]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[5]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[5]] <- subset(pbmc_seurat_filt[[5]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

#QC plots after fileration
options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[5]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[5]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[5]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")
print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

message("Number of cells:", {ncol(pbmc_seurat_filt[[5]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[5]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_3
pbmc_seurat_filt[[6]]
head(pbmc_seurat_filt[[6]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[6]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[6]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[6]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[6]] <- subset(pbmc_seurat_filt[[6]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[6]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[6]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[6]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

pbmc_seurat_filt[[6]]

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_4
pbmc_seurat_filt[[7]]
head(pbmc_seurat_filt[[7]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[7]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[7]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[7]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[7]] <- subset(pbmc_seurat_filt[[7]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[7]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[7]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[7]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))
dim(pbmc_seurat_filt[[7]])

pbmc_seurat_filt[[7]]

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_4
pbmc_seurat_filt[[8]]
head(pbmc_seurat_filt[[8]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[8]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[8]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[8]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[8]] <- subset(pbmc_seurat_filt[[8]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[8]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[8]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[8]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[8]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[8]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[8]])})

pbmc_seurat_filt[[9]]<-pbmc_seurat[[9]]

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_5
pbmc_seurat_filt[[9]]
head(pbmc_seurat_filt[[9]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[9]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[9]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[9]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[9]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[9]] <- subset(pbmc_seurat_filt[[9]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[9]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[9]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[9]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[9]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

pbmc_seurat_filt[[9]]

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_5
pbmc_seurat_filt[[10]]
head(pbmc_seurat_filt[[10]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[10]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[10]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[10]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[10]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[10]] <- subset(pbmc_seurat_filt[[10]], subset = nFeature_RNA > 200 &  nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[10]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[10]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[10]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[10]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[10]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[10]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_6
pbmc_seurat_filt[[11]]
head(pbmc_seurat_filt[[11]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[11]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[11]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[11]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[11]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[11]] <- subset(pbmc_seurat_filt[[11]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[11]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[11]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[11]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[11]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[11]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[11]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_6
pbmc_seurat_filt[[12]]
head(pbmc_seurat_filt[[12]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[12]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[12]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[12]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[12]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[12]] <- subset(pbmc_seurat_filt[[12]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[12]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[12]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[12]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[12]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[12]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[12]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_7
pbmc_seurat_filt[[13]]
head(pbmc_seurat_filt[[13]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[13]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[13]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[13]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[13]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[13]] <- subset(pbmc_seurat_filt[[13]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[13]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[13]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[13]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[13]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[13]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[13]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_7
pbmc_seurat_filt[[14]]
head(pbmc_seurat_filt[[14]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[14]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[14]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[14]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[14]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[14]] <- subset(pbmc_seurat_filt[[14]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[14]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[14]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[14]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[14]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[14]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[14]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_8
pbmc_seurat_filt[[15]]
head(pbmc_seurat_filt[[15]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[15]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[15]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[15]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[15]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[15]] <- subset(pbmc_seurat_filt[[15]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[15]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[15]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[15]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[15]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[15]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[15]])})

pbmc_seurat_filt[[16]]<-pbmc_seurat[[16]]

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_8
pbmc_seurat_filt[[16]]
head(pbmc_seurat_filt[[16]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[16]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[16]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[16]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[16]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[16]] <- subset(pbmc_seurat_filt[[16]], subset = nFeature_RNA > 200 &  nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[16]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[16]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[16]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[16]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[16]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[16]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_9
pbmc_seurat_filt[[17]]
head(pbmc_seurat_filt[[17]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[17]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[17]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[17]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[17]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[17]] <- subset(pbmc_seurat_filt[[17]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[17]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[17]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[17]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[17]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[17]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[17]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_9
pbmc_seurat_filt[[18]]
head(pbmc_seurat_filt[[18]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[18]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[18]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[18]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[18]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[18]] <- subset(pbmc_seurat_filt[[18]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[18]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[18]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[18]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[18]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[18]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[18]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_10
pbmc_seurat_filt[[19]]
head(pbmc_seurat_filt[[19]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[19]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[19]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[19]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[19]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[19]] <- subset(pbmc_seurat_filt[[19]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[19]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[19]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[19]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[19]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[19]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[19]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_10
pbmc_seurat_filt[[20]]
head(pbmc_seurat_filt[[20]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[20]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[20]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[20]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[20]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[20]] <- subset(pbmc_seurat_filt[[20]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[20]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[20]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[20]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[20]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[20]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[20]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_11
pbmc_seurat_filt[[21]]
head(pbmc_seurat_filt[[21]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[21]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[21]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[21]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[21]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[21]] <- subset(pbmc_seurat_filt[[21]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[21]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[21]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[21]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[21]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[21]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[21]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_11
pbmc_seurat_filt[[22]]
head(pbmc_seurat_filt[[22]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[22]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[22]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[22]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[22]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[22]] <- subset(pbmc_seurat_filt[[22]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[22]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[22]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[22]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[22]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[22]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[22]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_NT_12
pbmc_seurat_filt[[23]]
head(pbmc_seurat_filt[[23]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[23]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[23]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[23]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[23]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[23]] <- subset(pbmc_seurat_filt[[23]], subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[23]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[23]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[23]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[23]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[23]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[23]])})

#input, CreateSeuratObject, filter, save
## PBMC_scRNAseq_LL_12
pbmc_seurat_filt[[24]]
head(pbmc_seurat_filt[[24]]@meta.data)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[24]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[24]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[24]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[24]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))

#after filteration
options(repr.plot.width=25,repr.plot.height=15)
pbmc_seurat_filt[[24]] <- subset(pbmc_seurat_filt[[24]], subset =nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

options(repr.plot.width=25,repr.plot.height=15)
VlnPlot(object = pbmc_seurat_filt[[24]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=25,repr.plot.height=10)
cor_plot_umi_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[24]], feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans = 'log2')
cor_plot_genes_mito <- 
  FeatureScatter(object = pbmc_seurat_filt[[24]], feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
cor_plot_umi_genes <- 
  FeatureScatter(object = pbmc_seurat_filt[[24]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")+ theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(
    trans = "log2")

print(plot_grid(cor_plot_umi_mito))
print(plot_grid(cor_plot_genes_mito))
print(plot_grid(cor_plot_umi_genes))


message("Number of cells:", {ncol(pbmc_seurat_filt[[24]])})

message("Number of genes:", {nrow(pbmc_seurat_filt[[24]])})

pbmc_seurat_filt_obj <- Reduce(function(x, y) merge(x, y, do.normalize = F), pbmc_seurat_filt)

pbmc_seurat_filt_obj
pbmc_obj<-pbmc_seurat_filt_obj
pbmc_obj

### Number of cells and genes
message("Number of cells:", {ncol(pbmc_obj)})

message("Number of genes:", {nrow(pbmc_obj)})

message("filtered mean num genes:", {round(mean(pbmc_obj$nFeature_RNA), 3)})
message("filtered median num genes:", {median(pbmc_obj$nFeature_RNA)})
message("filtered mean num UMIs:", {round(mean(pbmc_obj$nCount_RNA), 3)})
message("filtered median num UMIs:", {median(pbmc_obj$nCount_RNA)})
message("filtered mean percent mito:", {round(mean(pbmc_obj$percent.mt), 3)})
message("filtered median percent mito:", {round(median(pbmc_obj$percent.mt),3)})

save(pbmc_obj,file="pbmc_obj.rda")

load("pbmc_obj.rda")

pbmc_obj_list<-SplitObject(pbmc_obj, split.by = "orig.ident")

save(pbmc_obj_list,file="pbmc_obj_list.rda")
