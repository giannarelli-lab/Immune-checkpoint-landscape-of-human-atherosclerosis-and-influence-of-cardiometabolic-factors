library(Seurat)

A = Read10X(data.dir = 'E:/HIMC_AWS_FastQ_Data/GSM7502475_Sample38/filtered_feature_bc_matrix')
A[["HTO"]] <- A$Custom[22:23,]
A$Custom <- A$Custom[1:21,]
citeseq <- CreateSeuratObject(A$`Gene Expression`, project = "CITE-Seq")
ADT_assay <- CreateAssayObject(counts = A$Custom)
HTO_assay <- CreateAssayObject(counts = A$HTO)
citeseq[["ADT"]] <- ADT_assay
citeseq[["HTO"]] <- HTO_assay
Assays(citeseq)
DefaultAssay(citeseq)

citeseq <- NormalizeData(citeseq)
citeseq <- FindVariableFeatures(citeseq, selection.method = "mean.var.plot")
citeseq <- ScaleData(citeseq, features = VariableFeatures(citeseq))

citeseq <- NormalizeData(citeseq, assay = "HTO", normalization.method = "CLR")
citeseq <- HTODemux(citeseq, assay = "HTO", positive.quantile = 0.89, kfunc = "clara")
table(citeseq$HTO_classification.global)
Idents(citeseq) <- "HTO_maxID"
RidgePlot(citeseq, assay = "HTO", features = rownames(citeseq[["HTO"]])[1:2], ncol = 2)
FeatureScatter(citeseq, feature1 = "HTO-7", feature2 = "HTO-8")
Idents(citeseq) <- "HTO_classification.global"
VlnPlot(citeseq, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
plot <- FeatureScatter(citeseq, feature1 = "HTO-7", feature2 = "HTO-8")
plot$coordinates$limits$x <- c(0,6)
plot$coordinates$limits$y <- c(0,6)
plot

citeseq.subset <- subset(citeseq, idents = "Negative", invert = TRUE)

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = citeseq.subset, assay = "HTO"))))

citeseq.subset <- RunTSNE(citeseq.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(citeseq.subset)

HTOHeatmap(citeseq, assay = "HTO")

citeseq.singlet <- subset(citeseq, idents = "Singlet")

citeseq.singlet <- FindVariableFeatures(citeseq.singlet, selection.method = "vst")

citeseq.singlet<- ScaleData(citeseq.singlet, features = VariableFeatures(citeseq.singlet))

citeseq.singlet <- RunPCA(citeseq.singlet, features = VariableFeatures(citeseq.singlet))
citeseq.singlet <- FindNeighbors(citeseq.singlet, reduction = "pca", dims = 1:10)
citeseq.singlet <- FindClusters(citeseq.singlet, resolution = 0.6, verbose = FALSE)
citeseq.singlet <- RunUMAP(citeseq.singlet, reduction = "pca", dims = 1:10)

DimPlot(citeseq.singlet, group.by = "HTO_classification")

Idents(citeseq.singlet)<-"hash.ID"
citeseq.plaque <- subset(citeseq.singlet, idents = "HTO-7")
saveRDS(citeseq.plaque, file = "GSM7502475_Sample38_palque_subset.rds")