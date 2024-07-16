### J.G.B.D

setwd("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/ICI/")

library(devtools)
library(dplyr)
library(remotes)
library(Seurat)
library(SeuratWrappers)
library(presto)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(sctransform)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(Matrix)
library(cowplot)
library(biomaRt)
library(pheatmap)
library(celldex)
library(scRNAseq)
library(CVRCFunc)
library(SingleR)
library(dittoSeq)
library(metap)
library(mclust)
library(RColorBrewer)
library(scales)
library(viridisLite)
library(CellChat)
library(slingshot)
library(tradeSeq)
library(Cairo)
library(VennDiagram)


library(devtools)
library(dplyr)
library(remotes)
library(Seurat)
library(SeuratWrappers)
library(presto)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(sctransform)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(Matrix)
library(cowplot)
library(biomaRt)
library(pheatmap)
library(celldex)
library(scRNAseq)
library(CVRCFunc)
library(SingleR)
library(dittoSeq)
library(metap)
library(mclust)
library(RColorBrewer)
library(scales)
library(viridisLite)
library(CellChat)
library(slingshot)
library(uwot)
library(tradeSeq)
library(Cairo)
library(VennDiagram)
library(speckle)

# Everything but the halted progression and baseline samples were removed from downstream analysis.

##### GSE253555 data
mg.m.bl <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Gourvest_2022/tenx_files/male_BL")
mg.m.bl <- CreateSeuratObject(counts = mg.m.bl$`Gene Expression`, project = "mg.m.bl", min.cells = 3, min.features = 200)

mg.f.bl <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Gourvest_2022/tenx_files/female_BL")
mg.f.bl <- CreateSeuratObject(counts = mg.f.bl$`Gene Expression`, project = "mg.f.bl", min.cells = 3, min.features = 200)


##### GSE246316 data
mg.f.hp <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Gourvest_2022/tenx_files/female_HP")
mg.f.hp <- CreateSeuratObject(counts = mg.f.hp$`Gene Expression`, project = "mg.f.hp", min.cells = 3, min.features = 200)

mg.m.hp <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Gourvest_2022/tenx_files/male_HP")
mg.m.hp <- CreateSeuratObject(counts = mg.m.hp$`Gene Expression`, project = "mg.m.hp", min.cells = 3, min.features = 200)

##### unpublished data
mg.f.reg <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Gourvest_2022/tenx_files/female_Reg")
mg.f.reg <- CreateSeuratObject(counts = mg.f.reg$`Gene Expression`, project = "mg.f.reg", min.cells = 3, min.features = 200)

mg.m.reg <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Gourvest_2022/tenx_files/male_Reg")
mg.m.reg <- CreateSeuratObject(counts = mg.m.reg$`Gene Expression`, project = "mg.m.reg", min.cells = 3, min.features = 200)

mg.combined <- merge(mg.f.bl, y = c(mg.f.hp, mg.f.reg, mg.m.bl, mg.m.hp, mg.m.reg))
mg.combined


##### GSE161494 data
ma.ctl <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Afonso_2021/tenx_files/count_C")
ma.ctl <- CreateSeuratObject(counts = ma.ctl, project = "ma.ctl", min.cells = 3, min.features = 200)

ma.mir33 <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Afonso_2021/tenx_files/count_A")
ma.mir33 <- CreateSeuratObject(counts = ma.mir33, project = "ma.mir33", min.cells = 3, min.features = 200)

ma.combined <- merge(ma.ctl, y = ma.mir33)
ma.combined


##### GSE168389 data
msc.wt.bl <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Schlegel_2021/tenx_files/WT_BL")
msc.wt.bl <- CreateSeuratObject(counts = msc.wt.bl, project = "msc.wt.bl", min.cells = 3, min.features = 200)

msc.wt.reg <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Schlegel_2021/tenx_files/WT_Reg")
msc.wt.reg <- CreateSeuratObject(counts = msc.wt.reg, project = "msc.wt.reg", min.cells = 3, min.features = 200)

msc.ko.bl <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Schlegel_2021/tenx_files/Ntn1_KO_BL")
msc.ko.bl <- CreateSeuratObject(counts = msc.ko.bl, project = "msc.ko.bl", min.cells = 3, min.features = 200)

msc.ko.reg <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Schlegel_2021/tenx_files/Ntn1_KO_Reg")
msc.ko.reg <- CreateSeuratObject(counts = msc.ko.reg, project = "msc.ko.reg", min.cells = 3, min.features = 200)

msc.combined <- merge(msc.wt.bl, y = c(msc.wt.reg, msc.ko.bl, msc.ko.reg))
msc.combined


##### GSE141038 data
msh.ctl <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Sharma_2020/tenx_files/Ctl_10X")
msh.ctl <- CreateSeuratObject(counts = msh.ctl, project = "msh.ctl", min.cells = 3, min.features = 200)

msh.bl <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Sharma_2020/tenx_files/BL_10X")
msh.bl <- CreateSeuratObject(counts = msh.bl, project = "msh.bl", min.cells = 3, min.features = 200)

msh.cd25 <- Read10X("/Volumes/Karsten/MooreLab/Sequencing/scRNAseq/Plaque/Mouse/Sharma_2020/tenx_files/Treat_10X")
msh.cd25 <- CreateSeuratObject(counts = msh.cd25, project = "msh.cd25", min.cells = 3, min.features = 200)

msh.combined <- merge(msh.ctl, y = c(msh.bl, msh.cd25))
msh.combined

#Add metadata


mg.df <- data.frame(mg.combined@meta.data)
mg.df <- mg.df %>% mutate(sex = case_when(startsWith(mg.df$orig.ident, "mg.f") ~ "Female", startsWith(mg.df$orig.ident, "mg.m") ~ "Male"))
mg.df <- mg.df %>% mutate(endpoint = case_when(endsWith(mg.df$orig.ident, "bl") ~ "Baseline", endsWith(mg.df$orig.ident, "hp") ~ "Chow", endsWith(mg.df$orig.ident, "reg") ~ "Chow + ApoB ASO"))
mg.df <- mg.df %>% mutate(athero = case_when(endsWith(mg.df$orig.ident, "bl") ~ "Baseline", endsWith(mg.df$orig.ident, "hp") ~ "Halted Progression", endsWith(mg.df$orig.ident, "reg") ~ "Regression"))
mg.df <- mg.df %>% mutate(background = case_when(startsWith(mg.df$orig.ident, "mg") ~ "LdlrKO")) 
mg.df <- mg.df %>% mutate(treatment = case_when(endsWith(mg.df$orig.ident, "bl") ~ "None", endsWith(mg.df$orig.ident, "hp") ~ "None", endsWith(mg.df$orig.ident, "reg") ~ "ApoB ASO"))
mg.df <- mg.df %>% mutate(investigator = case_when(startsWith(mg.df$orig.ident, "mg") ~ "Gourvest"))

mg.combined <- AddMetaData(object = mg.combined, metadata = mg.df$sex, col.name = 'Sex')
mg.combined <- AddMetaData(object = mg.combined, metadata = mg.df$endpoint, col.name = 'Endpoint')
mg.combined <- AddMetaData(object = mg.combined, metadata = mg.df$athero, col.name = 'Athero')
mg.combined <- AddMetaData(object = mg.combined, metadata = mg.df$background, col.name = 'Background')
mg.combined <- AddMetaData(object = mg.combined, metadata = mg.df$treatment, col.name = 'Treatment')
mg.combined <- AddMetaData(object = mg.combined, metadata = mg.df$investigator, col.name = 'Investigator')

table(mg.combined@meta.data$Sex)
table(mg.combined@meta.data$Endpoint)
table(mg.combined@meta.data$Athero)
table(mg.combined@meta.data$Background)
table(mg.combined@meta.data$Treatment)
table(mg.combined@meta.data$Investigator)


ma.df <- data.frame(ma.combined@meta.data)
ma.df <- ma.df %>% mutate(sex = case_when(startsWith(ma.df$orig.ident, "ma") ~ "Male"))
ma.df <- ma.df %>% mutate(endpoint = case_when(endsWith(ma.df$orig.ident, "ctl") ~ "Chow", endsWith(ma.df$orig.ident, "mir33") ~ "Chow + miR33 ASO"))
ma.df <- ma.df %>% mutate(athero = case_when(endsWith(ma.df$orig.ident, "ctl") ~ "Halted Progression", endsWith(ma.df$orig.ident, "mir33") ~ "Regression"))
ma.df <- ma.df %>% mutate(background = case_when(startsWith(ma.df$orig.ident, "ma") ~ "LdlrKO"))
ma.df <- ma.df %>% mutate(treatment = case_when(endsWith(ma.df$orig.ident, "mir33") ~ "miR33 ASO", endsWith(ma.df$orig.ident, "ctl") ~ "Control ASO"))
ma.df <- ma.df %>% mutate(investigator = case_when(startsWith(ma.df$orig.ident, "ma") ~ "Afonso"))

ma.combined <- AddMetaData(object = ma.combined, metadata = ma.df$sex, col.name = 'Sex')
ma.combined <- AddMetaData(object = ma.combined, metadata = ma.df$endpoint, col.name = 'Endpoint')
ma.combined <- AddMetaData(object = ma.combined, metadata = ma.df$athero, col.name = 'Athero')
ma.combined <- AddMetaData(object = ma.combined, metadata = ma.df$background, col.name = 'Background')
ma.combined <- AddMetaData(object = ma.combined, metadata = ma.df$treatment, col.name = 'Treatment')
ma.combined <- AddMetaData(object = ma.combined, metadata = ma.df$investigator, col.name = 'Investigator')

table(ma.combined@meta.data$Sex)
table(ma.combined@meta.data$Endpoint)
table(ma.combined@meta.data$Athero)
table(ma.combined@meta.data$Background)
table(ma.combined@meta.data$Treatment)
table(ma.combined@meta.data$Investigator)


msc.df <- data.frame(msc.combined@meta.data)
msc.df <- msc.df %>% mutate(sex = case_when(startsWith(msc.df$orig.ident, "msc") ~ "Male"))
msc.df <- msc.df %>% mutate(endpoint = case_when(endsWith(msc.df$orig.ident, "bl") ~ "Baseline", endsWith(msc.df$orig.ident, "ko.reg") ~ "Chow + Ntn1 KO", endsWith(msc.df$orig.ident, "wt.reg") ~ "Chow"))
msc.df <- msc.df %>% mutate(athero = case_when(endsWith(msc.df$orig.ident, "bl") ~ "Baseline", endsWith(msc.df$orig.ident, "ko.reg") ~ "HP + Ntn1 KO", endsWith(msc.df$orig.ident, "wt.reg") ~ "Halted Progression"))
msc.df <- msc.df %>% mutate(background = case_when(endsWith(msc.df$orig.ident, "bl") ~ "WT", startsWith(msc.df$orig.ident, "msc.wt.reg") ~ "WT", startsWith(msc.df$orig.ident, "msc.ko.reg") ~ "Ntn1KO"))
msc.df <- msc.df %>% mutate(treatment = case_when(endsWith(msc.df$orig.ident, "bl") ~ "None", endsWith(msc.df$orig.ident, "reg") ~ "Chow + tamoxifen"))
msc.df <- msc.df %>% mutate(investigator = case_when(startsWith(msc.df$orig.ident, "msc") ~ "Schlegel"))

msc.combined <- AddMetaData(object = msc.combined, metadata = msc.df$sex, col.name = 'Sex')
msc.combined <- AddMetaData(object = msc.combined, metadata = msc.df$endpoint, col.name = 'Endpoint')
msc.combined <- AddMetaData(object = msc.combined, metadata = msc.df$athero, col.name = 'Athero')
msc.combined <- AddMetaData(object = msc.combined, metadata = msc.df$background, col.name = 'Background')
msc.combined <- AddMetaData(object = msc.combined, metadata = msc.df$treatment, col.name = 'Treatment')
msc.combined <- AddMetaData(object = msc.combined, metadata = msc.df$investigator, col.name = 'Investigator')

table(msc.combined@meta.data$Sex)
table(msc.combined@meta.data$Endpoint)
table(msc.combined@meta.data$Athero)
table(msc.combined@meta.data$Background)
table(msc.combined@meta.data$Treatment)
table(msc.combined@meta.data$Investigator)


msh.df <- data.frame(msh.combined@meta.data)
msh.df <- msh.df %>% mutate(sex = case_when(startsWith(msh.df$orig.ident, "msh") ~ "M&F"))
msh.df <- msh.df %>% mutate(endpoint = case_when(endsWith(msh.df$orig.ident, "bl") ~ "Baseline", endsWith(msh.df$orig.ident, "ctl") ~ "Chow + ApoB ASO + IgG", endsWith(msh.df$orig.ident, "cd25") ~ "Chow + ApoB ASO + αCd25"))
msh.df <- msh.df %>% mutate(athero = case_when(endsWith(msh.df$orig.ident, "bl") ~ "Baseline", endsWith(msh.df$orig.ident, "ctl") ~ "Regression", endsWith(msh.df$orig.ident, "cd25") ~ "Regression + αCd25"))
msh.df <- msh.df %>% mutate(treatment = case_when(endsWith(msh.df$orig.ident, "cd25") ~ "ApoB ASO + αCD25", endsWith(msh.df$orig.ident, "ctl") ~ "ApoB ASO + IgG", endsWith(msh.df$orig.ident, "bl") ~ "None"))
msh.df <- msh.df %>% mutate(background = case_when(startsWith(msh.df$orig.ident, "msh") ~ "WT"))
msh.df <- msh.df %>% mutate(investigator = case_when(startsWith(msh.df$orig.ident, "msh") ~ "Sharma"))

msh.combined <- AddMetaData(object = msh.combined, metadata = msh.df$sex, col.name = 'Sex')
msh.combined <- AddMetaData(object = msh.combined, metadata = msh.df$endpoint, col.name = 'Endpoint')
msh.combined <- AddMetaData(object = msh.combined, metadata = msh.df$athero, col.name = 'Athero')
msh.combined <- AddMetaData(object = msh.combined, metadata = msh.df$background, col.name = 'Background')
msh.combined <- AddMetaData(object = msh.combined, metadata = msh.df$treatment, col.name = 'Treatment')
msh.combined <- AddMetaData(object = msh.combined, metadata = msh.df$investigator, col.name = 'Investigator')

#Deal with mitochondrial RNA

mg.combined[["percent.mt"]] <- PercentageFeatureSet(mg.combined, pattern = "^mt-")

VlnPlot(mg.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)
VlnPlot(object = mg.combined, features = "percent.mt",pt.size = 0.01) + ggtitle("Mitochondrial Content Per Cell") + scale_y_continuous(breaks=seq(0,100,5))

plot1 <- FeatureScatter(mg.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mg.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Trim mitochondrial genes based on plots
mg.combined <- subset(mg.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 7.5)


ma.combined[["percent.mt"]] <- PercentageFeatureSet(ma.combined, pattern = "^mt-")

VlnPlot(ma.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)
VlnPlot(object = ma.combined, features = "percent.mt",pt.size = 0.01) + ggtitle("Mitochondrial Content Per Cell") + scale_y_continuous(breaks=seq(0,100,5))

plot1 <- FeatureScatter(ma.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ma.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Trim mitochondrial genes based on plots
ma.combined <- subset(ma.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 7)


msc.combined[["percent.mt"]] <- PercentageFeatureSet(msc.combined, pattern = "^mt-")
VlnPlot(msc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)
VlnPlot(object = msc.combined, features = "percent.mt",pt.size = 0.01) + ggtitle("Mitochondrial Content Per Cell") + scale_y_continuous(breaks=seq(0,100,5))

plot1 <- FeatureScatter(msc.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(msc.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Trim mitochondrial genes based on plots
msc.combined <- subset(msc.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 6)


msh.combined[["percent.mt"]] <- PercentageFeatureSet(msh.combined, pattern = "^mt-")

plot1 <- FeatureScatter(msh.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(msh.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Trim mitochondrial genes based on plots
msh.combined <- subset(msh.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 8)

###

#Normalization
mg.combined <- NormalizeData(mg.combined, normalization.method = "LogNormalize", scale.factor = 10000)
ma.combined <- NormalizeData(ma.combined, normalization.method = "LogNormalize", scale.factor = 10000)
msh.combined <- NormalizeData(msh.combined, normalization.method = "LogNormalize", scale.factor = 10000)
msc.combined <- NormalizeData(msc.combined, normalization.method = "LogNormalize", scale.factor = 10000)

data.combined <- merge(mg.combined, y = c(ma.combined, msh.combined, msc.combined))


#Cell cycle
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

#computed in Cellcycle_genes.R
cc.genes.15 <- read.csv("putative_cc_genes_mouse.csv", header = F)
cc.genes.15 <- cc.genes.15$Gene

#INTEGRATION

data.combined[["RNAsub"]] <- CreateAssayObject(counts = data.combined@assays$RNA@counts[which(!row.names(data.combined@assays$RNA@counts) %in% cc.genes.15),])

data.combined.list <- SplitObject(data.combined, split.by = "orig.ident")

for (i in 1:length(data.combined.list)) {
  data.combined.list[[i]] <- SCTransform(data.combined.list[[i]],  assay = 'RNA', new.assay.name = 'SCT') 
}

for (i in 1:length(data.combined.list)) {
  data.combined.list[[i]] <- CellCycleScoring(data.combined.list[[i]], assay = 'SCT', s.features = m.s.genes, g2m.features = m.g2m.genes) 
}

for (i in 1:length(data.combined.list)) {
  data.combined.list[[i]] <- SCTransform(data.combined.list[[i]],  assay = 'RNAsub', new.assay.name = 'SCTsub') 
}

features <- SelectIntegrationFeatures(data.combined.list, nfeatures = 3000)
data.combined.list <- PrepSCTIntegration(data.combined.list, anchor.features = features)
data.combined.list <- lapply(X = data.combined.list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

all.genes <- rownames(data.combined)

anchors  <- FindIntegrationAnchors(data.combined.list, normalization.method = "SCT", anchor.features = all.genes, dims = 1:30, reduction = "rpca")
data.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)

data.integrated <- RunPCA(data.integrated, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, dims = 1:30)

data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
data.integrated <- FindClusters(data.integrated, resolution = 0.55)