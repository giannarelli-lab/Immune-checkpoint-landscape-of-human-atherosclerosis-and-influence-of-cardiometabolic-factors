library(Seurat)
library(ggplot2)
library(presto)
library(viridis)
library(pheatmap)
library(Hmisc)
library(CVRCFunc)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(CellChat)

Whole_dataset_annotated_human.rds <- readRDS(file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/Whole_dataset_annotated_human.rds')
Whole_dataset_annotated_human.rds <- NormalizeData(Whole_dataset_annotated_human.rds, assay = 'RNA')

cellchat <- createCellChat(object = Whole_dataset_annotated_human.rds@assays$RNA@data, meta = Whole_dataset_annotated_human.rds@meta.data, group.by = 'annotation_fine')
cellchat <- addMeta(cellchat, meta = Whole_dataset_annotated_human.rds@meta.data)
cellchat <- setIdent(cellchat, ident.use = "annotation_fine") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human
interaction_input <- read.csv(file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/cellchat/rerun_ICdb_issue/interaction_input_CellChatDB_ICI_editICname.csv', row.names = 1)
complex_input <- read.csv(file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/cellchat/complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/cellchat/cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/cellchat/geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 30)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
future::plan("multisession", workers = 1)
cellchat <- computeCommunProb(cellchat, nboot = 100, type = "truncatedMean", trim = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/cellchat/rerun_ICdb_issue/cellchat_fine_1.rds')


