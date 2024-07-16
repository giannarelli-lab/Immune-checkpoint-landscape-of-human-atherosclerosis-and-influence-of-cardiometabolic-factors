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

Whole_dataset_annotated_human.rds <- readRDS(file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/CHORD_PBMC/pbmc_obj_rm_doub_integrated_annotated.rds')
Whole_dataset_annotated_human.rds <- NormalizeData(Whole_dataset_annotated_human.rds, assay = 'RNA')
Idents(Whole_dataset_annotated_human.rds) <- 'Diagnosis'
Whole_dataset_annotated_human.rds_sub <- subset(Whole_dataset_annotated_human.rds, idents = 'CONTROL')
Idents(Whole_dataset_annotated_human.rds_sub) <- 'Condition'
Whole_dataset_annotated_human.rds_sub <- subset(Whole_dataset_annotated_human.rds_sub, idents = 'NT')

cellchat <- createCellChat(object = Whole_dataset_annotated_human.rds_sub@assays$RNA@data, meta = Whole_dataset_annotated_human.rds_sub@meta.data, group.by = 'annotation_fine')
cellchat <- addMeta(cellchat, meta = Whole_dataset_annotated_human.rds_sub@meta.data)
cellchat <- setIdent(cellchat, ident.use = "annotation_fine") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
interaction_input <- read.csv(file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/cellchat/interaction_input_CellChatDB_ICI_editICname.csv', row.names = 1)
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
future::plan("multiprocess", workers = 39) # do parallel (doesnt work via Rstudio)
cellchat <- identifyOverExpressedGenes(cellchat) # ADD CUSTOM GENES VIA PSEUDOBULKING DE ANALYSIS
cellchat <- identifyOverExpressedInteractions(cellchat)
future::plan("multiprocess", workers = 1)
cellchat <- computeCommunProb(cellchat, nboot = 100, type = "truncatedMean", trim = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, file = '/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/CHORD_PBMC/Cellchat/cellchat_fine_1_NT_CONTROL.rds')
