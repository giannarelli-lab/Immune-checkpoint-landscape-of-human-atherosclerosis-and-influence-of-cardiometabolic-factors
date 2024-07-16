#Now, on HPC:
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

#CellChatDB <- CellChatDB.mouse
#interaction_input <- CellChatDB$interaction
#complex_input <- CellChatDB$complex
#cofactor_input <- CellChatDB$cofactor
#geneInfo <- CellChatDB$geneInfo
#write.csv(interaction_input, file = "cellchat/interaction_input_CellChatDB.csv")
#write.csv(complex_input, file = "cellchat/complex_input_CellChatDB.csv")
#write.csv(cofactor_input, file = "cellchat/cofactor_input_CellChatDB.csv")
#write.csv(geneInfo, file = "cellchat/geneInfo_input_CellChatDB.csv")
#Now go on the Excel file and add the ICI gene interactions to interaction_input

interaction_input <- read.csv(file = 'interaction_input_CellChatDB_ICI.csv', row.names = 1)
complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo

Whole_dataset_annotated_mouse <- readRDS(file = '../Whole_dataset_annotated_mouse_08152023_BLHP.rds')
Baseline_dataset_annotated_mouse <- subset(x = Whole_dataset_annotated_mouse, subset = Athero == c('Baseline'))
HP_dataset_annotated_mouse <- subset(x = Whole_dataset_annotated_mouse, subset = Athero == c('Halted Progression'))

Whole_dataset_annotated_mouse <- NormalizeData(Whole_dataset_annotated_mouse, assay = 'RNA')
Baseline_dataset_annotated_mouse <- NormalizeData(Baseline_dataset_annotated_mouse, assay = 'RNA')
HP_dataset_annotated_mouse <- NormalizeData(HP_dataset_annotated_mouse, assay = 'RNA')

cellchat <- createCellChat(object = Whole_dataset_annotated_mouse@assays$RNA@data, meta = Whole_dataset_annotated_mouse@meta.data, group.by = 'annotation_fine')
cellchat.bl <- createCellChat(object = Baseline_dataset_annotated_mouse@assays$RNA@data, meta = Baseline_dataset_annotated_mouse@meta.data, group.by = 'annotation_fine')
cellchat.hp <- createCellChat(object = HP_dataset_annotated_mouse@assays$RNA@data, meta = HP_dataset_annotated_mouse@meta.data, group.by = 'annotation_fine')

cellchat <- addMeta(cellchat, meta = Whole_dataset_annotated_mouse@meta.data)
cellchat.bl <- addMeta(cellchat.bl, meta = Baseline_dataset_annotated_mouse@meta.data)
cellchat.hp <- addMeta(cellchat.hp, meta = HP_dataset_annotated_mouse@meta.data)

cellchat <- setIdent(cellchat, ident.use = "annotation_fine") # set "labels" as default cell identity
cellchat.bl <- setIdent(cellchat.bl, ident.use = "annotation_fine") # set "labels" as default cell identity
cellchat.hp <- setIdent(cellchat.hp, ident.use = "annotation_fine") # set "labels" as default cell identity

levels(cellchat@idents) # show factor levels of the cell labels
levels(cellchat.bl@idents) # show factor levels of the cell labels
levels(cellchat.hp@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize <- as.numeric(table(cellchat.bl@idents)) # number of cells in each cell group
groupSize <- as.numeric(table(cellchat.hp@idents)) # number of cells in each cell group

showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat.bl@DB <- CellChatDB.use
cellchat.hp@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat.bl <- subsetData(cellchat.bl) # This step is necessary even if using the whole database
cellchat.hp <- subsetData(cellchat.hp) # This step is necessary even if using the whole database

future::plan("multisession", workers = 30) # do parallel (doesnt work via Rstudio)

cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat.bl <- identifyOverExpressedGenes(cellchat.bl) 
cellchat.bl <- identifyOverExpressedInteractions(cellchat.bl)

cellchat.hp <- identifyOverExpressedGenes(cellchat.hp) 
cellchat.hp <- identifyOverExpressedInteractions(cellchat.hp)

future::plan("multisession", workers = 1)

cellchat <- computeCommunProb(cellchat, nboot = 100, type = "truncatedMean", trim = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat.bl <- computeCommunProb(cellchat.bl, nboot = 100, type = "truncatedMean", trim = 0.01)
cellchat.bl <- filterCommunication(cellchat.bl, min.cells = 10)
cellchat.bl <- computeCommunProbPathway(cellchat.bl)
cellchat.bl <- aggregateNet(cellchat.bl)

cellchat.hp <- computeCommunProb(cellchat.hp, nboot = 100, type = "truncatedMean", trim = 0.01)
cellchat.hp <- filterCommunication(cellchat.hp, min.cells = 10)
cellchat.hp <- computeCommunProbPathway(cellchat.hp)
cellchat.hp <- aggregateNet(cellchat.hp)


saveRDS(cellchat, file = 'cellchat_mouse_fine.rds')
saveRDS(cellchat.bl, file = 'cellchat.bl_mouse_fine.rds')
saveRDS(cellchat.hp, file = 'cellchat.hp_mouse_fine.rds')
