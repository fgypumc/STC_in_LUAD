library(CellChat)
library(tidyverse)
library(Seurat)
library(future)
################
dir='YourPath'
dir.create(paste0(dir,'/9_CellChat'))
setwd(paste0(dir,'/9_CellChat'))

############preparing files
sce=readRDS(paste0(dir,'/3.metacell/metacell_seurat_ec.rds'))
sce_meta=sce@meta.data
all_sce=readRDS(paste0(dir,'/1.load_data/GSE131907_tumor.rds'))
all_sce_meta=all_sce@meta.data
non_tumor=all_sce_meta %>% filter(!Cell_type%in%c('Oligodendrocytes','Undetermined','Epithelial cells'))
table(non_tumor$Cell_type)
Maligant=all_sce_meta[all_sce_meta$Cell_type=='Epithelial cells',]
Proliferative=rownames(Maligant)[!rownames(Maligant)%in%rownames(sce_meta)]
CellChat_cells=subset(all_sce,cells=c(rownames(sce_meta),Proliferative,rownames(non_tumor)))
a=data.frame(cell_type=CellChat_cells$Cell_type) 
a$cell_type=ifelse(rownames(a)%in%rownames(sce_meta),sce_meta$cluster0.25,a$cell_type)
a$cell_type=ifelse(rownames(a)%in%Proliferative,'Maligant cells',a$cell_type)
CellChat_cells=AddMetaData(CellChat_cells,a)
Seurat::Idents(CellChat_cells)=CellChat_cells$cell_type
table(Idents(CellChat_cells))
table(is.na(CellChat_cells$Cell_type))
saveRDS(CellChat_cells,file = paste0(dir,'/9_CellChat/CellChat_cells_allcluster.rds'))

#####################running cellchat
options(stringsAsFactors = FALSE)
cellchat <- createCellChat(object = CellChat_cells, group.by = "cell_type")
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 10) # do parallel
options(future.globals.maxSize= 891289600)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
save(cellchat,file=paste0(dir,'/9_CellChat/CellChat_result_signaling_all_cluster.Rdata'))

