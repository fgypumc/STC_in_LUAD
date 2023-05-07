library(ggplot2)
library(dplyr)
library(msigdbr)
library(Seurat)
library(pheatmap)
library(patchwork)
library(limma)
library(readxl)
library(AUCell) 
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
##############
dir='YourPath'
dir.create(paste0(dir,'/8.dorothea'))
setwd(paste0(dir,'/8.dorothea'))

##############dorothea
library(future)
plan()
plan("multiprocess", workers = 12)
plan()

################ Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

############ obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

##########running dorothea
sce <- run_viper(sce, regulon,
                 options = list(method = "scale", minsize = 4, 
                                eset.filter = FALSE, cores = 12,  
                                verbose = FALSE))
saveRDS(sce,paste0(dir,'/8.dorothea/CS_cells.rds'))
