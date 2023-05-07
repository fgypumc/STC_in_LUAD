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
##########
dir='YourPath'
dir.create(paste0(dir,'/6.CS_score'))
setwd(paste0(dir,'/6.CS_score'))

#######extracting genes constituting CS_score
WGCNA_gene=read.csv(paste0(dir,'/5.WGCNA/WGCNAgene.csv'))
mod10=WGCNA_gene[WGCNA_gene$module==10,]
mod12=WGCNA_gene[WGCNA_gene$module==12,]
sen=list()
sen$sen_down=mod12$`Gene symbol`
sen$sen_up=mod10$`Gene symbol`
save(sen,file = paste0(dir,'/6.CS_score/sen_signature.Rdata') )

########calculating CS_score for each cell
load(paste0(dir,'/6.CS_score/sen_signature.Rdata'))
sce=readRDS(paste0(dir,'/3.metacell/metacell_seurat_ec.rds'))
expr=sce@assays$RNA@data
cells_rankings <- AUCell_buildRankings(expr, nCores=30, plotStats=F)
cells_rankings
save(cells_rankings,file=paste0(dir,'/6.CS_score/cells_rankings.Rdata'))
sen_score <- AUCell_calcAUC(sen, cells_rankings)
save(sen_score,file=paste0(dir,'/6.CS_score/sen_score.Rdata'))
sen_score1=t(getAUC(sen_score)) 
sen_score1=as.data.frame(sen_score1)
sen_score1$sen_combined=sen_score1$sen_up - sen_score1$sen_down
sce=AddMetaData(sce,sen_score1)
saveRDS(sce,paste0(dir,'/3.metacell/metacell_seurat_ec.rds'))

################label tumors cells as CS cells and non-CS cells
sce=readRDS(paste0(dir,'/3.metacell/metacell_seurat_ec.rds'))
sen_combined=sce@meta.data[,c(24,25)]
sen_combined_high=head(arrange(sen_combined,desc(sen_combined_score)),1250) 
sen_combined_low=tail(arrange(sen_combined,desc(sen_combined_score)),1250)
meta=sce@meta.data
meta$sen_combined=ifelse(rownames(meta)%in%rownames(sen_combined_high),'high','intermidiate')
meta$sen_combined=ifelse(rownames(meta)%in%rownames(sen_combined_low),'low',meta$sen_combined)
sce=AddMetaData(sce,meta)
table(sce$sen_combined)
table(sce$Sample,sce$sen_combined)
saveRDS(sce,paste0(dir,'/3.metacell/metacell_seurat_ec.rds'))

####################DEGs
Idents(sce) <- sce$sen_combined
sen_combined_markers1 <- FindMarkers(sce, ident.1 = "high", ident.2 = "low", verbose = T, logfc.threshold = 0.25,only.pos =F)
save(sen_combined_markers1,file=paste0(dir,'/6.CS_score/sen_combined_markers_FC.Rdata'))

####################GSEA
sen_combined_markers1$gene<-rownames(sen_combined_markers1)
gsea_genes<-sen_combined_markers1 %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene,avg_log2FC)
ranks <- deframe(gsea_genes)
load(paste0(dir,'/6.CS_score/genesets_all.Rdata'))
fgseaRes <- fgsea(pathways = pathways,
                  stats = ranks ,
                  minSize=5,
                  maxSize=500,
                  nperm=10000)
fgseaRes_marker1=fgseaRes[fgseaRes$pval<0.05,]
save(fgseaRes,file=paste0(dir,'/6.CS_score/fgseaRes.Rdata'))
