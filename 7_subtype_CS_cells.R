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
dir.create(paste0(dir,'/7.subtype_CS_cells'))
setwd(paste0(dir,'/7.subtype_CS_cells'))

#############clustering
sce=readRDS(paste0(dir,'/3.metacell/metacell_seurat_ec.rds'))
sce=subset(sce,subset=sen_combined=='high')
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(object = sce)
sce <- RunPCA(object = sce,features = VariableFeatures(object = sce))
sce <- JackStraw(sce, num.replicate = 100)
sce <- ScoreJackStraw(sce, dims = 1:20)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:5)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:5)
saveRDS(sce,paste0(dir,'/7.subtype_CS_cells/CS_cells.rds'))

##########################DEGs
cluster_marker <- FindAllMarkers(sce, verbose = T, logfc.threshold = 0.25,only.pos =F)
save(cluster_marker ,file=paste0(dir,'/7.subtype_CS_cells/cluster_marker.Rdata'))

##########################fgsea
outTab=data.frame()
for (i in unique(cluster_marker$cluster)) {
  marker=cluster_marker[cluster_marker$cluster==i,]
  gsea_genes<- marker %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene,avg_log2FC)
  ranks <- deframe(gsea_genes)
  load(paste0(dir,'/5.GSVA/genesets_all.Rdata'))
  fgseaRes <- fgsea(pathways = pathways,
                    stats = ranks ,
                    minSize=5,
                    maxSize=500,
                    nperm=10000)
  fgseaRes_marker1=fgseaRes[fgseaRes$pval<0.05,]
  fgseaRes_marker1$cluster=i
  outTab=rbind(outTab,fgseaRes_marker1)
}
save(outTab,file=paste0(dir,'/7.subtype_CS_cells/fgsea.Rdata'))

##########################SASP factors
load(paste0(dir,'/6.CS_score/sen_signature.Rdata'))#27个基因
SASP1=read_excel(paste0(dir,'/6.CS_score/SASP_factors.xlsx'))###22个基因
marker=c('TGFB1','GDF15','VEGFA','ICAM1','IL6','CXCL3','CXCL2','CXCL1','CCL2','CXCL8')
sce=readRDS(paste0(dir,'/7.subtype_CS_cells/CS_cells.rds'))
table(sce$seurat_clusters)
Idents(sce)=sce$seurat_clusters
pdf(file=paste0(dir,'/7.subtype_CS_cells/SASP_factor_bubble.pdf'),width = 8,height = 5)
DotPlot(sce, features = marker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()

##########################APOPTOSIS factors
load(paste0(dir,'/6.CS_score/sen_signature.Rdata'))
anti_apo=sen_all$`ANTI- APOPTOSIS`
pro_apo=sen_all$`PRO-APOPTOSIS` 
apo=unique(c(sen_all$`ANTI- APOPTOSIS`,sen_all$`PRO-APOPTOSIS`))
Idents(sce)=sce$cluster0.25
pdf(file=paste0(dir,'/7.subtype_CS_cells/APOPTOSIS_bubble.pdf'),width = 8,height = 5)
DotPlot(sce, features = apo)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()
