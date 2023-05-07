###############
library(Seurat)
library(devtools)
library(clustree)
library(tidyverse)
library(gridExtra)
library(ggridges)
library(ggplot2)
library(ggExtra)
library(clustree)
library(DoubletFinder)
###############
dir='YourPath'
dir.create(paste0(dir,'/2.inferCNV'))
setwd(paste0(dir,'/2.inferCNV'))

###############selecting epithelial and stromal cells
sce=readRDS(paste0(dir,'/1.load_data/GSE131907_tumor.rds'))

ec <- row.names(sce@meta.data)[which(sce@meta.data$Cell_type.refined %in% c('Epithelial cells','Endothelial cells','Fibroblasts'))]
sce <- subset(sce, cells=ec)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(object = sce)
sce <- RunPCA(object = sce,features = VariableFeatures(object = sce))
saveRDS(sce,paste0(dir,"/2.inferCNV/NonImmune.rds"))

###############preparing files for running inferCNV
NonImmune=readRDS(paste0(dir,"/2.inferCNV/NonImmune.rds"))
table(NonImmune$Cell_type)###"Epithelial cells"  "Endothelial cells" "Fibroblasts"    
name=data.frame(NonImmune$Cell_type)
Endothelial=name[name$NonImmune.Cell_type=="Endothelial cells",,drop=F]
Fibroblasts=name[name$NonImmune.Cell_type=="Fibroblasts",,drop=F]
Epithelial=name[name$NonImmune.Cell_type=="Epithelial cells",,drop=F]
endo1=rownames(Endothelial)[1:300]
endo2=rownames(Endothelial)[301:600]
fibro1=rownames(Fibroblasts)[1:2000]
fibro2=rownames(Fibroblasts)[2001:2600]
sampled=c(endo1,fibro1,endo2,fibro2,rownames(Epithelial))
sampled_sce=subset(NonImmune,cells=sampled)
df=as.matrix(sampled_sce@assays$RNA@counts)
write.table(df,paste0(dir,"/2.inferCNV/df.txt"),sep='\t')

metadata=sampled_sce@meta.data
metadata=metadata[,8,drop=F]
a=metadata[1:300,,drop=F]
b=metadata[301:600,,drop=F]
a$Cell_type='Endothelial_cells_reference'
b$Cell_type='Fibroblasts_reference'
metadata=rbind(a,b, metadata[601:nrow(metadata),,drop=F])
table(metadata$Cell_type)
write.table(metadata,'cellAnnotations.txt',sep = '\t',col.names = F)
table(colnames(df)==rownames(metadata)) 
identical(rownames(dat),rownames(geneInfor))

library(AnnoProbe)   
geneInfor=annoGene(rownames(df),"SYMBOL",'human')    
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
write.table(geneInfor,'genePos.txt',sep = '\t',row.names = F,col.names = F)

############### running inferCNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix='df.txt', # 可以直接提供矩阵对象
                                    annotations_file="cellAnnotations.txt",
                                    delim="\t",
                                    gene_order_file="genePos.txt",
                                    ref_group_names=c('Endothelial_cells_reference','Fibroblasts_reference'))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  
                             cluster_by_groups=T,  
                             denoise=T, 
                             HMM=F) 