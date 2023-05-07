library(SingleCellExperiment)
library(ExperimentHub)
library(genomicInstability)
library(readr)
library(dplyr)
library(mixtools)
library(Seurat)
library(clusterProfiler)
############
dir='YourPath'
dir.create(paste0(dir,'/4.genomicInstability'))
setwd(paste0(dir,'/4.genomicInstability'))
############
sce=readRDS(paste0(dir,"/3.metacell/metacell_seurat_ec.rds"))
matrix= as.matrix(GetAssayData(sce,slot = 'data')) 
a=bitr(rownames(matrix), fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db', drop = TRUE)
matrix1=matrix[a$SYMBOL,]
rownames(matrix1)=a$ENTREZID
cnv <- inferCNV(matrix1)
class(cnv)
dim(cnv$nes)
cnv <- genomicInstabilityScore(cnv)
names(cnv)
length(cnv$gis)
par(mai=c(.8, .8, .2, .2))
plot(density(cnv$gis), lwd=2, xlab="GIS", main="")
cnv <- genomicInstabilityScore(cnv)
names(cnv)
length(cnv$gis)
giDensityPlot(cnv)
cnv <- giLikelihood(cnv, recompute=FALSE, normal=1:2, tumor=3)
save(cnv,file = paste0(dir,'/4.genomicInstability/cnv.Rdata'))

load(file = paste0(dir,'/genomicInstability/cnv.Rdata'))
cnv_norm <- inferCNV(matrix1, nullmat=matrix1[, cnv$gi_likelihood<.25, drop=FALSE])
names(cnv_norm)  
cnv_norm <- genomicInstabilityScore(cnv_norm, likelihood=TRUE)
save(cnv_norm,file = paste0(dir,'/genomicInstability/cnv_norm.Rdata'))

############caculating GIS score for metacell
load( paste0(dir,'/3.metacell/testdb/mc.pbmc_mc_f.Rda'))
cluster=data.frame(cluster=object@mc)
load(paste0(dir,'/genomicInstability/cnv_norm.Rdata'))
mat=data.frame(cnv_norm$gis) 
sameID=intersect(rownames(cluster),rownames(mat))
rt=data.frame(cluster=cluster[sameID,],GIS=mat[sameID,])
rownames(rt)=sameID
phenotype<-aggregate(list(rt[,2:length(colnames(rt))]),
                     list(rt[,1]),FUN=mean)
colnames(phenotype)[2]='GIS_mean'
phenotype$Group.1=paste0('C', phenotype$Group.1)
rownames(phenotype)=phenotype$Group.1
phenotype$GIS_group=ifelse(phenotype$GIS_mean>median(phenotype$GIS_mean),'high','low')
save(phenotype,file=paste0(dir,'/4.genomicInstability/metacell_GIS.Rdata'))
