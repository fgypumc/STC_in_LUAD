####################
library(data.table)
library(useful)
library(Seurat)
library(dplyr)
library(Matrix)
library(R.utils)
###################
dir='YourPath'
dir.create(paste0(dir,'/1.load_data'))
setwd(paste0(dir,'/1.load_data'))

############
gzfile ='GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz'
gunzip(gzfile, remove = F)
rawdata=readRDS(paste(dir,"/1.load_data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds", sep=""))
rawdata<- CreateSeuratObject(counts = rawdata, project = "LUAD", min.cells = 3, min.features = 200)
head(rawdata)[1:6,1:6]
metadata <- read.table(paste(dir,"/1.load_data/GSE131907_Lung_Cancer_cell_annotation (1).txt", sep=""),sep = '\t',header = T)
rownames(metadata)=metadata$Index
rawdata <- AddMetaData(object = rawdata, metadata = metadata)
head(rawdata@meta.data)
saveRDS(rawdata,'GSE131907.rds')

#######selecting tumor samples
rawdata=readRDS('GSE131907.rds')
metadata=rawdata@meta.data
table(metadata$Sample_Origin)
metadata$work=ifelse(metadata$Sample_Origin%in%c('tL/B','mBrain','tLung','mLN'),1,0)
table(metadata$work)
rawdata$work=metadata$work
tumor=SplitObject(rawdata,'work')
tumor=tumor[[2]]
saveRDS(tumor,'GSE131907_tumor.rds')

###############quality check and clustering
sce=readRDS(paste0(dir,'/1.load_data/GSE131907_tumor.rds'))
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
sce<- subset(x=sce, subset = nFeature_RNA > 200 & percent.mt < 20)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(object = sce)
sce <- RunPCA(object = sce,features = VariableFeatures(object = sce))
sce <- JackStraw(sce, num.replicate = 100)
sce <- ScoreJackStraw(sce, dims = 1:20)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.7)
sce <- RunUMAP(sce, dims = 1:10)
saveRDS(sce,paste0(dir,'/1.load_data/GSE131907_tumor.rds'))
