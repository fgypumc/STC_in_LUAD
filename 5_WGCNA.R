###############
library(Seurat)
library(devtools)
library(WGCNA)
library(Seurat)
library(tidyverse)
library(reshape2)
library(stringr)
library(readxl)
###############
dir='YourPath'
dir.create(paste0(dir,'/5.WGCNA'))
setwd(paste0(dir,'/5.WGCNA'))

#############preparing the matrix for WGCNA
load(paste0(dir,'/metacell/AllPatient/metacell_exp.Rdata'))
sig=read_excel(paste0(dir,'/5.WGCNA/senescence_signature.xlsx'))
df=meta_df[,colnames(meta_df)%in%sig$`Gene symbol`]
head(df)[1:6,1:6]
df=t(df)
saveRDS(df,paste0(dir,'/5.WGCNA/WGCNA_df.rds'))

load(file=paste0(dir,'/metacell/AllPatient/metacell_GIS.Rdata'))
phenotype=phenotype[colnames(df),]
table(rownames(phenotype)==colnames(df))####TRUE

############ running WGCNA
############checking sample
datExpr_filted = t(df)
gsg = goodSamplesGenes(datExpr_filted, verbose = 3)
gsg$allOK

##########picking up the power
powers = c(c(1:10,by = 1), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr_filted, powerVector = powers, verbose = 5)
power=sft$powerEstimate  ######I use power = 4

############constructing the network
nGenes = ncol(datExpr_filted)
nSamples = nrow(datExpr_filted)
corType = "pearson"
maxPOutliers = ifelse(corType=="pearson",1,0.05) 

cor <- WGCNA::cor
net = blockwiseModules(datExpr_filted, power = power, maxBlockSize = nGenes,#nGenes
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("dataExpr", ".tom"),
                       verbose = 3)

cor<-stats::cor
mergedColors = labels2colors(net$colors)
pdf(file=paste0(dir,'/5.WGCNA/WGCNA.pdf'),height = 8,width = 8)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##################genes in each module
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors) 
color=data.frame(moduleLabels)
color$gene=rownames(color)
write_csv(color,paste0(dir,'/5.WGCNA/WGCNAgene.csv'))

##################the intercorrelation of modules
MEs = net$MEs
MEs_col = MEs
library(stringr)
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf(file=paste0(dir,'/5.WGCNA/intercorrelation.pdf'),height = 8,width = 8)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

###########selecting the most correlated module with GIS
phenotype$GIS_group=ifelse(phenotype$GIS_mean>median(phenotype$GIS_mean),'high','low')
datTraits=phenotype
datTraits <- as.data.frame(datTraits)
# colnames(datTraits) = "condition"
design = model.matrix(~0+ datTraits$GIS_group)
c = as.factor(datTraits$GIS_group)
levels(c)
colnames(design) = levels(c)

moduleColors <- labels2colors(net$colors)
#moduleColors <- net$colors
MEs0 = moduleEigengenes(datExpr_filted, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
nSamples = nrow(datExpr_filted)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot (correlation between module and conditions)
pdf(file=paste0(dir,'/5.WGCNA/GIS_group_color.pdf'),height = 8,width = 8)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
