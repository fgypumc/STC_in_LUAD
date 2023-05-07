library(dplyr)
library(GSVA)
################
dir='YourPath'
dir.create(paste0(dir,'/11_immune_cell_infiltration'))
setwd(paste0(dir,'/11_immune_cell_infiltration'))
#################caculate immune_cell infiltration
load(paste0(dir,'/11_immune_cell_infiltration/immune.gene.lists.Rdata'))
GSVA=function(x){
  gsva(as.matrix(x),path, kcdf="Gaussian",method = "ssgsea",parallel.sz=30) 
}
pan_cancer_immune_cell_infiltration=lapply(pan_cancer_data, GSVA)
save(pan_cancer_immune_cell_infiltration,file=paste0(dir,'/11_immune_cell_infiltration/pan_cancer_immune_cell_infiltration.Rdata'))

############# correlation between immune_cell infiltration and Cluster10

load(file=paste0(dir,'/11_immune_cell_infiltration/pan_cancer_immune_cell_infiltration.Rdata'))
load(paste0(dir,'/10_PanCancer_signature/pan_cancer_cluster_marker.Rdata'))
cor_function=function(expr1,expr2){
  sameID=intersect(colnames(expr1),colnames(expr2))
  expr1=expr1[,sameID]
  expr2=expr2[,sameID]
  corFilter=0            #相关系数过滤标准
  pvalueFilter=1
  outTab=data.frame()
  for(i in rownames(expr1)){
    for(j in rownames(expr2)){
      x=as.numeric(expr1[i,])
      y=as.numeric(expr2[j,])
      corT=cor.test(x,y,method = "spearm")
      cor=corT$estimate
      pvalue=corT$p.value
      outTab=rbind(outTab,cbind(gene=i,pathways=j,cor,pvalue))
    }
  }
  return(outTab)
}
pan_immune_cor_result=mapply(cor_function,pan_cancer_immune,pan_cancer_cluster_marker,SIMPLIFY = F)
pan_immune_cor_result=data.table::rbindlist(pan_immune_cor_result,use.names=TRUE, fill=TRUE, idcol="dataset") 
pan_immune_cor_result=as.data.frame(pan_immune_cor_result)
save(pan_immune_cor_result,file = paste0(dir,'/11_immune_cell_infiltration/pan_immune_cor_result.Rdata'))
