library(dplyr)
library(GSVA)
################
dir='YourPath'
dir.create(paste0(dir,'/12_pathway_activities'))
setwd(paste0(dir,'/12_pathway_activities'))
#################caculate immune_pathway_activities
load(paste0(dir,'/12_pathway_activities/pathway.Rdata'))

GSVA=function(x){
  gsva(as.matrix(x),path, kcdf="Gaussian",method = "ssgsea",parallel.sz=30) 
}
pan_cancer_pathway=lapply(pan_cancer_data, GSVA)
save(pan_cancer_pathway,file=paste0(dir,'/12_pathway_activities/pan_cancer_pathway.Rdata'))

##############correlation between pathway_activities and Cluster10

load(paste0(dir,'/12_pathway_activities/pan_cancer_pathway.Rdata'))
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
pan_pathway_cor_result=mapply(cor_function,pan_pathway,pan_cancer_cluster_marker,SIMPLIFY = F)
pan_pathway_cor_result=data.table::rbindlist(pan_pathway_cor_result,use.names=TRUE, fill=TRUE, idcol="dataset") 

save(pan_pathway_cor_result,file=paste0(dir,'/12_pathway_activities/pan_pathway_cor_result.Rdata'))