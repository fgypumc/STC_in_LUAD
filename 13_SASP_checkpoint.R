library(dplyr)
library(GSVA)
################
dir='YourPath'
dir.create(paste0(dir,'/13_SASP_checkpoint'))
setwd(paste0(dir,'/13_SASP_checkpoint'))

###SASP factors
load(file=paste0(dir,'/10_PanCancer_signature/pan_cancer_data.Rdata'))
load(file=paste0(dir,'/10_PanCancer_signature/cluster_marker.Rdata'))
SASP_genes=function(x){
  load(paste0(dir,'/13_1.CS_signature/sen_signature.Rdata'))#27个基因
  SASP=read_excel(paste0(dir,'/13_1.CS_signature/SASP_factors.xlsx'))
  expr=x[rownames(x)%in%SASP,]
  return(expr)
}
SASP_expr=lapply(pan_cancer_data, SASP_genes)
SASP_expr=lapply(SASP_expr, na.omit)
pan_SASP_cor_result=mapply(cor_function,SASP_expr,pan_cancer_cluster_marker,SIMPLIFY = F)
pan_SASP_cor_result=data.table::rbindlist(pan_SASP_cor_result,use.names=TRUE, fill=TRUE, idcol="dataset") 
save(pan_SASP_cor_result,file = paste0(dir,'/13_SASP_checkpoint/pan_SASP_cor_result.Rdata'))

####check points
load(file=paste0(dir,'/10_PanCancer_signature/pan_cancer_data.Rdata'))
load(file=paste0(dir,'/10_PanCancer_signature/cluster_marker.Rdata'))
checkpoint_genes=function(x){
  ImmuneCheckpoint=read_excel(paste0(dir,'/13_SASP_checkpoint/ImmuneCheckpoint.xlsx'))
  expr=x[rownames(x)%in%ImmuneCheckpoint$Symbol,]
  return(expr)
}
checkpoint_expr=lapply(pan_cancer_data, checkpoint_genes) 
checkpoint_expr=lapply(checkpoint_expr, na.omit)
cor_function=function(expr1,expr2){
  corFilter=0            
  pvalueFilter=1
  outTab=data.frame()
  for(i in rownames(expr1)){
    for(j in rownames(expr2)){
      x=as.numeric(expr1[i,])
      y=as.numeric(expr2[j,])
      corT=cor.test(x,y,method = "spearm")
      cor=corT$estimate
      pvalue=corT$p.value
      if(pvalue<pvalueFilter){
        outTab=rbind(outTab,cbind(gene=i,pathways=j,cor,pvalue))
      }
    }
  }
  return(outTab)
}
pan_checkpoint_cor_result=mapply(cor_function,checkpoint_expr,pan_cancer_cluster_marker,SIMPLIFY = F)
pan_checkpoint_cor_result=data.table::rbindlist(pan_checkpoint_cor_result,use.names=TRUE, fill=TRUE, idcol="dataset") 
save(pan_checkpoint_cor_result,file = paste0(dir,'/13_SASP_checkpoint/pan_checkpoint_cor_result.Rdata'))

