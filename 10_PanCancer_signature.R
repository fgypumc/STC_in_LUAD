library(dplyr)
library(GSVA)
################
dir='YourPath'
dir.create(paste0(dir,'/10_PanCancer_signature'))
setwd(paste0(dir,'/10_PanCancer_signature'))

######caculating signature
load(paste0(dir,'/10_PanCancer_signature/pan_cancer_cli.Rdata'))
load(file=paste0(dir,'/10_PanCancer_signature/pan_cancer_data_all.Rdata'))
load(file=paste0(dir,'/10_PanCancer_signature/pan_cancer_data.Rdata'))
load(file=paste0(dir,'/10_PanCancer_signature/cluster_marker.Rdata'))
GSVA=function(x){
  gsva(as.matrix(x),cluster_marker_list, kcdf="Gaussian",method = "ssgsea",parallel.sz=30) 
}
pan_cancer_cluster_marker_all=lapply(pan_cancer_data_all, GSVA)
pan_cancer_cluster_marker=lapply(pan_cancer_cluster_marker_all, function(x){x[,colnames(x)%in%rownames(pan_cancer_cli)]})
save(pan_cancer_cluster_marker_all,file=paste0(dir,'/10_PanCancer_signature/pan_cancer_cluster_marker_all.Rdata'))
save(pan_cancer_cluster_marker,file=paste0(dir,'/10_PanCancer_signature/pan_cancer_cluster_marker.Rdata'))
