###############
library(parallel)
library(tidyverse)
cl.cores <- detectCores()
cl <- makeCluster(30)
library(metacell)
library(Seurat)
####################
dir='YourPath'
dir.create(paste0(dir,'/3.metacell'))
setwd(paste0(dir,'/3.metacell'))
####################
if(!dir.exists("testdb")){
  dir.create("testdb/")
}
scdb_init("testdb/", force_reinit=T) 
if(!dir.exists("figs")){
  dir.create("figs/")
}
scfigs_init("figs/") 

####################running metacell analysis
sce <- readRDS(paste0(dir,"/2.inferCNV/NonImmune.rds"))
sce=subset(sce,Cell_type.refined %in% c('Epithelial cells'))
matrix <- R@assays$RNA@counts
cellid <- colnames(matrix)
cellmeta <-data.frame(cellid,type='10x',batch_set_id=paste0('test_',cellid),amp_batch_id=paste0('test_',cellid),seq_batch_id=paste0('test_',cellid),spike_count=0,row.names=1)
mat <- scm_new_matrix(matrix, cellmeta, stat_type = "umi")
scdb_add_mat('pbmc', mat) 
nms <- c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes <- c(grep("^IGJ", nms, v=T),grep("^IGH",nms,v=T),grep("^IGK", nms, v=T),grep("^IGL", nms, v=T))
bad_genes <- unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10",ig_genes))
mcell_mat_ignore_genes(new_mat_id='pbmc', mat_id='pbmc', bad_genes, reverse=F)
mcell_mat_ignore_small_cells('pbmc', 'pbmc', 800)

mcell_add_gene_stat(gstat_id='pbmc', mat_id='pbmc', force=T)
mcell_gset_filter_varmean(gset_id='pbmc_feats', gstat_id='pbmc', T_vm=0.08, force_new=T) #(会在database中生成gset.pbmc_feats.Rda）
mcell_gset_filter_cov(gset_id='pbmc_feats', gstat_id='pbmc', T_tot=100, T_top3=2)
mcell_plot_gstats(gstat_id='pbmc', gset_id='pbmc_feats') #(会在database中生成gstat.pbmc.Rda,在figs中生成三个图片pbmc.szcor.png 、pbmc.top3.png、pbmc.varmin.png）


mcell_add_cgraph_from_mat_bknn(mat_id='pbmc', gset_id='pbmc_feats', graph_id='pbmc_graph', K=100, dsamp=T) #(生成cgraph.pbmc_graph.Rda)
mcell_coclust_from_graph_resamp(coc_id='pbmc_coclust', graph_id='pbmc_graph', min_mc_size=20, p_resamp=0.75, n_resamp=500)
mcell_mc_from_coclust_balanced(coc_id='pbmc_coclust',mat_id='pbmc', mc_id='pbmc_mc',K=30, min_mc_size=30,alpha=2) #(生成coclust.pbmc_coclust.Rda)


mcell_plot_outlier_heatmap(mc_id='pbmc_mc', mat_id = 'pbmc', T_lfc=3) #(会在figs中生成pbmc_mc.outlier.png)
mcell_mc_split_filt(new_mc_id='pbmc_mc_f', mc_id='pbmc_mc', mat_id='pbmc', T_lfc=3, plot_mats=F) #（会在database中生成mc.pbmc_mc_f.Rda）
mcell_gset_from_mc_markers(gset_id='pbmc_markers',mc_id='pbmc_mc_f')


mc_colorize_default('pbmc_mc_f')
mcell_mc2d_force_knn(mc2d_id='pbmc_2dproj',mc_id='pbmc_mc_f',graph_id='pbmc_graph')
tgconfig::set_param("mcell_mc2d_cex",1, "metacell")  #设置metacell在图中的点大小
tgconfig::set_param("mcell_mc2d_height",1100, "metacell") #设置图片的高度
tgconfig::set_param("mcell_mc2d_width",1100, "metacell") #设置图片的宽度
mcell_mc2d_plot('pbmc_2dproj') #(会在figs中生成pbmc_2dproj.2d_graph_proj.png)

mc_hc = mcell_mc_hclust_confu(mc_id="pbmc_mc_f",
                              graph_id="pbmc_graph")
sapply(.scdb,length)

mc_sup = mcell_mc_hierarchy(mc_id="pbmc_mc_f",
                            mc_hc=mc_hc, T_gap=0.04)
sapply(.scdb,length)

mcell_mc_plot_hierarchy(mc_id="pbmc_mc_f",
                        graph_id="pbmc_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800, heigh=2000, min_nmc=2)


#####extracting new matrix for infiltrated cells
load( paste0(dir,'/3.metacell/testdb/mc.pbmc_mc_f.Rda'))
cluster=data.frame(cluster=object@mc)
sce <- readRDS(paste0(dir,"/2.inferCNV/NonImmune.rds"))
sce=subset(sce, cells=as.character(rownames(cluster)))
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(object = sce)
sce <- RunPCA(object = sce,features = VariableFeatures(object = sce))
#sce <- JackStraw(sce, num.replicate = 100)
#sce <- ScoreJackStraw(sce, dims = 1:20)
#ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.7)
sce <- RunUMAP(sce, dims = 1:10)
saveRDS(sce,paste0(dir,'/3.metacell/metacell_seurat_ec.rds'))

############extracting new matrix for metacell
mat=as.matrix(sce@assays$RNA@data)  
sameID=intersect(rownames(cluster),colnames(mat))
rt=data.frame(cluster=cluster[sameID,],t(mat[,sameID]))
new.data<-aggregate(list(rt[,2:length(colnames(rt))]),
                    list(rt[,1]),FUN=mean)
head(new.data)[1:6,1:6]
meta_df=as.data.frame(new.data)
meta_df$Group.1=paste0('C', meta_df$Group.1)
rownames(meta_df)=meta_df$Group.1
meta_df=meta_df[,-1]
head(meta_df)[1:6,1:6]
save(meta_df,file=paste0(dir,'/3.metacell/metacell_exp.Rdata'))
