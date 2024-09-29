# AUCell
library(AUCell)
library(anndata)
library(future)
library("Matrix")
library(ggplot2)
library(dplyr)
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
plan("multisession", workers =20)
plan()
set.seed(123)
h5ad=read_h5ad("output/data/adata_B_TLS.h5ad")
meta.data=h5ad$obs
X=t(h5ad$X)
dat2= Seurat::CreateSeuratObject(X[,rownames(h5ad$obs)],meta.data = h5ad$obs)
dat2@meta.data=meta.data
umap=data.frame(umap_1=h5ad$obsm["X_umap"]$X_umap[,1],umap_2=h5ad$obsm["X_umap"]$X_umap[,2])
rownames(umap)=colnames(dat2)


# AUC model run 
cells_rankings <- AUCell_buildRankings(X[,rownames(h5ad$obs)]%>%as.matrix, nCores=30, plotStats=TRUE)
gene=clusterProfiler::read.gmt("raw_data/c7.all.v2023.2.Hs.symbols.gmt")
gene=gene[gene$gene %in%rownames(cells_rankings), ]
cells_AUC <- AUCell_calcAUC(split(gene$gene,gene$term), cells_rankings,nCores = 40, aucMaxRank=nrow(cells_rankings)*0.1)
rownames(auc)=rownames(h5ad$obs)
saveRDS(cells_AUC,"/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/Bcell_AUC.RDS")
# read AUC
cells_AUC=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/Bcell_AUC.RDS")
auc_data=getAUC(cells_AUC) %>% as.data.frame()
auc_data_cell_type=data.frame(TLS=meta.data$TLS,cell_type_sub=meta.data$cell_type_sub,t(auc_data))

auc_CD74_B_cells_cell_type=auc_data_cell_type[auc_data_cell_type$cell_type_sub=="CD74 B cells",]


auc_CD74_B_cells_cell_type 