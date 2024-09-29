library(hdf5r)
library(dplyr)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
library(anndata)
library(dplyr)
library("Matrix")
library(monocle)
library(progressr)

# 创建进度条
future::plan("multisession", workers =10)  # Reduce the number of workers
options(future.globals.maxSize = 200 * 1024^3)
h5ad=read_h5ad("output/data/adata_B_TLS.h5ad")
meta.data=h5ad$obs
X=t(h5ad$X)
dat2= Seurat::CreateSeuratObject(X[,rownames(h5ad$obs)],meta.data = h5ad$obs)
dat2@meta.data=meta.data
umap=data.frame(umap_1=h5ad$obsm["X_umap"]$X_umap[,1],umap_2=h5ad$obsm["X_umap"]$X_umap[,2])
rownames(umap)=colnames(dat2)


# 创建新的 CellDataSet 对象
cds <- newCellDataSet(cellData  = as.matrix(X[Matrix::rowMeans(X[,rownames(meta.data)])>0.01,rownames(meta.data)]),
                      phenoData = new("AnnotatedDataFrame", data = meta.data))
# 筛选表达至少在3个细胞中非零的基因
cds <- detectGenes(cds, min_expr = 0.1)

cds <- estimateSizeFactors(cds)

#cds <- estimateDispersions(cds, minDisp = 1e-8, maxDisp = 1e8, verbose = FALSE)
cds <- estimateDispersions(cds)
#cds <- estimateDispersions(cds, fitType = "mean", verbose = FALSE)

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, 
                                   mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

cds <- reduceDimension(cds, max_components =2,reduction_method='DDRTree',num_paths=1)# ncenter

cds <- orderCells(cds)





