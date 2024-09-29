library(ggplot2)
library(dplyr)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/home/dell/.miniconda/envs/py_scRNA/bin/python')
library(anndata)
library(dplyr)
library(ggplot2)
#reticulate::use_python('/home/dell/.miniconda/envs/py_scRNA/bin/python')
future::plan("multisession", workers =10)  # Reduce the number of workers
options(future.globals.maxSize = 10000 * 1024^2) # 将最大值设为600 MiB
future::plan("multisession", workers =10)  # Reduce the number of workers

meta.data=read_h5ad("output/data/adata.h5ad")$obs

color_dict = c(
"Cancer" ="#AC8BBF",
"Endothelial"= "#91D1C2CC",
"Myeloid"= "#8491B4CC",
"Fibroblast" ="#7E6148CC",
"NK & T"= "blue",
"B cell"= "red",
"DC"= "#F6D151",
"SMC"= "#D2568C"
)

a=ggplot(meta.data ,aes(x=x_slide_px,y=y_slide_px,color=cell_type))+geom_point(size=0.1)+
scale_color_manual(values=color_dict)+ ggplot2::coord_fixed()+theme_classic()
ggplot2::ggsave("output/fig1/ggsave.png",a,height=20*2,width=15*2)




#  基因原位表达 
library(dplyr)
library(Seurat)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
library(anndata)
library(dplyr)
library(ggplot2)
library("Matrix")
library(ggalluvial)
library(cowplot)
library(tidyr)
library(future)

reticulate::use_python('/home/dell/.miniconda/envs/py_scRNA/bin/python')
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
adata=read_h5ad("output/data/adata.h5ad")
X=t(adata$X)
TLS_meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
TLS_meta.data$TLS=ifelse(TLS_meta.data$dbscan_labels!=-1 ,"TLS","nTLS")  
TLS_meta.data=TLS_meta.data[TLS_meta.data$TLS=="TLS",]
meta=adata$obs
meta$TLS="nTLS"
meta$TLS[rownames(meta)%in% rownames(TLS_meta.data)]="TLS"
gene="CD19,MS4A1,CD79A,CD79B,CD8A,CD8B,CD4,CXCR5,IL21,IL6,CD3D,CD3E,CD3G,CD40,ITGAX,GNLY,MKI67,CXCL12,CXCL13,CCL21,CCL19,CCL2,CCL3,CCL4,CCL5,CCL8,CCL18"
gene=gene%>% stringr::str_split(",| ") %>% unlist %>% unique
gene=gene[gene!=""]
gene_exp=data.frame(X[gene,]%>%t %>% as.matrix,meta)
gene_exp$x=gene_exp$x_slide_px                   
gene_exp$y=gene_exp$y_slide_px                   

for(gene_name in gene){
ggplot_obj=ggplot()+
geom_point(data=gene_exp,aes_string(x="x",y="y"),color="#DCDCDC",size=0.3)+
geom_point(data=gene_exp[gene_exp$TLS %in% "TLS" ,],aes_string(x="x",y="y"),color="#A9A9A9",size=1)+
geom_point(data=gene_exp[gene_exp[,gene_name]>0.0001 &gene_exp$TLS %in% "TLS" ,],
aes_string(x="x",y="y",color=gene_name),size=3)+
  # 第一种颜色渐变
  scale_color_gradient(low =  "#FBDFE2", high ="#B83945", name = gene_name)+
   ggplot2::coord_fixed()+ ggplot2::labs(title=paste(gene_name))+ theme_classic()
ggplot2::ggsave("output/fig1/gene/"%>% paste0(gene_name,".png"),ggplot_obj,height=20*2,width=15*2)
}
gene_name="GNLY"

