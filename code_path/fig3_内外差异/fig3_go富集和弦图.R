
library(dplyr)
#library(SeuratDisk)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
library(anndata)
library("Matrix")
#library(Seurat)
library(future.apply)
library(ggplot2)
library(clusterProfiler)
library(purrr)
library(rstatix)
library(ggpubr)
library(tidyr)
library(org.Hs.eg.db)

reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")

meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs

X=t(read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$X)
data= Seurat::CreateSeuratObject(X %>% as.matrix(),meta.data = meta.data)
data <- Seurat::NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data=data[,data$cell_type_sub=="CD74 B cells"]
Seurat::Idents(data)=data$TLS

diff_exp_data <- Seurat::FindMarkers(data, ident.1 = "Yes",
                                    ident.2 = "No" ,group.by="TLS",
                                    logfc.threshold=0.001)

AverageExpress=Seurat::AverageExpression(data,
                                         group.by = c( 'TLS'))$RNA
diff_exp_data$TLS_mean=AverageExpress[rownames(diff_exp_data),"Yes"]
diff_exp_data$No_TLS_mean=AverageExpress[rownames(diff_exp_data),"No"]
diff_exp_data$gene=rownames(diff_exp_data)

GO=clusterProfiler::enrichGO(
rownames(diff_exp_data)[diff_exp_data$avg_log2FC >0.1&diff_exp_data$ p_val_adj<0.05 ],
keyType = "SYMBOL",OrgDb="org.Hs.eg.db",
ont = "ALL",pvalueCutoff = 1)
a=dotplot(GO)
dfTerm =  read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/model/bioladder1/GOplot/demoData1.txt")
dfFC =  read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/model/bioladder1/GOplot/demoData2.txt")
library(GOplot)       # 用于绘图  安装方式 devtools::install_github("wencke/wencke.github.io")
library(RColorBrewer) # 颜色
l=a$data %>% group_by(Description) %>%  
dplyr::summarize(gene=unlist(geneID%>% stringr::str_split("/")))
l$logFC=diff_exp_data$avg_log2FC[ match(l$gene,diff_exp_data$gene)]
path=c("regulation of T cell activation",
"B cell activation",
"lymphocyte differentiation",
"immune response-activating signaling pathway")
l=l[l$Description%in% path,]
l=l[rank(-l$logFC) <30,]
# 整理数据
dfClean = dfTerm %>% 
  separate_rows(Genes,sep = ",") %>%
  left_join(dfFC,by=c("Genes"="Gene"))

dfPlot = chord_dat(data.frame(dfClean),  process=unique(dfClean$Term))
colnames(l)=c("Term" ,  "Genes" ,  "logFC")
l$Term=as.character(l$Term)
l[,2:3] %>%{.[!duplicated(.),]}
dfPlot = chord_dat(l[,1:2] %>%as.data.frame,l[,2:3] %>%{.[!duplicated(.),] } %>%
as.data.frame,  process=unique(l$Term))

# 绘图
pdf("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig3/CD74内外差异图.pdf")
GOChord(dfPlot,
        title="GOChord plot",      # 标题名称
        space = 0.02,              # 基因方格之间的空隙
        gene.order = "logFC",      #  基因的排序方式，可以按照"logFC", "alphabetical", "none", 
        gene.space = 0.25,         # 基因标签离图案的距离
        gene.size = 5,             # 基因标签的大小
        lfc.col=c('firebrick3', 'white','royalblue3'), # logFC图例的颜色
        ribbon.col=brewer.pal(ncol(dfPlot)-1, "Set2"),   # 条带颜色设置
        process.label = 8          # 图例标签的大小
)
dev.off()