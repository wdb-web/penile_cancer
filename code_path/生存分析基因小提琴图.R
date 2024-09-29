life_gene= "ALOX5 AMFR CD40 CD74 CD83 EVL HLA-DMA HLA-DMB HLA-DPB1 HLA-DQA1 HLA-DRA LAPTM5 MKNK2 NCF1 NFATC1 PARVG PSTPIP1 PTPRC" %>% 
stringr::str_split(" ")%>% unlist
# AUCell
library(AUCell)
library(anndata)
library(future)
library("Matrix")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
reticulate::use_python('//data/project/wangdongbin//.conda/envs/mistyr/bin/python')
plan("multisession", workers =20)
set.seed(123)
h5ad=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")
meta.data=h5ad$obs
X=t(h5ad$X)
meta.data=meta.data[meta.data$cell_type=="B cell",]


X=X[life_gene,rownames(meta.data)] %>% as.matrix%>%t
B_cell=data.frame(meta.data[,c("TLS"  ,"cell_type_sub")],  X)



for(gene in colnames(B_cell)[-c(1:2)]){
    bxp <- ggplot(B_cell,aes_string(x = "TLS", y = gene , color = "TLS"),
)+geom_violin()+ ggsci::scale_color_aaas()+
 # ,facet.by = "cell_type_sub"
    ggplot2::facet_grid(~cell_type_sub)
stat.test <-   B_cell %>%
  group_by(cell_type_sub) %>%
  t_test( as.formula(paste0(gene,"~TLS ")) ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "TLS")

stat.test <- stat.test %>% add_xy_position(x = "TLS")
boxplot_plot= bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}", tip.length = 0
  ) + theme_classic()+
  scale_y_continuous(expand = expansion(mult = c(0, .5)))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/fig2_violin"%>% paste0(gene,".pdf"),boxplot_plot)

    bxp <- ggplot(B_cell,aes_string(x = "TLS", y = gene , color = "TLS"),
)+geom_boxplot()+ ggsci::scale_color_aaas()+
 # ,facet.by = "cell_type_sub"
    ggplot2::facet_grid(~cell_type_sub)
stat.test <-   B_cell %>%
  group_by(cell_type_sub) %>%
  t_test( as.formula(paste0(gene,"~TLS ")) ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "TLS")
stat.test <- stat.test %>% add_xy_position(x = "TLS")
boxplot_plot= bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}", tip.length = 0
  ) + theme_classic()+
  scale_y_continuous(expand = expansion(mult = c(0, .5)))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/fig2_boxplot"%>% paste0(gene,".pdf"),boxplot_plot)
}


# AUCell
library(AUCell)
library(anndata)
library(future)
library("Matrix")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")
diff_exp_data=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/diff_exp_data.RDS")

diff_exp_data=diff_exp_data%>% filter(avg_log2FC >0.1 & p_val<0.05)  %>%
  mutate(cell_type_sub_gene= paste(gene,cluster ))
diff_exp_data=diff_exp_data[diff_exp_data$pct.1>0.2,]
reticulate::use_python('//data/project/wangdongbin//.conda/envs/mistyr/bin/python')
plan("multisession", workers =20)
set.seed(123)
h5ad=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")
meta.data=h5ad$obs
X=t(h5ad$X)
meta.data=meta.data[meta.data$cell_type=="B cell",]
Seurat_data=Seurat::CreateSeuratObject(X[,rownames(meta.data)],meta.data = meta.data)
Seurat_data$cell_type_sub


Seurat_data <- NormalizeData(Seurat_data)

AverageExpress=Seurat::AverageExpression(Seurat_data,group.by = c( 'cell_type_sub'),features=diff_exp_data$gene %>%unique
)$RNA
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/B_cell_heatmap.pdf",width = 3,height = 15)

pheatmap::pheatmap(AverageExpress%>%t%>% scale %>%t, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100) %>%rev
                   )

dev.off()


