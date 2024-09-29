

#/home/wangdongbin/.conda/envs/mistyr/bin/R
library(dplyr)
#library(SeuratDisk)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
library(anndata)
library("Matrix")
#library(Seurat)
library(future.apply)
library(ggplot2)

library(purrr)
library(rstatix)
library(ggpubr)
library(tibble)
library(tidyr)
library(clusterProfiler)
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")
#/data/project/wangdongbin/.miniconda/envs/R_scRNA/bin/R

library(dplyr)
library(anndata)
library("Matrix")
library(ggplot2)
h5ad=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")
B_meta.data=h5ad$obs
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
Tcell=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
meta.data<- Tcell$obs
x=t(Tcell$X)

a.meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  

gene_exp_pie=meta.data

gene_exp_pie=gene_exp_pie %>%
  group_by(cell_subtype, TLS) %>%
  summarise(n = dplyr::n())%>% ungroup%>%#group_by(sample_id) %>%
  mutate(rate=n/sum(n)) %>% ungroup %>%group_by(cell_subtype,TLS#, TNM_stage
  ) %>% summarise(rate=mean(rate)) %>%na.omit

gene_exp_pie=gene_exp_pie[gene_exp_pie$TLS=="TLS",]
gene_exp_pie$cell_subtype <- factor(gene_exp_pie$cell_subtype, levels = gene_exp_pie$cell_subtype[order(gene_exp_pie$rate,decreasing =T)])
levels= gene_exp_pie$cell_subtype %>% levels
color_dict=c(
 'CTL'= '#8c564b',
 'Naive T'= '#e377c2',
 'NK'= '#bcbd22',
 'Treg'= '#17becf',
 'Proliferative T'= '#1f77b4',
 'γδ T'= '#ff7f0e',
 'Th1'= '#2ca02c',
 "CD74 B cells"="red",
 'Th17'= 'blue'
 
 )
b=gene_exp_pie%>%ggplot( aes(x =cell_subtype , y = rate,fill=cell_subtype)) +
   geom_bar(stat = "identity")+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)+
  labs(
       x = "TNM Stage",
       y = "Mean Rate")

ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig6/fig6c_TLS不同细胞比例.pdf",
       b,width=8,height=6)

b=meta.data%>% na.omit%>%ggplot( aes(x =x_slide_mm , y = y_slide_mm,color=cell_subtype)) +
   geom_point(data=B_meta.data[B_meta.data$cell_type_sub=="CD74 B cells" & B_meta.data$TLS=="Yes",],aes(x =x_slide_mm , y = y_slide_mm,color="CD74 B cells"))+
 geom_point()+
  theme_classic()+

  #ggpubr::stat_compare_means(method = "t.test")+
  scale_color_manual(values=color_dict)+ggplot2::coord_fixed()

ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig6/fig6b_TLS原位.pdf",
       b,width=8*3,height=12*3,limitsize = FALSE)

# %% 
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) 
data= Seurat::CreateSeuratObject(x %>% as.matrix(),meta.data = meta.data)
data <- Seurat::NormalizeData(data[ ,data$cell_subtype=="Naive T"], normalization.method = "LogNormalize", scale.factor = 10000)

diff_exp_data <- Seurat::FindMarkers(data, ident.1 = "TLS",
                                    ident.2 = "nTLS" ,group.by="TLS",
                                    logfc.threshold=0.001)

                                    
library(rlang)
scale_color <- function(.x) {
 .x
}
diff_exp_data$p_val_adj[diff_exp_data$p_val_adj==0]=min(diff_exp_data$p_val_adj[diff_exp_data$p_val_adj!=0])
diff_exp_data$gene=""
diff_exp_data$gene[diff_exp_data$avg_log2FC >0.6&diff_exp_data$ p_val_adj<0.05 ]=rownames(diff_exp_data)[diff_exp_data$avg_log2FC >0.6&diff_exp_data$ p_val_adj<0.05 ]

a=ggplot(diff_exp_data,aes(x=avg_log2FC,y=pct.1-pct.2 ,label=gene ))+
 geom_point(aes(size=-log10(p_val_adj),
                                  color=-log10(p_val_adj)))+
  theme_classic()+ #geom_text()+
ggrepel::geom_text_repel()+
  scale_colour_gradientn(
      colours= c( "#9FBBD5","#B8474D","#B8474D"),values =c(0,25,50),
      rescaler=as_function(~ scale_color(.x))
    )+#scale_size(breaks =c(0,40), range = c(1, 6))+
  scale_x_continuous(limits = c(-1,1))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig6/fig6e_dotplot.pdf",a, width = 7, height = 6)

# 富集分析
GO=gseGO(
setNames(diff_exp_data$avg_log2FC[diff_exp_data$avg_log2FC >0.1&diff_exp_data$ p_val_adj<0.05 ],
rownames(diff_exp_data)[diff_exp_data$avg_log2FC >0.1&diff_exp_data$ p_val_adj<0.05 ]) %>% sort(decreasing=T),keyType = "SYMBOL",OrgDb="org.Hs.eg.db",
ont = "ALL",pvalueCutoff = 1)
GO=clusterProfiler::enrichGO(
rownames(diff_exp_data)[diff_exp_data$avg_log2FC >0.1&diff_exp_data$ p_val_adj<0.05 ],keyType = "SYMBOL",OrgDb="org.Hs.eg.db",
ont = "ALL",pvalueCutoff = 1)

ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig6/fig6f_富集.pdf",dotplot(GO, showCategory = 15), width = 7, height = 6)

# AUCell 类型

cells_rankings <- AUCell_buildRankings(mixture_matrix, nCores=5, plotStats=TRUE)


cells_AUC <- AUCell_calcAUC(marker_list, cells_rankings,nCores = 50, aucMaxRank=nrow(cells_rankings)*0.1)



writexl::write_xlsx(list(cibersort_results=cibersort_results %>% as.data.frame %>% tibble::rownames_to_column("sample"),
ssgsea_results=ssgsea_results %>% as.data.frame%>% tibble::rownames_to_column("sample"),
AUCell=AUCell::getAUC(cells_AUC) %>% as.matrix %>% as.data.frame%>% tibble::rownames_to_column("sample")),
"/home/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/output/CTL_score/score.xlsx")
# %%
# AUcell


#/home/wangdongbin/.conda/envs/mistyr/bin/R
library(dplyr)
#library(SeuratDisk)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
library(anndata)
library("Matrix")
#library(Seurat)
library(future.apply)
library(ggplot2)
library(purrr)
library(rstatix)
library(ggpubr)
library(tibble)
library(tidyr)
library(dplyr)
library(plyr)
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")

future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
Tcell=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
meta.data<- Tcell$obs
x=t(Tcell$X)

a.meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  

library(AUCell)
library(GSEABase)
library(rlang)
data=as.matrix(x)[,meta.data$cell_subtype=="Naive T"]
cells_rankings <- AUCell_buildRankings(data[,!is.na(colnames(data))], nCores=30, plotStats=TRUE) 
c2 <- read.gmt("/home/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/output/fig2/ssgsea/pathway_gmt.gmt") 
c2=c2[c2$gene %in%rownames(cells_rankings),]

geneSets <- split(c2$gene,c2$term)#  lapply(unique(c2$term), function(x){print(x);c2$gene[c2$term == x]})

# Build the rankings without specifying nCores

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores = 50, aucMaxRank=nrow(cells_rankings)*0.1)


getAUC(cells_AUC)%>%rownames
TLS=meta.data$TLS[match(colnames(getAUC(cells_AUC)),rownames(meta.data))]

auc_data=getAUC(cells_AUC)%>%data.frame %>%t %>%data.frame(TLS ,.)%>%tibble::rownames_to_column("cell")


mean_df <- auc_data %>%
  group_by(TLS) %>%
  summarise_if(is.numeric, mean)

# 计算 p 值
p_values <- auc_data[,-1] %>%
  pivot_longer(!c(TLS), names_to = "gene_set", values_to = "expression") %>%
  group_by(gene_set) %>%
  summarise(p_value = t.test(expression ~ TLS, data = .)$p.value)
auc_data_long= auc_data[,-1] %>%
  pivot_longer(!c(TLS), names_to = "gene_set", values_to = "expression") %>%
  group_by(gene_set) 
aucell_test=ddply(auc_data_long, .(gene_set), 
function(x) { wilcox.test(expression~TLS, data=x, paired=FALSE)$p.value })

aucell_test$TLS= mean_df[1,-1] %>% unlist
aucell_test$nTLS= mean_df[2,-1] %>% unlist
aucell_test$fdr=p.adjust(aucell_test$V1 )
aucell_test=aucell_test[aucell_test$fdr<0.05,]
aucell_test$gene_set[aucell_test$TLS>aucell_test$nTLS]

color_dict=c("#FA7F73","#8DD1C6")
ggplot_list=list()
ggplot_list[[1]]=auc_data %>% ggplot(aes(x=TLS,y=BIOCARTA_TCR_PATHWAY,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)

ggplot_list[[2]]=auc_data %>% ggplot(aes(x=TLS,y=KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)


ggplot_list[[3]]=auc_data %>% ggplot(aes(x=TLS,y=GOBP_ACTIVATION_OF_IMMUNE_RESPONSE,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)

ggplot_list[[4]]=auc_data %>% ggplot(aes(x=TLS,y=GOBP_IMMUNE_RESPONSE,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)

ggplot_list[[5]]=auc_data %>% ggplot(aes(x=TLS,y=GOBP_B_CELL_MEDIATED_IMMUNITY,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)


# 下调
 

ggplot_list[[6]]=auc_data %>% ggplot(aes(x=TLS,y=GSE36476_CTRL_VS_TSST_ACT_40H_MEMORY_CD4_TCELL_OLD_DN,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)
ggplot_list[[7]]=auc_data %>% ggplot(aes(x=TLS,y=GSE36476_CTRL_VS_TSST_ACT_72H_MEMORY_CD4_TCELL_OLD_DN ,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)
ggplot_list[[8]]=auc_data %>% ggplot(aes(x=TLS,y=GSE3982_DC_VS_MAC_DN ,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)
ggplot_list[[9]]=auc_data %>% ggplot(aes(x=TLS,y=GSE3982_MEMORY_CD4_TCELL_VS_TH1_DN ,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)
ggplot_list[[10]]=auc_data %>% ggplot(aes(x=TLS,y=GSE3982_MEMORY_CD4_TCELL_VS_TH2_DN  ,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)
 
ggplot_list[[11]]=auc_data %>% ggplot(aes(x=TLS,y=GSE3982_EFF_MEMORY_CD4_TCELL_VS_TH1_DN  ,fill=TLS))+geom_boxplot()+
stat_compare_means()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)


cowplot::ggsave2("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig6/fig6g_AUcell.pdf",
 cowplot::plot_grid(plotlist =ggplot_list), width = 12, height =12 )

# %% 
# 原位
auc_data=getAUC(cells_AUC)%>%data.frame %>%t %>%
data.frame(TLS ,x=meta.data$x_slide_mm[match(colnames(getAUC(cells_AUC)),rownames(meta.data))],
y=meta.data$y_slide_mm[match(colnames(getAUC(cells_AUC)),rownames(meta.data))]
,.)%>%tibble::rownames_to_column("cell")
a.meta.data=read_h5ad("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/data/adata.h5ad")$obs

a=auc_data %>% ggplot(aes(x=x,y=y,
color=TLS))+
geom_point()+
  theme_classic()
a=auc_data %>% ggplot(aes(x=x,y=y,
color=GOBP_IMMUNE_RESPONSE))+
geom_point(data=a.meta.data,aes(x=x_slide_mm,y=y_slide_mm),color="grey")+
geom_point()+

  theme_classic()+
  viridis::scale_color_viridis(option="E",
  limits = c(quantile(auc_data$GOBP_IMMUNE_RESPONSE,0.05),
  quantile(auc_data$GOBP_IMMUNE_RESPONSE,0.95))) +
  ggplot2::coord_fixed()+ ggplot2::labs(title="GOBP_IMMUNE_RESPONSE")
ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig6/AUCell.pdf",
a,width = 25, height =17)
pdf("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig3/CD74内外差异图.pdf",width = 25/2, height =28/2)
 GOChord(dfPlot,
        title="GOChord plot",      # 标题名称
        space = 0.02,              # 基因方格之间的空隙
        gene.order = "logFC",      #  基因的排序方式，可以按照"logFC", "alphabetical", "none", 
        gene.space = 0.25,         # 基因标签离图案的距离
        gene.size = 5,             # 基因标签的大小
        lfc.col=c( 'white','firebrick3')%>%rev , # logFC图例的颜色
        ribbon.col=brewer.pal(ncol(dfPlot)-1, "Set2"),   # 条带颜色设置
        process.label = 8          # 图例标签的大小
)
dev.off()