library(dplyr)
library(Seurat)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/data/project/wangdongbin/.conda/envs/mistyr/bin/python')
library(anndata)
library(dplyr)
library("Matrix")
library(ggplot2)
library(ggalluvial)
library("Matrix")
library(cowplot)
library(scatterpie)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(survival)
library(scatterpie)
library(clusterProfiler)
library(enrichplot)
library(AUCell)
library(ggpubr)
library(reticulate)
Sys.setenv(PATH = paste("/data/project/wangdongbin/.miniconda/bin", Sys.getenv("PATH"), sep = ":"))
library(reticulate)
use_condaenv("/data/project/wangdongbin/.miniconda/envs/py_scRNA", conda = "/data/project/wangdongbin/.miniconda/bin/conda")
.libPaths("/home/dell/.miniconda/envs/R_scRNA/lib")

Sys.setenv(RETICULATE_PYTHON = '/data/project/wangdongbin/.conda/envs/mistyr/bin/python')
reticulate::use_python('/data/project/wangdongbin/.conda/envs/mistyr/bin/python')
reticulate::use_python('/data/project/wangdongbin/.miniconda/envs/py_scRNA/bin/python')



adata=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")


gene_exp_pie=adata$obs

gene_exp_pie=gene_exp_pie%>% group_by(cell_type_sub,TLS,
  TNM_stage=if_else(TNM_stage=="I","Early","Late"),sample
) %>%
  summarise(n=n())%>% ungroup%>%group_by(sample) %>%
  mutate(rate=n/sum(n)) %>% ungroup %>%group_by(cell_type_sub,TLS#, TNM_stage
  ) %>% summarise(rate=mean(rate))

tidyr::pivot_wider(values_from=n,names_from=TNM_stage)
gene_exp_pie[is.na(gene_exp_pie)]=0
# 比例柱状图
color_dict=c("CD74 B cells"="#FA7F73","IGKV4 Plasma cells"="#8DD1C6","IGLV4 Plasma cells"="#B1DD6D",
            "IGHV1 Plasma cells"="#BBB8D9","IGHM Plasma cells"="#FCB264")
color_dict=c("#FA7F73","#8DD1C6")

gene_exp_pie=gene_exp_pie[gene_exp_pie$TLS=="Yes",]
gene_exp_pie$cell_type_sub <- factor(gene_exp_pie$cell_type_sub, levels = gene_exp_pie$cell_type_sub[order(gene_exp_pie$rate,decreasing =T)])
levels= gene_exp_pie$cell_type_sub %>% levels
b=gene_exp_pie%>%ggplot( aes(x =cell_type_sub , y = rate,fill=cell_type_sub)) +
   geom_bar(stat = "identity")+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_fill_manual(values=color_dict)+
  labs(
       x = "TNM Stage",
       y = "Mean Rate")

ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig2/fig2c_TLS不同细胞比例.pdf",
       b,width=8,height=6)


gene_exp_pie=adata$obs

gene_exp_pie=gene_exp_pie%>% group_by(cell_type_sub,TLS,
  TNM_stage=if_else(TNM_stage=="I","Early","Late"),sample
) %>%
  summarise(n=n())%>% ungroup%>%group_by(sample) %>%
  mutate(rate=n/sum(n)) %>% ungroup %>%group_by(cell_type_sub,TLS#, TNM_stage
  ) 

tidyr::pivot_wider(values_from=n,names_from=TNM_stage)
gene_exp_pie[is.na(gene_exp_pie)]=0
# 比例柱状图
color_dict=c("CD74 B cells"="#FA7F73","IGKV4 Plasma cells"="#8DD1C6","IGLV4 Plasma cells"="#B1DD6D",
            "IGHV1 Plasma cells"="#BBB8D9","IGHM Plasma cells"="#FCB264")

gene_exp_pie=gene_exp_pie[gene_exp_pie$TLS=="Yes",]
gene_exp_pie$cell_type_sub <- factor(gene_exp_pie$cell_type_sub, levels = c("CD74 B cells","IGKV4 Plasma cells","IGLV4 Plasma cells","IGHV1 Plasma cells","IGHM Plasma cells"))

b=gene_exp_pie%>%ggplot( aes(x =cell_type_sub , y = rate,color=cell_type_sub)) +
   geom_boxplot()+
  theme_classic()+#ggpubr::stat_compare_means(method = "t.test")+
  scale_color_manual(values=color_dict)+
  labs(
       x = "TNM Stage",
       y = "Mean Rate")

ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig2/fig2c_TLS不同细胞比例_boxplot.pdf",
       b,width=8,height=6)


#  细胞浸润

library(dplyr)
library(Seurat)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/data/project/wangdongbin/.conda/envs/mistyr/bin/python')
reticulate::use_python('/data/project/wangdongbin/.miniconda/envs/py_scRNA/bin/python')

library(anndata)
library(dplyr)
library("Matrix")
library(ggplot2)
library(ggalluvial)
library("Matrix")
library(cowplot)
library(scatterpie)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(survival)
library(scatterpie)
library(clusterProfiler)
library(enrichplot)
library(AUCell)
library(ggpubr)
library(reticulate)
#https://github.com/Moonerss/CIBERSORT
library(CIBERSORT)
library(GSVA)
#library(ImmuCellAI)
#library(xCell)
setwd( "/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS")
adata=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")

#Sys.setenv(PATH = paste("/data/project/wangdongbin/.miniconda/bin", Sys.getenv("PATH"), sep = ":"))
library(reticulate)
mixture_data=read.csv("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/raw_data/RNA_OS.csv") %>% na.omit
mixture_data_OS=read.csv("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/raw_data/RNA_OS.csv")[,c(1:4)] %>% na.omit

rownames(mixture_data)=mixture_data$X
mixture_matrix=mixture_data[,-c(1:4)] %>%t %>% as.matrix
#devtools::install_local("/home/wangdongbin/data/data/r-immunedeconv-2.1.2-r43hdfd78af_2.tar.bz2")
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB


exp=adata$X %>% as.matrix
sig_matrix=data.frame(cell_type= adata$obs[,'cell_type_sub'],exp)
#　CIBERSORT
sig_matrix_mean=sig_matrix%>%group_by(cell_type) %>% summarise_if(is.numeric,mean) %>% data.frame(chech.name=FALSE)
rownames(sig_matrix_mean)=sig_matrix_mean[,1]
sig_matrix_mean=sig_matrix_mean[,-1] %>%t
sig_matrix_mean_matrix <- sig_matrix_mean %>% as.matrix

#sig_matrix, mixture_file, 
options(mc.cores =10)
cibersort_results <- cibersort(sig_matrix=sig_matrix_mean_matrix, mixture_file=mixture_matrix, QN = FALSE)

# ssgsea 类型
x=t(adata$X) %>%as.matrix
Seurat_Object=Seurat::CreateSeuratObject(counts=x,meta.data=adata$obs)
Seurat_Object=Seurat::NormalizeData(Seurat_Object)
Seurat_Object <- NormalizeData(Seurat_Object)

Seurat::Idents(Seurat_Object)=adata$obs[,'cell_type_sub']
marker=Seurat::FindAllMarkers(Seurat_Object)
marker_list=marker[marker$avg_log2FC>0.8 &marker$p_val_adj <0.01, ] %>% group_by(cluster) %>% 
    dplyr::arrange(avg_log2FC) %>% top_n(n=10,wt=-p_val_adj) %>%{split(.$gene ,.$cluster)}  #group_by(cluster) %>% group_split()
#ssgsea_results <- gsva(mixture_matrix, marker_list, method = "ssgsea")
gsvaPar <- gsvaParam(mixture_matrix, marker_list, kcdf = "Gaussian")
ssgsea_results <- gsva(gsvaPar, verbose=FALSE)%>%as.data.frame %>%t

# AUCell 类型

cells_rankings <- AUCell_buildRankings(mixture_matrix, nCores=5, plotStats=TRUE)


cells_AUC <- AUCell_calcAUC(marker_list, cells_rankings,nCores = 5, aucMaxRank=nrow(cells_rankings)*0.1)
cells_AUC_data=AUCell::getAUC(cells_AUC) %>% as.data.frame%>%t
cells_AUC_data_OS=data.frame(mixture_data_OS[,c(2,3)],cells_AUC_data)

for(gene in colnames(cells_AUC_data_OS)[3:7]){
  cat(gene)
  try({
OS=cells_AUC_data_OS
st2=OS[,c("futime","fustat",gene)] %>% data.frame(  check.names=F) %>% na.omit
st2$OS_STATUS=st2$fustat %>% as.factor()%>% as.numeric()
#st2=st2[!stringr::str_detect(st2$group,"PM"),]
#将结果变量中缺失数据删除，读者根据自己数据特点决定是否需要此命令。
#Fit survival ROC model with method of "KM"
library(survival)
library(survminer)
# SROC= survivalROC(Stime = st2$futime,
#                   status = st2$OS_STATUS %>% as.factor()%>% as.numeric(),
#                   marker = st2[,gene],
#                   predict.time = 150, method= "KM" ) #构建生存函数
# cut.op= SROC$cut.values[which.max(SROC$TP-SROC$FP)] #计算最佳截断值
# cut.op#输出最佳截断值

st2[,gene]=if_else(st2[,gene]>mean(st2[,gene]),"High","Low")
# 绘制Kaplan-Meier生存曲线
km_fit <- survfit(Surv(st2$futime/30,  st2$OS_STATUS) ~ st2[,gene], data = st2)
names(km_fit$strata)=c(paste(gene,"High" ,table(st2[,gene])[1]),paste(gene,"Low" ,table(st2[,gene])[2]) )
# 绘制图形

#pdf(paste0(gene,".pdf"))
surv_plot=ggsurvplot(km_fit,       pval = TRUE,
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme =  cowplot::theme_cowplot(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
)

cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/免疫浸润/AUC_",gene,".pdf"), 
                 plot =cowplot::plot_grid(surv_plot$plot,surv_plot[["table"]], 
                                          nrow = 2,rel_heights = c(1,0.5)), width = 8, height = 6)
})
} 


ssgsea_results_OS=data.frame(mixture_data_OS[,c(2,3)],ssgsea_results)

for(gene in colnames(ssgsea_results_OS)[3:7]){
  cat(gene)
  try({
OS=ssgsea_results_OS
st2=OS[,c("futime","fustat",gene)] %>% data.frame(  check.names=F) %>% na.omit
st2$OS_STATUS=st2$fustat %>% as.factor()%>% as.numeric()
#st2=st2[!stringr::str_detect(st2$group,"PM"),]
#将结果变量中缺失数据删除，读者根据自己数据特点决定是否需要此命令。
#Fit survival ROC model with method of "KM"
library(survival)
library(survminer)
# SROC= survivalROC(Stime = st2$futime,
#                   status = st2$OS_STATUS %>% as.factor()%>% as.numeric(),
#                   marker = st2[,gene],
#                   predict.time = 150, method= "KM" ) #构建生存函数
# cut.op= SROC$cut.values[which.max(SROC$TP-SROC$FP)] #计算最佳截断值
# cut.op#输出最佳截断值

st2[,gene]=if_else(st2[,gene]>mean(st2[,gene]),"High","Low")
# 绘制Kaplan-Meier生存曲线
km_fit <- survfit(Surv(st2$futime/30,  st2$OS_STATUS) ~ st2[,gene], data = st2)
names(km_fit$strata)=c(paste(gene,"High" ,table(st2[,gene])[1]),paste(gene,"Low" ,table(st2[,gene])[2]) )
# 绘制图形

#pdf(paste0(gene,".pdf"))
surv_plot=ggsurvplot(km_fit,       pval = TRUE,
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme =  cowplot::theme_cowplot(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
)

cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/免疫浸润/ssgsea_",gene,".pdf"), 
                 plot =cowplot::plot_grid(surv_plot$plot,surv_plot[["table"]], 
                                          nrow = 2,rel_heights = c(1,0.5)), width = 8, height = 6)
})
} 



cibersort_results_OS=data.frame(mixture_data_OS[,c(2,3)],cibersort_results)

for(gene in colnames(cibersort_results_OS)[3:7]){
  cat(gene)
  try({
OS=cibersort_results_OS
st2=OS[,c("futime","fustat",gene)] %>% data.frame(  check.names=F) %>% na.omit
st2$OS_STATUS=st2$fustat %>% as.factor()%>% as.numeric()
#st2=st2[!stringr::str_detect(st2$group,"PM"),]
#将结果变量中缺失数据删除，读者根据自己数据特点决定是否需要此命令。
#Fit survival ROC model with method of "KM"
library(survival)
library(survminer)
# SROC= survivalROC(Stime = st2$futime,
#                   status = st2$OS_STATUS %>% as.factor()%>% as.numeric(),
#                   marker = st2[,gene],
#                   predict.time = 150, method= "KM" ) #构建生存函数
# cut.op= SROC$cut.values[which.max(SROC$TP-SROC$FP)] #计算最佳截断值
# cut.op#输出最佳截断值

st2[,gene]=if_else(st2[,gene]>mean(st2[,gene]),"High","Low")
# 绘制Kaplan-Meier生存曲线
km_fit <- survfit(Surv(st2$futime/30,  st2$OS_STATUS) ~ st2[,gene], data = st2)
names(km_fit$strata)=c(paste(gene,"High" ,table(st2[,gene])[1]),paste(gene,"Low" ,table(st2[,gene])[2]) )
# 绘制图形

#pdf(paste0(gene,".pdf"))
surv_plot=ggsurvplot(km_fit,       pval = TRUE,
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme =  cowplot::theme_cowplot(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
)

cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/免疫浸润/cibersort_",gene,".pdf"), 
                 plot =cowplot::plot_grid(surv_plot$plot,surv_plot[["table"]], 
                                          nrow = 2,rel_heights = c(1,0.5)), width = 8, height = 6)
})

} 


# AUCell 注释
library(clusterProfiler)
library(AUCell)
library(dplyr)
reticulate::use_python('/data/project/wangdongbin/.conda/envs/mistyr/bin/python')
library(anndata)
library("Matrix")
library(ggplot2)
h5ad=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")
exprMatrix=t(h5ad$X)[,rownames(h5ad$obs)]%>%as.matrix
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=30, plotStats=TRUE)
gene=clusterProfiler::read.gmt("raw_data/genesets.v2024.1.Hs.gmt")
gene=gene[gene$term%>% stringr::str_detect("GERMINAL|_GC_"), ]

gene=gene[gene$gene %in%rownames(cells_rankings), ]
cells_AUC <- AUCell_calcAUC(split(gene$gene,gene$term), cells_rankings,nCores = 40, aucMaxRank=nrow(cells_rankings)*0.1)
umap=h5ad$obsm$X_umap %>% as.matrix
rownames(umap)=rownames(h5ad$obs)
#cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=T, nCores=36, assign=F)
#selectedThresholds <- getThresholdSelected(cells_assignment)
AUCell_data=getAUC(cells_AUC)%>% as.data.frame %>%t %>%data.frame(umap,.)
#scale_color_gradientn(colors =c ("#FF5010","#800080",   "#FFFF00"), limits = c(0, quantile(AUCell_data$exp,0.95)))

for(path in  rownames(cells_AUC)){
a=ggplot(data=AUCell_data,aes(x=X1,y=X2))+geom_point(aes_string(color=path))+
scale_color_gradientn(colors =c ("#800080",   "#FFFF00"), 
limits = c(0, quantile(AUCell_data[,path],0.95)))+
labs(title=path)+theme_classic()
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",
path%>% make.names,"umap.pdf" ),a, width = 24, height = 18)
}



path="GSE12366_GC_BCELL_VS_PLASMA_CELL_DN"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.998))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)
path="GSE12366_GC_BCELL_VS_PLASMA_CELL_UP"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.998))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)


path="GSE12845_IGD_NEG_BLOOD_VS_DARKZONE_GC_TONSIL_BCELL_DN"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.99))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)
path="GSE12845_IGD_NEG_BLOOD_VS_DARKZONE_GC_TONSIL_BCELL_DN"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.997))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 16)


path="GSE12366_GC_VS_MEMORY_BCELL_UP"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.00), quantile(AUCell_data[,path],1))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)
path="GSE12366_GC_VS_NAIVE_BCELL_UP"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.995))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)
path="GSE12845_PRE_GC_VS_DARKZONE_GC_TONSIL_BCELL_UP"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.998))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)


path="GSE12845_IGD_NEG_BLOOD_VS_PRE_GC_TONSIL_BCELL_DN"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.998))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)
GSE12845_PRE_GC_VS_DARKZONE_GC_TONSIL_BCELL_UP

path="GSE12845_PRE_GC_VS_DARKZONE_GC_TONSIL_BCELL_UP"

a=ggplot(AUCell_data, aes(X1, X2))  +
  geom_point(aes_string(color=path),size=0.5)+ 
  viridis::scale_color_viridis(option="E",limits = c(quantile(AUCell_data[,path],0.05), quantile(AUCell_data[,path],0.91))) +
  theme(legend.position = "none") + theme_bw()+ #ggplot2::facet_wrap(~HPV)+
  xlab(path)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/GC_AUC/",path,".pdf"),a,width= 14)
GSE12845_PRE_GC_VS_DARKZONE_GC_TONSIL_BCELL_UP
