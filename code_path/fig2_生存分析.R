library(dplyr)
library(Seurat)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
library(anndata)
library(dplyr)
library(ggplot2)
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
library(survminer)
library(scatterpie)
library(clusterProfiler)
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")
diff_exp_data=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/diff_exp_data.RDS")

diff_exp_data=diff_exp_data%>% filter(avg_log2FC >0.1 & p_val<0.05)  %>%
  mutate(cell_type_sub_gene= paste(gene,cluster ))




diff_exp_data=diff_exp_data[diff_exp_data$pct.1>0.2,]

diff_exp_data=diff_exp_data%>% filter(avg_log2FC >0.1 & p_val<0.05,cluster=="CD74 B cells")  %>%
  mutate(cell_type_sub_gene= paste(gene,cluster ))
diff_exp_data=diff_exp_data[diff_exp_data$pct.1>0.2,]
OS=read.csv("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/raw_data/RNA_OS.csv")
# 生存分析 


source("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/code_path/KM.R")

life_p_list=list()

for(row in 1:nrow(diff_exp_data)){
  try({ gene=diff_exp_data$gene[row] %>% make.names()
  cat(gene)
  try({

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

st2[,gene]=if_else(st2[,gene]>median(st2[,gene]),"High","Low")
# 绘制Kaplan-Meier生存曲线
km_fit <- survfit(Surv(st2$futime,  st2$OS_STATUS) ~ st2[,gene], data = st2)
names(km_fit$strata)=c(paste(gene,"High" ,table(st2[,gene])[1]),paste(gene,"Low" ,table(st2[,gene])[2]) )
# 绘制图形
p=ggsurvplot2(km_fit,       pval = TRUE,
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme =  cowplot::theme_cowplot(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
)$pval

#pdf(paste0(gene,".pdf"))
life_p_list[[row]]= data.frame(gene=diff_exp_data$gene[row],p=p)

  #pdf(paste0(gene,".pdf"))
    if(p<0.05){
      surv_plot=ggsurvplot(km_fit,       pval = TRUE,
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme =  cowplot::theme_cowplot(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
)

cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/life/",diff_exp_data$cell_type_sub_gene[row],".pdf"),
                 plot =cowplot::plot_grid(surv_plot$plot,surv_plot[["table"]],
                                          nrow = 2,rel_heights = c(1,0.5)), width = 8, height = 6)

    }
}
  )
})

}

gene=do.call(rbind,life_p_list)$gene[do.call(rbind,life_p_list)$p<0.05]
write.csv(do.call(rbind,life_p_list),
  "/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig2/life/生存p值.csv")
colnames(OS)=make.names(colnames(OS))
for(gene in make.names(do.call(rbind,life_p_list)$gene[do.call(rbind,life_p_list)$p<0.05])){
  cat(gene)
  try({

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

st2[,gene]=if_else(st2[,gene]>median(st2[,gene]),"High","Low")
# 绘制Kaplan-Meier生存曲线
km_fit <- survfit(Surv(st2$futime/30/12,  st2$OS_STATUS) ~ st2[,gene], data = st2)
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

cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/life/",gene,".pdf"), 
                 plot =cowplot::plot_grid(surv_plot$plot,surv_plot[["table"]], 
                                          nrow = 2,rel_heights = c(1,0.5)), width = 8, height = 6)
})

} 
# 生存有关的基因为：

life_gene=do.call(rbind,life_p_list)$gene[do.call(rbind,life_p_list)$p<0.05]
color_dict=c("CD74 B cells"="#FA7F73","IGKV4 Plasma cells"="#8DD1C6","IGLV4 Plasma cells"="#B1DD6D",
            "IGHV1 Plasma cells"="#BBB8D9","IGHM Plasma cells"="#FCB264")
# AUCell
library(AUCell)
library(anndata)
library(future)
library("Matrix")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
plan("multisession", workers =20)
plan()
set.seed(123)
h5ad=read_h5ad("output/data/adata_B_TLS.h5ad")
meta.data=h5ad$obs
X=t(h5ad$X)
meta.data=meta.data[meta.data$cell_type=="B cell",]
X=X[life_gene,rownames(meta.data)] %>% as.matrix%>%t
B_cell=data.frame(meta.data[,c("TLS"  ,"cell_type_sub")],  X)

B_cell <- B_cell[complete.cases(B_cell), ]

# 或者仅移除 ABLIM1 列中的缺失值
# add p value 并且确定xy的位置

for(gene in colnames(B_cell)[-c(1:2)]){
    bxp <- ggplot(B_cell,aes_string(x = "cell_type_sub", y = gene , color = "cell_type_sub"),
)+geom_violin()+ scale_color_manual(values=color_dict)
stat.test <-   B_cell %>%
  t_test( as.formula(paste0(gene,"~cell_type_sub ")) ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "cell_type_sub")

stat.test <- stat.test %>% add_xy_position(x = "cell_type_sub")
boxplot_plot= bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj.signif}", tip.length = 0
  ) + theme_classic()+
  scale_y_continuous(expand = expansion(mult = c(0, .5)))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/life_gene_exp/fig2_violin"%>% paste0(gene,".pdf"),boxplot_plot,width=6)

    bxp <- ggplot(B_cell,aes_string(x = "cell_type_sub", y = gene , color = "cell_type_sub"),
)+geom_boxplot()+scale_color_manual(values=color_dict)
stat.test <-   B_cell %>%
  t_test( as.formula(paste0(gene,"~cell_type_sub ")) ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "cell_type_sub")
boxplot_plot= bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj.signif}", tip.length = 0
  ) + theme_classic()+
  scale_y_continuous(expand = expansion(mult = c(0, .5)))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig2/life_gene_exp/fig2_boxplot"%>% paste0(gene,".pdf"),boxplot_plot,width=6)
}


