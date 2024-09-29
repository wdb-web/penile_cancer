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
library(data.table)
library(survival)
library(CellChat)
library(scatterpie)
library(CellChat)
library(anndata)
library(dplyr)
library(Seurat)
library(tidyr)
library(future)

reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")

color_dict=c()


adata=read_h5ad("output/data/adata.h5ad")$obs
adata_meta.data=read_h5ad("output/data/adata.h5ad")$obs
TLS_meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
TLS_meta.data$TLS=ifelse(TLS_meta.data$dbscan_labels!=-1 ,"TLS","nTLS")  
B_meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
TLS_meta.data$cell_type=as.character(TLS_meta.data$cell_type)
TLS_meta.data$cell_type[match(rownames(B_meta.data),rownames(TLS_meta.data))]=#[meta.data$cell_type_sub=="CD74 B cells"]
as.character(B_meta.data$cell_type_sub)
T_meta.data=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")$obs
TLS_meta.data$cell_type[match(rownames(T_meta.data),rownames(TLS_meta.data))]=#[meta.data$cell_type_sub=="CD74 B cells"]
as.character(T_meta.data$cell_subtype)


TF_AUC=read.csv("output/fig7/TF/naive_T_Cell_ST_auc_mtx.csv",check.names=F)
naive_T_Cell=TLS_meta.data[TLS_meta.data$cell_type=="Naive T",] %>% na.omit

AUC_exp=data.frame(naive_T_Cell,TF_AUC[match(TF_AUC$Cell,rownames(naive_T_Cell)),],check.names=F)

AUC_exp_mean=AUC_exp%>% group_by(TLS)%>%  summarise_at(vars(67:655), mean, na.rm = TRUE)
colnames(AUC_exp_mean)[-1]  [(AUC_exp_mean[1,-1]/AUC_exp_mean[2,-1] ) >1.04]
colnames(AUC_exp_mean)
  # T 细胞增殖相关的转录因子

  # 与 T 细胞增殖、分化、激活相关的正向转录因子
positive_tf_list <- c("MYC", "STAT5A", "STAT5B", "GATA3", "TBX21", 
                      "RORC", "BCL6", "NFATC1", "NFATC2", "FOS", "JUN", 
                      "NFKB1", "NFKB2")

AUC_exp=AUC_exp%>% filter(! sample%in%c("S12","S13","S16","S17","S4","S5","S9"))
adata_meta.data=adata_meta.data%>% filter(! sample_id %in%c("S12","S13","S16","S17","S4","S5","S9"))
TLS_meta.data=TLS_meta.data[TLS_meta.data$TLS=="TLS",]%>% filter(! sample%in%c("S12","S13","S16","S17","S4","S5","S9")) 
color_dict=list(bg = "#D3D3D3",TLS="#696969" ,sample=c("#9FBBD5","#B8474D"),
sample2=c("#1DBDC6","#E4C671"),
sample3=c("#6EBDB7","#9C4084")
)
ggplot_list=list()

ggplot_list[["NFKB1"]]=ggplot(AUC_exp%>%{.[order(.$`NFKB1(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`NFKB1(+)`>0.068))+

geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$bg)+

geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",],aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


ggplot_list[["NFKB2"]]=ggplot(AUC_exp%>%{.[order(.$`NFKB2(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`NFKB2(+)`>0.06))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+

geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()



ggplot_list[["RUNX1"]]=ggplot(AUC_exp%>%{.[order(.$`RUNX1(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`RUNX1(+)`>0.08))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample2),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


ggplot_list[["NFATC1"]]=ggplot(AUC_exp%>%{.[order(.$`NFATC1(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`NFATC1(+)`>0.07))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample3),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()

ggplot_list[["NFATC2"]]=ggplot(AUC_exp%>%{.[order(.$`NFATC2(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`NFATC2(+)`>0.11))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample3),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


ggplot_list[["STAT5A"]]=ggplot(AUC_exp%>%{.[order(.$`STAT5A(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`STAT5A(+)`>0.06))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


ggplot_list[["STAT5B"]]=ggplot(AUC_exp%>%{.[order(.$`STAT5B(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`STAT5B(+)`>0.09))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()



ggplot_list[["FOS"]]=ggplot(AUC_exp%>%{.[order(.$`FOS(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`FOS(+)`>0.06))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample2),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()





ggplot_list[["BCL6"]]=ggplot(AUC_exp%>%{.[order(.$`BCL6(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`BCL6(+)`>0.056))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample3),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()




ggplot_list[["KLF2"]]=ggplot(AUC_exp%>%{.[order(.$`KLF2(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`KLF2(+)`>0.07))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


ggplot_list[["BCL11B"]]=ggplot(AUC_exp%>%{.[order(.$`BCL11B(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`BCL11B(+)`>0.10))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


ggplot_list[["LEF1"]]=ggplot(AUC_exp%>%{.[order(.$`LEF1(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`LEF1(+)`>0.07))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample3),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


ggplot_list[["RUNX3"]]=ggplot(AUC_exp%>%{.[order(.$`RUNX3(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`RUNX3(+)`>0.12))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()

ggplot_list[["EP300"]]=ggplot(AUC_exp%>%{.[order(.$`EP300(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`EP300(+)`>0.08))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()

ggplot_list[["XBP1"]]=ggplot(AUC_exp%>%{.[order(.$`XBP1(+)`),]},aes(x=x_xlide_px,y=y_xlide_px,color=`XBP1(+)`>0.10))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
scale_colour_manual(values =c( color_dict$sample),
            na.value = "white")+ggplot2::coord_fixed()+ ggplot2::theme_classic()


for(ggplot_names  in names(ggplot_list)){
  ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/转录因子原位/",ggplot_names,".png"),ggplot_list[[ggplot_names]]+ggplot2::labs(title=ggplot_names),
  ,height=15*1.5,width=20*1.5)
    ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/转录因子原位/",ggplot_names,".pdf"),ggplot_list[[ggplot_names]]+ggplot2::labs(title=ggplot_names),
  ,height=15*1.5,width=20*1.5)
}



# AUCELL 可视化 

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
c2 <- clusterProfiler::read.gmt("/home/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/output/fig2/ssgsea/pathway_gmt.gmt") 
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
library(plyr )

aucell_test=ddply(auc_data_long, .(gene_set), function(x) { wilcox.test(expression~TLS, data=x, paired=FALSE)$p.value })
aucell_test$FDR= p.adjust(aucell_test$V1)
auc_data_long_mean=auc_data_long%>%group_by( TLS ,gene_set ) %>% dplyr::summarize(expression=mean(expression))
auc_data_wider_mean=auc_data_long_mean%>% tidyr::pivot_wider(names_from=TLS,values_from=expression) %>% mutate(diff= TLS/nTLS)
auc_data_wider_mean_add_FDR= full_join(auc_data_wider_mean, aucell_test)
auc_data_wider_mean_add_FDR=auc_data_wider_mean_add_FDR[auc_data_wider_mean_add_FDR$FDR<0.01,]
auc_data_wider_mean_add_FDR=auc_data_wider_mean_add_FDR[auc_data_wider_mean_add_FDR$diff>1.1,]


TLS=meta.data$TLS[match(colnames(getAUC(cells_AUC)),rownames(meta.data))]
x=meta.data$x_xlide_px[match(colnames(getAUC(cells_AUC)),rownames(meta.data))]
y=meta.data$y_xlide_px[match(colnames(getAUC(cells_AUC)),rownames(meta.data))]
y=meta.data$y_xlide_px[match(colnames(getAUC(cells_AUC)),rownames(meta.data))]
meta.data=meta.data%>% filter(! sample_id%in%c("S12","S13","S16","S17","S4","S5","S9")) 
sample_id=meta.data$sample_id[match(colnames(getAUC(cells_AUC)),rownames(meta.data))]

auc_data=getAUC(cells_AUC)%>%data.frame %>%t %>%data.frame(TLS,x,y ,sample_id,.)%>% na.omit %>%tibble::rownames_to_column("cell")

ggplot_list[["BIOCARTA_TCR_PATHWAY"]]=ggplot(auc_data%>%{.[order(.$BIOCARTA_TCR_PATHWAY),]},aes(x=x,y=y,color=BIOCARTA_TCR_PATHWAY))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
viridis::scale_color_viridis(option="E",limits = c(quantile(auc_data$BIOCARTA_TCR_PATHWAY,0.05), quantile(auc_data$BIOCARTA_TCR_PATHWAY,0.95))) +
ggplot2::coord_fixed()+ ggplot2::theme_classic()
ggplot_names="BIOCARTA_TCR_PATHWAY"
# 相关的正向通路和对应的功能
naive_t_cell_pathways <- list(
  BIOCARTA_TCR_PATHWAY = "T 细胞受体（TCR）信号通路，T 细胞活化的关键步骤，激活下游信号分子，导致 T 细胞的增殖、分化和激活。",
  BIOCARTA_IL2RB_PATHWAY = "IL-2 受体信号通路，促进 T 细胞增殖和存活。",
  BIOCARTA_THELPER_PATHWAY = "辅助性 T 细胞（Th 细胞）分化通路，指导 naive T 细胞分化为不同 Th 细胞亚型。",
  BIOCARTA_CD40_PATHWAY = "CD40 信号通路，调控 T 细胞与 B 细胞、树突状细胞等抗原呈递细胞之间的相互作用。",
  BIOCARTA_IL7_PATHWAY = "IL-7 信号通路，维持 naive T 细胞存活，防止 T 细胞的凋亡。",
  BIOCARTA_TOB1_PATHWAY = "TOB1 是 T 细胞活化的负调控因子，有助于维持免疫反应的平衡。",
  KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY = "T 细胞受体信号通路，通过识别抗原激活 T 细胞，促进其增殖和分化。",
  GSE13738_RESTING_VS_TCR_ACTIVATED_CD4_TCELL_UP = "展示休眠 CD4+ T 细胞与 TCR 激活后的基因表达差异，体现 T 细胞的激活。",
  GSE22886_NAIVE_CD4_TCELL_VS_12H_ACT_TH1_UP = "展示 naive CD4+ T 细胞与 12 小时激活后 Th1 细胞的基因表达差异，体现 Th1 细胞的分化过程。",
  GSE22886_NAIVE_TCELL_VS_DC_UP = "展示 naive T 细胞与树突状细胞之间的基因表达差异，反映 T 细胞与抗原呈递细胞的相互作用。",
  GSE22886_NAIVE_TCELL_VS_NKCELL_UP = "展示 naive T 细胞与自然杀伤细胞（NK 细胞）之间的基因表达差异，涉及 T 细胞激活后的功能变化。"
)

# 打印路径及其功能
names(naive_t_cell_pathways)

for(ggplot_names  in names(naive_t_cell_pathways)){
ggplot_list[[ggplot_names]]=ggplot(auc_data%>%{.[order(.[,ggplot_names]),]},aes_string(x="x",y="y",color=ggplot_names))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
viridis::scale_color_viridis(option="E",limits = c(quantile(auc_data[,ggplot_names],0.05), quantile(auc_data[,ggplot_names],0.95))) +
ggplot2::coord_fixed()+ ggplot2::theme_classic()


  ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/AUCell_原位/",ggplot_names,".png"),ggplot_list[[ggplot_names]]+ggplot2::labs(title=ggplot_names),
  ,height=15*1.5,width=20*1.5)
    ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/AUCell_原位/",ggplot_names,".pdf"),ggplot_list[[ggplot_names]]+ggplot2::labs(title=ggplot_names),
  ,height=15*1.5,width=20*1.5)
}
pathway_exp=list("BIOCARTA_TCR_PATHWAY"="A",
"GSE13738_RESTING_VS_TCR_ACTIVATED_CD4_TCELL_UP"="G",
"GSE22886_NAIVE_CD4_TCELL_VS_12H_ACT_TH1_UP"="C",
"GSE22886_NAIVE_TCELL_VS_NKCELL_UP"="D",
"KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"="E"
)



for(ggplot_names  in names(pathway_exp)){
ggplot_list[[ggplot_names]]=ggplot(auc_data%>%{.[order(.[,ggplot_names]),]},aes_string(x="x",y="y",color=ggplot_names))+
geom_point(data=adata_meta.data  ,aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3")+
geom_point(data=TLS_meta.data[TLS_meta.data$TLS=="TLS",] ,aes(x=x_xlide_px,y=y_xlide_px),color=color_dict$TLS)+
geom_point()+
viridis::scale_color_viridis(option=pathway_exp[[ggplot_names]],limits = c(quantile(auc_data[,ggplot_names],0.05), quantile(auc_data[,ggplot_names],0.95))) +
ggplot2::coord_fixed()+ ggplot2::theme_classic()
  ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/AUCell_原位/",ggplot_names,".png"),ggplot_list[[ggplot_names]]+ggplot2::labs(title=ggplot_names),
  ,height=15*1.5,width=20*1.5)
    ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/AUCell_原位/",ggplot_names,".pdf"),ggplot_list[[ggplot_names]]+ggplot2::labs(title=ggplot_names),
  ,height=15*1.5,width=20*1.5)
}