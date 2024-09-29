celltype_color_dict =c('Cancer'= '#1f77b4',
 'Other T cell'= '#ff7f0e',
 'Fibroblast'= '#2ca02c',
 'Myeloid'='#8c564b',
 'Endothelial'= '#9467bd',
 'CD74 B cells'=  "#db4d57",
 'naive_CD74 B cells'= '#e377c2',
 'DC'= '#7f7f7f',
 'Naive T'= '#bcbd22',
 "IGHV1 Plasma cells"=  "#689eac",
  "IGHM Plasma cells" = '#ADB467',
 "IGLV4 Plasma cells"  ='#CECF84' ,
  "IGKV4 Plasma cells"  = '#5DA05A',
  "Other"="#D3D3D3"
  )

color_dict=c("CD74 B cells"="#FA7F73","IGKV4 Plasma cells"="#8DD1C6","IGLV4 Plasma cells"="#B1DD6D",
            "IGHV1 Plasma cells"="#BBB8D9","IGHM Plasma cells"="#FCB264",
            'Naive T'='#1f77b4'
            )
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
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")

B_meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
adata=read_h5ad("output/data/adata.h5ad")
Tcell=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
#x=t(Tcell$X)

a.meta.data=read_h5ad("output/data/dbscan.h5ad")$obs

X=t(adata$X)
adata_meta.data=adata$obs
adata_meta.data$cell_type=as.character(adata_meta.data$cell_type)


meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  
adata_meta.data$TLS="nTLS"
adata_meta.data$TLS[rownames(adata_meta.data)%in% (rownames(meta.data)[meta.data$TLS=="TLS"])]="TLS"


adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))]=meta.data$cell_type_sub %>% as.character

Tcell=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
Tcell_meta.data<- Tcell$obs
Tcell_meta.data$cell_subtype=Tcell_meta.data$cell_subtype %>% as.character
# TLS 
a.meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
adata_meta.data$TLS[match(rownames(a.meta.data),rownames(adata_meta.data))]=ifelse(a.meta.data$dbscan_labels!=-1 ,"TLS","nTLS")  
adata_meta.data$TLS[match(rownames(a.meta.data),rownames(adata_meta.data))]=ifelse(a.meta.data$dbscan_labels!=-1 ,"TLS","nTLS")  
adata_meta.data$TLS="nTLS"
adata_meta.data$TLS[rownames(adata_meta.data)%in% (rownames(meta.data)[meta.data$TLS=="TLS"])]="TLS"



adata_meta.data$cell_type[match(rownames(Tcell_meta.data),rownames(adata_meta.data))]=Tcell_meta.data$cell_subtype# %>% as.character
adata_meta.data$cell_type[!adata_meta.data$cell_type%in%c("Naive T","CD74 B cells","IGKV4 Plasma cells","IGLV4 Plasma cells","IGHV1 Plasma cells","IGHM Plasma cells")]="Other"
adata_meta.data=adata_meta.data %>% filter(! sample_id %in%c("S12","S13","S16","S17","S4","S5","S9")) 
a=ggplot()+ 
geom_point(data=adata_meta.data[adata_meta.data$cell_type=="Other",],mapping=aes(x=x_xlide_px,y=y_xlide_px),color="#D3D3D3",size=0.1)+
geom_point(data=adata_meta.data[adata_meta.data$TLS=="TLS",],mapping=aes(x=x_xlide_px,y=y_xlide_px),color="#696969",size=0.2)+
#geom_point(data=adata_meta.data[adata_meta.data$TLS=="TLS",],mapping=aes(x=x_xlide_px,y=y_xlide_px),color="#696969",size=0.2	)+
geom_point(data=adata_meta.data[adata_meta.data$cell_type!="Other",],mapping=aes(x=x_xlide_px,y=y_xlide_px,color=cell_type),size=0.1)+
scale_color_manual(values=celltype_color_dict)+
theme_classic()+ggplot2::coord_fixed()
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/细胞通讯原位可视化/细胞通讯原位可视化.png",a, width = 15*2, height =20*2 )
# 使用 ggplot2 绘图并明确顺序
adata_meta.data$plot_color <- ifelse(adata_meta.data$TLS == "TLS", "#696969",
                              ifelse(adata_meta.data$cell_type != "Other", 
                                     celltype_color_dict[adata_meta.data$cell_type], 
                                     NA))
adata_meta.data %>% case_when(
adata_meta.data$cell_type == "Other" ~ "Other",

adata_meta.data$cell_type != "Other" ~ adata_meta.data$cell_type,


)
adata_meta.data$plot_color="Other"
adata_meta.data$plot_color[adata_meta.data$TLS=="TLS"]="TLS"
adata_meta.data$plot_color[adata_meta.data$cell_type != "Other"]=adata_meta.data$cell_type[adata_meta.data$cell_type != "Other"]



a = ggplot() 
 a=a+ geom_point(data = adata_meta.data[adata_meta.data$cell_type == "Other", ], mapping = aes(x = x_xlide_px, y = y_xlide_px), color = "#D3D3D3", 
 size = 0.1, alpha = 1) 
  
  a=a+ geom_point(data = adata_meta.data[adata_meta.data$TLS == "TLS", ], mapping = aes(x = x_xlide_px, y = y_xlide_px),  color = "#A9A9A9", size = 2) 
  
  a=a+ geom_point(data = adata_meta.data[adata_meta.data$cell_type != "Other", ], mapping = aes(x = x_xlide_px, y = y_xlide_px, color = cell_type), size = 0.1) 
  
  a=a+scale_color_manual(values = color_dict) 
  a=a+theme_classic() 
  a=a+ ggplot2::coord_fixed()

# 保存图像
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/细胞通讯原位可视化/细胞通讯原位可视化.png", 
       a, width = 15 * 1.5, height = 20 * 1.5)
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/细胞通讯原位可视化/细胞通讯原位可视化.pdf", 
       a, width = 15 * 1.5, height = 20 * 1.5)



library()
t_cell_related_interactions <- c(
  "AREG_EGFR",
  "CD40LG_ITGA5_ITGB1",
  "CD40LG_ITGAM_ITGB2",
  "CD86_CD28",
  "CD86_CTLA4",
  "HLA-A_CD8A",
  "HLA-A_CD8B",
  "HLA-B_CD8A",
  "HLA-B_CD8B",
  "HLA-C_CD8A",
  "HLA-C_CD8B",
  "HLA-DMA_CD4",
  "HLA-DMB_CD4",
  "HLA-DOA_CD4",
  "HLA-DOB_CD4",
  "HLA-DPA1_CD4",
  "HLA-DPB1_CD4",
  "HLA-DQA1_CD4",
  "HLA-DQA2_CD4",
  "HLA-DQB1_CD4",
  "HLA-DRA_CD4",
  "HLA-DRB1_CD4",
  "HLA-DRB5_CD4",
  "IL2_IL2RA_IL2RB",
  "IL2_IL2RB_IL2RG",
  "IL7_IL7R_IL2RG",
  "LTA_TNFRSF1A",
  "TNF_TNFRSF1A"
)
split_sp_raw <- list()  # initialize an empty list to store plots
draw_empty_ggplot <- function(title) {
  ggplot() + 
    geom_blank() +  # 保持图表空白
    labs(title = title) +  # 添加标题
    theme_void() +  # 移除所有背景、坐标轴、标签等
    theme(plot.title = element_text(hjust = 0.5))  # 标题居中
}
for(sample in names(cellchat_run)){
    plot_title <- paste(sample, "no form B cell to navie T", sep=" - ")
    split_sp_raw[[sample]] <- draw_empty_ggplot(plot_title)
}


source("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/细胞通讯画图定制版本.R")

for(signals in ligand_receptor_pairs ){
    system(paste( "echo ",signals ,">>run.sh"))
split_sp=split_sp_raw
sub_line=c("Naive T","CD74 B cells","IGHM Plasma cells","IGHV1 Plasma cells","IGKV4 Plasma cells", "IGLV4 Plasma cells")
for(sample in c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])
print(paste(sample,signals,sep="-"))
try({split_sp[[paste(signals,sample,sep="-")]]=netVisual_aggregate_check(
                    object=cellchat_run[[sample]],chech.name=signals,
                    alpha.edge=1,
                    #signaling.name =paste(sample,signaling,sep="-"),
                    layout = "spatial",
                    sources.use = c("CD74 B cells","IGHM Plasma cells","IGHV1 Plasma cells","IGKV4 Plasma cells", "IGLV4 Plasma cells"),
                    targets.use="Naive T",
 edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.09, vertex.label.cex = 0,subcell=sub_line
                    )+
 ggplot2::labs(title=paste(sample,signals,sep="-"))+scale_color_manual(values=color_dict)
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/细胞通讯原位可视化/",signals,".png"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE,dpi=300)
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/细胞通讯原位可视化/",signals,".pdf"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
system(paste( "echo ",signals ,"is ok >>run.sh"))
}

# 转录因子 原位可视化 
