library(Seurat)
library(tidyverse)
library(magrittr)
library(monocle)
library(anndata)
library(dplyr)
#library(SeuratDisk)
library(ggplot2)
library(ggalluvial)
library("Matrix")
library(monocle)
plot_genes_in_pseudotime=monocle::plot_genes_in_pseudotime
reticulate::use_python('/data/project/wangdongbin/.conda/envs/mistyr/bin/python')


meta.data_map=read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/map_cells_to_space.h5ad")$X
cell_type=read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/map_cells_to_space.h5ad")$obs$cell_type_sub
cell_type=cell_type[meta.data_map%>%apply(2,which.max)]


colors20=c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

x=t(read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/sc_B_cell_type.h5ad")$X)
meta.data=read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/sc_B_cell_type.h5ad")$obs

data2=Seurat::CreateSeuratObject(x,meta.data=meta.data)

# 创建新的 CellDataSet 对象
cds <- newCellDataSet(cellData  = as.matrix(data2@assays$RNA$counts),
                      phenoData = new("AnnotatedDataFrame", data = data2@meta.data))

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
# 筛选表达至少在3个细胞中非零的基因
cds <- detectGenes(cds, min_expr = 0.1)
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, 
                                   mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id%>%unique)
cds_subset <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree'
,num_paths=0)
#cds_subset <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
# 规范化数据
# 计算分化时间
cds <- orderCells(cds_subset,reverse=F)
saveRDS(cds,"/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_cell_type_B_2.RDS")
cds <- orderCells(cds,reverse=F,root_state=2)
saveRDS(cds,"/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_cell_type_B_root_state_2.RDS")


color_dict=c("CD74 B cells"="#FA7F73","IGKV4 Plasma cells"="#8DD1C6","IGLV4 Plasma cells"="#B1DD6D",
            "IGHV1 Plasma cells"="#BBB8D9","IGHM Plasma cells"="#FCB264")


cds=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_cell_type_B_root_state_2.RDS")

#saveRDS(cds,"/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_cell_type_B.RDS")
#cds=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_cell_type_B_2.RDS")
cds$CD74_exp= cds@assayData$exprs["CD74",]
cds$sample2=cds$sample %>% as.character
cds$sample2[stringr::str_detect(cds$sample,"\\(") ] =stringr::str_remove_all(cds$sample,".*\\(|\\)")[stringr::str_detect(cds$sample,"\\(") ]

x3=data.table::fread("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/raw_data/单细胞测序患者临床数据汇总3.29.txt")
colors20=c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

cds$sample2= stringr::str_extract(cds$sample2,"\\d+")
x3$V1= stringr::str_extract(x3$V1,"\\d+")

x3$V1[!x3$V1%in%cds$sample2] %>%table

cds$第8版T分期= x3$第8版T分期[match(cds$sample2,x3$V1)]
#cds=cds[,!is.na(cds$第8版T分期)]
cell_type=read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/map_cells_to_space.h5ad")$obs$cell_type_sub
cell_type=cell_type[meta.data_map%>%apply(2,which.max)]
names(cell_type)=colnames(read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/map_cells_to_space.h5ad")$X)
cds$cell_subtype=cell_type[colnames(cds)]

x3$V1%in%cds$sample2
cds$T_state= x3$分期[match(cds$sample2,x3$V1)]
#cds$T_state[is.na(cds$T_state)]="x"
cds$T_state =dplyr::case_when(
  cds$T_state %in%c( "I","IIA") ~ "Early",
  cds$T_state %in%c( "IIIA","IIIB","IV","IIA") ~ "Late"
)
cds$T_state= "Early"




#subset(cds,!is.na(cds$T_state))
cds$T_state= x3$分期[match(cds$sample2,x3$V1)]
#cds$T_state[is.na(cds$T_state)]="x"
cds$T_state =dplyr::case_when(
  cds$T_state %in%c( "I","IIA") ~ "Early",
  cds$T_state %in%c( "IIIA","IIIB","IV","IIA") ~ "Late"
)
#source("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/code_path/拟时序/monocle_plot.R")
cds_filtered <- cds
cds_filtered@assayData$exprs=cds_filtered@assayData$exprs[,colnames(cds)[( !is.na(pData(cds)$T_state))]]
cds_filtered=cds_filtered[,colnames(cds_filtered@assayData$exprs)]


cds_filtered@assayData$exprs=cds_filtered@assayData$exprs[,which(colnames(cds_filtered@assayData$exprs)%in% rownames(pData(cds)))]
pData(cds_filtered)%>%dim
cds_filtered@reducedDimS=cds_filtered@reducedDimS[,colnames(cds_filtered@assayData$exprs)]
cds_filtered$T_state= x3$分期[match(cds_filtered$sample2,x3$V1)]
cds_filtered$M= x3$M[match(cds_filtered$sample2,x3$V1)] %>%factor( )
cds_filtered$N= x3$N[match(cds_filtered$sample2,x3$V1)]%>%factor()
cds_filtered$第8版T分期= x3$第8版T分期[match(cds_filtered$sample2,x3$V1)]%>%factor()
cds_filtered$TNM_state= x3$分期[match(cds_filtered$sample2,x3$V1)]%>%factor()
M=plot_cell_trajectory(cds_filtered, color_by = "M") +  scale_color_manual(values = colors20)
M$data$M=M$data$M%>% factor(levels=c(0,1))
M$data=M$data[order(M$data$M),]


plot_cell_trajectory(cds_filtered, color_by = "T_state") +  scale_fill_manual(values = colors20)
plot_cell_trajectory(cds_filtered, color_by = "M") +  scale_fill_manual(values = colors20)
plot_cell_trajectory(cds_filtered, color_by = "第8版T分期")   +  scale_fill_manual(values = colors20)
combined_plot <- gridExtra::marrangeGrob(list(
plot_cell_trajectory(cds_filtered, color_by = "cell_subtype")+ scale_color_manual(values=color_dict),
plot_cell_trajectory(cds_filtered, color_by = "Pseudotime")+     viridis::scale_color_viridis() ,
plot_cell_trajectory(cds_filtered, color_by = "CD74_exp")+   scale_colour_gradient(low = "white", high = "#FA7F73")
,#viridis::scale_fill_viridis(option = "C"),
plot_cell_trajectory(cds_filtered, color_by = "N") +  scale_color_manual(values = colors20),
#plot_cell_trajectory(cds_filtered, color_by = "M") +  scale_color_manual(values = colors20),
M,
plot_cell_trajectory(cds_filtered, color_by = "第8版T分期")   +  scale_color_manual(values = colors20),
plot_cell_trajectory(cds_filtered, color_by = "TNM_state")+ scale_color_manual(values=colors20),
plot_cell_trajectory(cds, color_by = "State")+ scale_color_manual(values=colors20)

)
, nrow =7, ncol =1)
# 显示合并后的图形
pdf("output/fig3/monicol/fig3d+e_trajectory_greg.pdf",width=8, height=26)
grid::grid.draw(combined_plot)
dev.off()

plot_cell_trajectory(cds_filtered, color_by = "cell_subtype")+ scale_color_manual(values=color_dict)+gplot2::facet_grid(~cell_subtype)

# 显示合并后的图形
pdf("output/fig3/monicol/fig3d+e_trajectory_greg_cell_subtype.pdf",width=8, height=26)
plot_cell_trajectory(cds_filtered, color_by = "cell_subtype")+ scale_color_manual(values=color_dict)+ggplot2::facet_grid(cell_subtype~.)
dev.off()





gene=c("CD74","LTB","HLA-DRA","HLA-A"  ,   "HLA-B" )
HSMM_filtered <- cds[gene,]
cds_subset <- HSMM_filtered[gene,]
cds@phenoData@data=cds@phenoData@data[colnames(cds@assayData$exprs),]
png("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_gene.png")
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type_sub",trend_formula="~ sm.ns(Pseudotime, df=3)")
dev.off()
exp=data.frame(as.matrix(exprs(cds_subset))%>%t,time=cds_subset$Pseudotime,cell_type_all_Treg=cds_subset$cell_subtype,check.names=F) 





genes="PTEN
RASA4
TBC1D10C
ZFP36L1
MED15
VPS29
SYVN1
IGHM
SH3BGRL3
CD74
HLA-DRA
TMSB10
HLA-DPA1
TMSB4X
RAC2
SLC25A6
HLA-F
CORO1A
LTB
HLA-DRB
CD37
PFN1
UBB
CYBA
SRSF2
B2M
PSAP
POU2AF1
H2BC5
IGHA1
CBFB
DUSP2
JCHAIN
SEC31A
SLC25A17
AQP7
ECI2
CCND2
PDPK1
ELP1
PSMD10
PVR
HSP90B1
NPFF
XBP1
CPNE1
EPB41L5
OGT
SNRNP70
TXNIP
IGLV1-40
COX1
COX2
MCL1
BTG2
IRF4
P4HB
HERPUD1
IGKC
RRBP1
SDC1
MZB1
PDK1
DERL3
SEC11C
PLPP5
COL1A1
IGHG1/2
IGHV1-69
IGKV4-1
GCDH
PIM2
IGLL1
IGLL5
IGLC1/2
H2AX
NEAT1
JUND
IGLV6-57
PHKG2
CTF1
SELENOS
TP53INP1

" %>% stringr::str_split("\n")%>%.[[1]]
proliferation_genes <- c(
  "CCND1",  # Cyclin D1
  "CDK4",   # Cyclin-dependent kinase 4
  "CDK6",   # Cyclin-dependent kinase 6
  "CCNE1",  # Cyclin E1
  "CDK2",   # Cyclin-dependent kinase 2
  "EGFR",   # Epidermal Growth Factor Receptor
  "VEGFA",  # Vascular Endothelial Growth Factor A
  "FGFR1",  # Fibroblast Growth Factor Receptor 1
  "KRAS",   # KRAS oncogene
  "PIK3CA", # Phosphoinositide-3-kinase catalytic subunit alpha
  "AKT1",   # AKT serine/threonine kinase 1
  "MYC",    # MYC proto-oncogene
  "FOS",    # FOS proto-oncogene, AP-1 transcription factor subunit
  "JUN",    # JUN proto-oncogene, AP-1 transcription factor subunit
  "TP53",   # Tumor protein p53
  "RB1",    # Retinoblastoma protein
  "PTEN",   # Phosphatase and tensin homolog
  "CDKN1A", # Cyclin-dependent kinase inhibitor 1A (p21)
  "CDKN2A", # Cyclin-dependent kinase inhibitor 2A (p16)
  "MDM2",   # MDM2 proto-oncogene
  "BCL2",   # BCL2 apoptosis regulator
  "E2F1",   # E2F transcription factor 1
  "MKI67"   # Marker of proliferation Ki-67
)

gene=c("CD74","LTB","HLA-DRA","HLA-A"  ,   "HLA-B" ,proliferation_genes,genes[genes%in% rownames(cds)]) %>% unique
HSMM_filtered <- cds[gene,]
cds_subset <- HSMM_filtered[gene,]


pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_gene.pdf")
plot_pseudotime_heatmap(cds_subset, 
                             num_cluster =2, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = viridis::viridis(256))
dev.off()

gene=c("HLA-DRA","HLA-A"  ,   "HLA-B" )
HSMM_filtered <- cds[gene,]
cds_subset <- HSMM_filtered[gene,]
source("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/code_path/拟时序/monocle_plot.R")
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/monocle_gene.pdf")
plot_genes_in_pseudotime(cds_subset, color_by = "cell_subtype")
dev.off()


exp=data.frame(as.matrix(exprs(cds_subset))%>%t,time=cds_subset$Pseudotime,cell_subtype=cds_subset$cell_subtype,check.names=F) 

for(i in gene){
a=ggplot(exp,aes_string(x="time",y=paste0("`",i,"`")
))+geom_point(aes(color=cell_subtype),size=2)+geom_smooth(method = "gam", formula = y ~ sm.ns(x, df = 3),color="#000000") + 
#stat_smooth(formula="y ~ sm.ns(x,  df=3)") +
#scale_y_log10()+ 
scale_color_manual(values=color_dict)+ 
monocle:::monocle_theme_opts()
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/gene/",i%>% make.names,".pdf"),a ,width = 10,height =3)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/monicol/gene/",i%>% make.names,".png"),a ,width = 10,height = 4)

}