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

adata=read_h5ad("output/data/adata.h5ad")
X=t(adata$X)
TLS_meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
TLS_meta.data$TLS=ifelse(TLS_meta.data$dbscan_labels!=-1 ,"TLS","nTLS")  
B_meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
TLS_meta.data$cell_type=as.character(TLS_meta.data$cell_type)
TLS_meta.data$cell_type[match(rownames(B_meta.data),rownames(TLS_meta.data))]=#[meta.data$cell_type_sub=="CD74 B cells"]
as.character(B_meta.data$cell_type_sub)
T_meta.data=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")$obs
TLS_meta.data$cell_type[match(rownames(T_meta.data),rownames(TLS_meta.data))]=#[meta.data$cell_type_sub=="CD74 B cells"]
as.character(T_meta.data$cell_subtype)

#rownames(adata_meta.data) [adata_meta.data$cell_type%in%c("TLS_Naive T","nTLS_Naive T")]%>% data.frame %>% write.csv("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Naive_T_cell_names.csv")

TLS_meta.data=TLS_meta.data[TLS_meta.data$TLS %in%"TLS",]
X=X[,rownames(TLS_meta.data)]
#Seurat_data=CreateSeuratObject(X,meta.data =TLS_meta.data )


CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB
) # use Secreted Signaling

conversion.factor = 0.18
TLS_meta.data=TLS_meta.data[!is.na(TLS_meta.data$cell_type),]
# 创建进度条
options(future.globals.maxSize = 50 * 1024^3)
future::plan("multisession", workers =20)  # Reduce the number of workers
library(progressr)
handlers(global = TRUE)
handlers("progress")
with_progress({
p <- progressr::progressor(steps = length(TLS_meta.data$sample%>%unique))
# 使用 future_apply 并传递进度条
 cellchat_run=future.apply::future_lapply(TLS_meta.data$sample%>%unique,\(.x){

.x=as.character(.x)
meta=TLS_meta.data[TLS_meta.data$sample==.x,]
data.input=X[,rownames(meta)]
spatial.locs=meta[,c("x_slide_mm","y_slide_mm")]
d = computeCellDistance(spatial.locs)
spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
#meta$samples=meta$group
meta$cell_type=as.character(meta$cell_type)
cellchat <- createCellChat(object =data.input, meta = meta,group.by = "cell_type",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat@DB=CellChatDB.use
# set the used database in the object
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
#future::plan("multisession", workers =2)  # Reduce the number of workers


cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 422
#>

# Project the data
cellchat <- smoothData(cellchat, adj = PPI.human)
#cellchat <- projectData(cellchat, PPI.human)
#cellchat2 <- computeCommunProb(cellchat, type = "truncatedMean")
cellchat2 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE,raw.use = FALSE,
                              interaction.range = 250,
                              scale.distance =1000,k.min = 5, nboot = 200,
                              contact.dependent = TRUE, contact.range = 100)
gc()

cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
#saveRDS(cellchat2,paste0("/home/wangdongbin/项目/2024/hongkong_chinese_smi/hongkong_chinese_smi/output/cell_chat/ZFP36/",i,".RDS"))
return(cellchat2)
    }, future.seed = TRUE)

})
names(cellchat_run)=TLS_meta.data$sample%>%unique
saveRDS(cellchat_run,"/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/cellchat_k.min_5.RDS")  

#cellchat=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cellchat.RDS")  
#names(cellchat_run)=TLS_meta.data$sample%>%unique
cellchat_run=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/cellchat_k.min_5.RDS")  
#table(meta.data$sample_id,meta.data$TLS)


#cellchat_run=cellchat_run[-which(names(cellchat_run)%in%c("S12","S13","S16","S17","S4","S5","S9"))]
a=list()
for (i in names(cellchat_run)) {
  # Extract the CellChat object for the current sample
cellchat_obj <- cellchat_run[[i]]
a[[i]]=cellchat_obj@net$count#[rownames(cellchat_run[["S7"]]@net$count),rownames(cellchat_run[["S7"]]@net$count)]
  # You can add more analyses or visualizations here
  # For example, netVisual_bubble or netAnalysis_contribution
}

a <- list()

# Iterate over each CellChat object in cellchat_run
for (i in names(cellchat_run)) {
  # Extract the CellChat object for the current sample
  cellchat_obj <- cellchat_run[[i]]
  # Store the count matrix in the list a
 count<- cellchat_obj@net$count
 a[[i]]=count %>% data.frame(check.names=F)%>%  tibble::rownames_to_column("form") %>% tidyr::pivot_longer(!form)
 a[[i]]$sample=i
}


cellchat_run_data_long= a%>% purrr::reduce(rbind) %>% group_by(form , name )%>% dplyr::summarise(exp=sum(value)) 
table(cellchat_run_data_long$form)
cell_chat_count=cellchat_run_data_long%>% 
    tidyr::pivot_wider(names_from=name,values_from=exp) %>% 
    data.frame(check.names=F)


rownames(cell_chat_count)=cell_chat_count[,1]
cell_rank=c("DC"            ,  "CD74 B cells"                 ,
   "IGHM Plasma cells"  ,"IGHV1 Plasma cells",
"IGKV4 Plasma cells" ,"IGLV4 Plasma cells" ,         
"Naive T"           , "Proliferative T"  ,  "Th1"     ,          
"Th17"          ,     "Treg"             ,  "γδ T", "CTL"    ,"NK"       )
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/net.pdf",width=12, height=12)
netVisual_circle(cell_chat_count[cell_rank,cell_rank] %>% as.matrix, weight.scale = T, label.edge= F)
dev.off()

merged_cellchat <- mergeCellChat(cellchat_run, add.names = names(cellchat_run),
                                 merge.data = FALSE,
    cell.prefix = T)


bind_net=lapply(names(merged_cellchat@net),\(i){

 .x= merged_cellchat@net[[i]]
 pval_data <- as.data.frame(.x$pval)
# 筛选 p
pval_data=pval_data%>%tibble::rownames_to_column("from")%>% tidyr::pivot_longer(!from)
pval_data=pval_data[pval_data$value<0.05,]
  colnames(pval_data)[3]="p"
#
# 筛选 prob
   prob_data <- as.data.frame(.x$prob)
prob_data=prob_data%>%tibble::rownames_to_column("from")%>% tidyr::pivot_longer(!from)
   colnames(prob_data)[3]="prob"
try({ 
    d=left_join(pval_data,prob_data)
  d$group=i
  return(d)})
})
bind_net=do.call(rbind,bind_net)
bind_net_count=bind_net %>% #[stringr::str_detect(bind_net$from ,"ETS_CD74 B cells")|stringr::str_detect(bind_net$from ,"TLS_CD74 B cells"),]
                        filter(p<0.05)%>% dplyr::mutate(to=stringr::str_remove(name,"\\..*"),path=stringr::str_remove(name,".*\\.")) %>%
                        group_by(from ,to,path ) %>% 
                        summarise(count=n(),prob=mean(prob))
bind_net_count=bind_net_count %>% filter(from %in% c("CD74 B cells","IGHM Plasma cells","IGHV1 Plasma cells","IGKV4 Plasma cells","IGLV4 Plasma cells") & 
to %in% c("Naive T"))


bind_net_count_to_wider=bind_net_count[,-5] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = count, values_fill = list(count = 0))



#get_max_count=bind_net_count%>%group_by(path)%>%summarise(max_count=max(count)>5)
#get_max_count=get_max_count$path[get_max_count$max_count]
#bind_net_count=bind_net_count[bind_net_count$path %in%No_rm_path, ]
bind_net_count$prob[bind_net_count$prob == 0] <- NA
bind_net_count$prob.original <- bind_net_count$prob
bind_net_count $prob <- -1/log(bind_net_count$prob)
min.quantile = 0 ;max.quantile = 0.99

min.cutoff <- quantile(bind_net_count$prob, min.quantile, na.rm = T)
max.cutoff <- quantile(bind_net_count$prob, max.quantile, na.rm = T)
bind_net_count$prob[bind_net_count$prob < min.cutoff] <- min.cutoff
bind_net_count$prob[bind_net_count$prob > max.cutoff] <- max.cutoff
bind_net_count=bind_net_count%>%na.omit
# 标准化
bind_net_count$rate=#[bind_net_count$from=="ETS_CD74 B cells"]
bind_net_count$count#[bind_net_count$from=="ETS_CD74 B cells"]/18
bind_net_count$rate=#[bind_net_count$from=="TLS_CD74 B cells"]
bind_net_count$count#[bind_net_count$from=="TLS_CD74 B cells"]/13


bind_net_count_to_wider=bind_net_count[,c(1,2,3,7)] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = rate, values_fill = list(rate = 0))


#adata.write_h5ad("/data/project/wangdongbin/data/2024//2024_7_13_E007_TLS/output/data/adata_B_cell_type_dbscan_TLS.h5ad")
ligand_receptor_pairs <- c(
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
  "HLA-E_CD8A",
  "HLA-E_CD8B",
  "HLA-F_CD8A",
  "HLA-F_CD8B",
  "HLA-G_CD8A",
  "HLA-G_CD8B",
  "CD40LG_CD40",
  "CD80_CD28",
  "CD80_CTLA4",
  "CD86_CD28",
  "CD86_CTLA4",
  "MIF_CD74",
  "MIF_CD44",
  "ICAM1_ITGAL",
  "ICAM1_ITGAM_ITGB2",
  "ICAM1_ITGAX_ITGB2",
  "JAG1_NOTCH1",
  "JAG1_NOTCH2",
  "JAG1_NOTCH3",
  "IL2_IL2RA_IL2RB",
  "IL2_IL2RB_IL2RG",
  "IFNG_IFNGR1_IFNGR2"
)


bind_net_count$to=stringr::str_remove(bind_net_count$to," Plasma cells")
combos <- expand.grid(from = bind_net_count$from %>% unique, to = bind_net_count$to %>% unique,  stringsAsFactors = FALSE)
a=ggplot(bind_net_count#[bind_net_count$from%in%c("CD74 B cells","CD74 B cells"),]
, #& bind_net_count$path%in% ligand_receptor_pairs
         aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/cell_chat_netVisual_bubble_from.pdf",
       a,height=10,width=7,limitsize = FALSE)

ggsave("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/cell_chat_netVisual_bubble_from.pdf",
       a,height=50,width=5,limitsize = FALSE)

ligand_receptor_pairs="TNF_TNFRSF1B
LTA_TNFRSF1B
LTA_TNFRSF14
IL16_CD4
HLA-F_CD8B
HLA-E_KLRK1
HLA-DRA_CD4
HLA-DQB1_CD4
HLA-DQA2_CD4
HLA-DMB_CD4
HLA-B_CD8A
HLA-A_CD8A
CXCL12_CXCR4
CCL19_CCR7
IL6_IL6R_IL6ST
COL4A4_CD44
CLEC2D_KLRB1
CD55_ADGRE5
CD40LG_ITGAM_ITGB2
CCL21_CCR7
" %>% stringr::str_split("\n")%>% unlist
a=ggplot(bind_net_count[bind_net_count$from%in%c("TLS_CD74 B cells","nTLS_CD74 B cells")&
 bind_net_count$path%in% ligand_receptor_pairs ,],
         aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/细胞通讯/cell_chat_netVisual_bubble_from.pdf",
       a,height=6,width=4,limitsize = FALSE)
# 与T细胞增殖、分化、激活和抗原呈递相关的受体-配体对
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

a=ggplot(bind_net_count[bind_net_count$path%in% t_cell_related_interactions,]
, #& bind_net_count$path%in% ligand_receptor_pairs
         aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/cell_chat_netVisual_bubble_from.pdf",
       a,height=7.5,width=5,limitsize = FALSE)

ggsave("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/cell_chat_netVisual_bubble_from.pdf",
       a,height=50,width=5,limitsize = FALSE)


library(ggplot2)
library(ggnewscale)





for(signals in t_cell_related_interactions){
    
    
try({
system(paste( "echo ",signals ,">>run.sh"))
split_st=stringr::str_split(signals,"_")[[1]]
ggplot_obj=ggplot()+
geom_point(data=gene_exp,aes_string(x="x",y="y"),color="#DCDCDC")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" ,],aes_string(x="x",y="y"),color="#A9A9A9")+


geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" & gene_exp[,split_st[1]]>0,],
aes_string(x="x",y="y",color=split_st[1]),size=5)+
  # 第一种颜色渐变
  scale_color_gradient(low =  "#FBDFE2", high ="#B83945", name = split_st[1]) +
  # 添加新的颜色尺度
  ggnewscale::new_scale("color") +
  # 第二种颜色渐变
  geom_point(data=gene_exp[gene_exp$cell_type %in% "NK & T"& gene_exp[,split_st[2]]>0,], aes_string(x="x",y="y",color=split_st[2]),size=5) +
  scale_color_gradient(low = "#CFE7C4", high = "#4F845C", name = split_st[2]) +
  theme_minimal()+
 ggplot2::labs(title=paste(signals))+
 ggplot2::coord_fixed()
#+scale_color_manual(values=celltype_color_dict)
  if(length(split_st)>2){
ggplot_obj=ggplot_obj+    geom_point(data=gene_exp[gene_exp$cell_type %in% "NK & T"& gene_exp[,split_st[3]]>0,],aes_string(x="x",y="y",color=split_st[3]),size=5)
  scale_color_gradient(low = "#C7DFF0", high = "#377483", name = split_st[3]) 
  }
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/原位/",signals,".pdf"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
})
}
##

#  基因原位 
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

adata=read_h5ad("output/data/adata.h5ad")
X=t(adata$X)
TLS_meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
TLS_meta.data$TLS=ifelse(TLS_meta.data$dbscan_labels!=-1 ,"TLS","nTLS")  
B_meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
TLS_meta.data$cell_type=as.character(TLS_meta.data$cell_type)
TLS_meta.data$cell_type[match(rownames(B_meta.data),rownames(TLS_meta.data))]=#[meta.data$cell_type_sub=="CD74 B cells"]
as.character(B_meta.data$cell_type_sub)
T_meta.data=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")$obs
TLS_meta.data$cell_type[match(rownames(T_meta.data),rownames(TLS_meta.data))]=#[meta.data$cell_type_sub=="CD74 B cells"]
as.character(T_meta.data$cell_subtype)

#rownames(adata_meta.data) [adata_meta.data$cell_type%in%c("TLS_Naive T","nTLS_Naive T")]%>% data.frame %>% write.csv("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Naive_T_cell_names.csv")

TLS_meta.data=TLS_meta.data[TLS_meta.data$TLS %in%"TLS",]
TLS_adata=read_h5ad("output/data/dbscan.h5ad")
TLS_data=t(TLS_adata$X)

gene_exp=data.frame(x=TLS_meta.data$x_slide_mm,
y=TLS_meta.data$y_slide_mm,
cell_type=TLS_meta.data$cell_type,
TLS_data[,rownames(TLS_meta.data)]%>%t)

   

ligand_receptor_pairs <- c(
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
for(signals in ligand_receptor_pairs%>% make.names ){
try({
  cat(signals)
  system(paste( "echo ",signals ,">>run.sh"))
split_st=stringr::str_split(signals,"_")[[1]]
ggplot_obj=ggplot()+
geom_point(data=gene_exp,aes_string(x="x",y="y"),color="#DCDCDC")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "CD74 B cells" ,],aes_string(x="x",y="y"),color="#A9A9A9")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "CD74 B cells" & gene_exp[,split_st[1]]>0,]%>%.[order(.[,split_st[1]]),],
aes_string(x="x",y="y",color=split_st[1]),size=5)+
  # 第一种颜色渐变
scale_color_gradient(low =  "#FBDFE2", high ="#B83945", name = split_st[1]) +theme_classic()+ggplot2::coord_fixed()
  # 添加新的颜色尺度
  # 第二种颜色渐变
if(length(split_st)==2){
ggplot_obj=ggplot_obj+ ggnewscale::new_scale("color") + 
geom_point(data=gene_exp[gene_exp$cell_type %in% "Naive T"& gene_exp[,split_st[2]]>0,]%>%.[order(.[,split_st[2]]),], 
aes_string(x="x",y="y",color=split_st[2]),size=5) +
  scale_color_gradient(low = "#CFE7C4", high = "#4F845C", name = split_st[2]) +
 ggplot2::labs(title=signals)   
  }
#+scale_color_manual(values=celltype_color_dict)
  if(length(split_st)>2){
     gene_exp$ bing_gene_exp=0
 gene_exp$bing_gene_exp[gene_exp$cell_type %in% "Naive T"]= 
  gene_exp[gene_exp$cell_type %in% "Naive T",split_st[2]]+gene_exp[gene_exp$cell_type %in% "Naive T",split_st[3]]

ggplot_obj=ggplot_obj+      ggnewscale::new_scale("color") +
geom_point(data=gene_exp[gene_exp$cell_type %in% "Naive T"&  gene_exp$bing_gene_exp>0,]%>%.[order(.[,"bing_gene_exp"]),],aes_string(x="x",y="y",color="bing_gene_exp"),size=5)+
  scale_color_gradient(low = "#C7DFF0", high = "#377483", name =paste0(split_st[2],"_", split_st[3])) 
  }
#ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".pdf"),ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/原位/",signals,".png"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/TLS_cell_chat/原位/",signals,".pdf"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
})
}

