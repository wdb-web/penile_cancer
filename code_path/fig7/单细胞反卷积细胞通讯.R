path=c(
"/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Endothelial_cell_map_cells_to_space_cell_to_sp.h5ad"  ,
"/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Myeloid_cell_map_cells_to_space_merged.h5ad"  ,
"/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/be_B_cell_map_cells_to_space_cell_to_sp.h5ad" ,               
"/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/be_T_cell_map_cells_to_space_cell_to_sp.h5ad" 
)
color_dict=c("CD74 B cells"="#FA7F73","IGKV4 Plasma cells"="#8DD1C6","IGLV4 Plasma cells"="#B1DD6D",
            "IGHV1 Plasma cells"="#BBB8D9","IGHM Plasma cells"="#FCB264",
            
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
data_sc_to_sp=lapply(path,\(read_path){
B_meta.data=read_h5ad(read_path)$X
data.frame(sc=rownames(B_meta.data),sp= colnames(B_meta.data)[apply(B_meta.data,1,which.max)])  
})

data_sc_to_sp_bind_data=do.call(rbind,data_sc_to_sp)

adata=read_h5ad("output/data/adata.h5ad")
adata_sc =read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")


B_meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
adata=read_h5ad("output/data/adata.h5ad")
Tcell=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
meta.data<- Tcell$obs
a.meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  
adata_meta.data=adata$obs
adata_meta.data$cell_type=as.character(adata_meta.data$cell_type)
meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub=="CD74 B cells"]=paste0(
  ifelse(meta.data$TLS=="Yes","TLS","nTLS"),
"_",meta.data$cell_type_sub)[meta.data$cell_type_sub=="CD74 B cells"]
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub!="CD74 B cells"]=paste0(
meta.data$cell_type_sub)[meta.data$cell_type_sub!="CD74 B cells"] %>% stringr::str_remove("Plasma cells")
Tcell=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
meta.data<- Tcell$obs
meta.data=meta.data[!is.na(meta.data$cell_subtype),]
meta.data$TLS="nTLS"
a.meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
meta.data$TLS[rownames(meta.data) %in%(rownames(a.meta.data)[a.meta.data$dbscan_labels!=-1]) ]="TLS"
table(meta.data$TLS,meta.data$cell_subtype)
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_subtype=="Naive T"]=paste0(
meta.data$TLS,
"_",meta.data$cell_subtype)[meta.data$cell_subtype=="Naive T"]
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_subtype!="Naive T"]=paste0(
meta.data$cell_subtype)[meta.data$cell_subtype!="Naive T"]



exp=t(adata_sc$X[data_sc_to_sp_bind_data$sc,])
colnames(exp)=data_sc_to_sp_bind_data$sc
meta_data=adata_meta.data[data_sc_to_sp_bind_data$sp,]#%>%t
rownames(meta_data)=data_sc_to_sp_bind_data$sc





Seurat_data=CreateSeuratObject(exp,meta.data =meta_data )




CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB
) # use Secreted Signaling

conversion.factor = 0.18

# 创建进度条
options(future.globals.maxSize = 50 * 1024^3)
future::plan("multisession", workers =20)  # Reduce the number of workers
library(progressr)
handlers(global = TRUE)
handlers("progress")
with_progress({
p <- progressr::progressor(steps = length(meta_data$sample%>%unique))
# 使用 future_apply 并传递进度条
 cellchat_run=future.apply::future_lapply(meta_data$sample%>%unique,\(.x){

.x=as.character(.x)
meta=meta_data[meta_data$sample==.x,]
data.input=exp[,rownames(meta)]
meta[,c("x_slide_mm","y_slide_mm")]
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
names(cellchat_run)=meta.data$sample_id%>%unique
saveRDS(cellchat_run,"/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/sc_to_sp_cellchat_k.min_5.RDS")  

cellchat_run=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/sc_to_sp_cellchat_k.min_5.RDS")
cellchat_run=cellchat_run[-which(names(cellchat_run)%in%c("S12","S13","S16","S17","S4","S5","S9"))]
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
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/code_path/fig7/scRNA/net.pdf",width=12, height=12)
netVisual_circle(cell_chat_count[,-1] %>% as.matrix, weight.scale = T, label.edge= F)
dev.off()

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

cellchat_run_data_long= a%>% purrr::reduce(rbind)%>%
  filter(! sample %in%c("S12","S13","S16","S17","S4","S5","S9")) %>%
  group_by(form , name )%>% dplyr::summarise(exp=mean(value))
table(cellchat_run_data_long$form)




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
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/scRNA/net.pdf",width=12, height=12)
netVisual_circle(cell_chat_count[c("TLS_CD74 B cells","nTLS_CD74 B cells",#"Fibroblast"   ,
"Myeloid" ,"NK" ,"Th17","Th1","CTL" ,
"TLS_Naive T","nTLS_Naive T"),c("TLS_CD74 B cells","nTLS_CD74 B cells",#"Fibroblast"   ,
"Myeloid" ,"NK" ,"Th17","Th1","CTL" ,

"TLS_Naive T","nTLS_Naive T")] %>% as.matrix, weight.scale = T, label.edge= F)
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
bind_net_count=bind_net_count %>% filter(from %in% c("TLS_CD74 B cells","nTLS_CD74 B cells") & to %in% c("TLS_Naive T","nTLS_Naive T","DC"))


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

bind_net_count$from=stringr::str_replace(bind_net_count$from,"No","NOTLS")%>%stringr::str_replace("Yes","TLS")
bind_net_count$to=stringr::str_replace(bind_net_count$to,"No","NOTLS")%>%stringr::str_replace("Yes","TLS")
bind_net_count_to_wider=bind_net_count[,c(1,2,3,7)] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = rate, values_fill = list(rate = 0))

bind_net_count$to=stringr::str_remove(bind_net_count$to," Plasma cells")
combos <- expand.grid(from = bind_net_count$from %>% unique, to = bind_net_count$to %>% unique,  stringsAsFactors = FALSE)
a=ggplot(bind_net_count[bind_net_count$from%in%c("TLS_CD74 B cells","nTLS_CD74 B cells") & bind_net_count$path%in% ligand_receptor_pairs,],
         aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/scRNA/cell_chat_netVisual_bubble_from.pdf",
       a,height=10,width=7,limitsize = FALSE)


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
a=ggplot(bind_net_count[bind_net_count$from%in%c("TLS_CD74 B cells","nTLS_CD74 B cells")& bind_net_count$path%in% ligand_receptor_pairs ,],
         aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig7/scRNA/cell_chat_netVisual_bubble_from.pdf",
       a,height=7,width=4,limitsize = FALSE)
