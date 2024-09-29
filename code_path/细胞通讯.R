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

meta.data=adata=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
adata=read_h5ad("output/data/adata.h5ad")

X=t(adata$X)
adata_meta.data=adata$obs
adata_meta.data$cell_type=as.character(adata_meta.data$cell_type)

adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub=="CD74 B cells"]=paste0(
  ifelse(meta.data$TLS=="Yes","TLS","ETS"),
"_",meta.data$cell_type_sub)[meta.data$cell_type_sub=="CD74 B cells"]
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub!="CD74 B cells"]=paste0(
meta.data$cell_type_sub)[meta.data$cell_type_sub!="CD74 B cells"] %>% stringr::str_remove("Plasma cells")
table(adata_meta.data$cell_type,adata_meta.data$cell_type)
meta.data
table(meta.data$cell_type_sub,meta.data$TLS)

Seurat_data=CreateSeuratObject(X,meta.data =adata_meta.data )




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
p <- progressr::progressor(steps = length(meta.data$sample_id%>%unique))
# 使用 future_apply 并传递进度条
 cellchat_run=future.apply::future_lapply(meta.data$sample_id%>%unique,\(.x){

.x=as.character(.x)
meta=adata_meta.data[adata_meta.data$sample_id==.x,]
data.input=X[,rownames(meta)]
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
cellchat2 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05,
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
saveRDS(cellchat_run,"/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/2024_8_27_cellchat_k.min_5.RDS")  

#cellchat=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cellchat.RDS")  
names(cellchat_run)=meta.data$sample_id%>%unique
cellchat_run=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cellchat_k.min_3.RDS")  
#table(meta.data$sample_id,meta.data$TLS)


a=list()
for (i in names(cellchat_run)) {
  # Extract the CellChat object for tshe current sample
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
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/net.pdf",width=12, height=12)
netVisual_circle(cell_chat_count[,-1] %>% as.matrix, weight.scale = T, label.edge= F)
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

bind_net_count_to_wider=bind_net_count[,-5] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = count, values_fill = list(count = 0))

No_rm_path=bind_net_count_to_wider$path[bind_net_count_to_wider$`TLS_CD74 B cells_NK & T`/13 >
bind_net_count_to_wider$`ETS_CD74 B cells_NK & T`/18  ] %>%
 unique
#get_max_count=bind_net_count%>%group_by(path)%>%summarise(max_count=max(count)>5)
#get_max_count=get_max_count$path[get_max_count$max_count]
bind_net_count=bind_net_count[bind_net_count$path %in%No_rm_path, ]
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
bind_net_count$rate[bind_net_count$from=="ETS_CD74 B cells"]=bind_net_count$count[bind_net_count$from=="ETS_CD74 B cells"]/18
bind_net_count$rate[bind_net_count$from=="TLS_CD74 B cells"]=bind_net_count$count[bind_net_count$from=="TLS_CD74 B cells"]/13
bind_net_count$rate=#[bind_net_count$from=="ETS_CD74 B cells"]
bind_net_count$count#[bind_net_count$from=="ETS_CD74 B cells"]/18
bind_net_count$rate=#[bind_net_count$from=="TLS_CD74 B cells"]
bind_net_count$count#[bind_net_count$from=="TLS_CD74 B cells"]/13

bind_net_count$from=stringr::str_replace(bind_net_count$from,"No","NOTLS")%>%stringr::str_replace("Yes","TLS")
bind_net_count$to=stringr::str_replace(bind_net_count$to,"No","NOTLS")%>%stringr::str_replace("Yes","TLS")
bind_net_count_to_wider=bind_net_count[,c(1,2,3,7)] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = rate, values_fill = list(rate = 0))

bind_net_count=bind_net_count[bind_net_count$path %in% bind_net_count_to_wider$path[
  bind_net_count_to_wider$`TLS_CD74 B cells_NK & T`/bind_net_count_to_wider$`NOTLS_CD74 B cells_NK & T` >1.4],]

# 

path="1. CXCL13->CXCR3
2. CXCL13->CXCR5
3. CXCL8->CXCR2
4. CCL19->CCR7
5. CD40LG->ITGA5_ITGB1
6. CD40LG->ITGAM_ITGB2
7. CD99->CD99
8. CRTAM->CADM1
9. LGALS9->CD44
10. LGALS9->HAVCR2
11. SELL->PODXL
"%>% stringr::str_remove_all("\\d+\\. ")%>% stringr::str_split("\n")%>% unlist
c("CX3CL1->CX3CR1"      ,
"CXCL1->CXCR1"       ,  "CXCL13->CXCR3"   ,     "CXCL13->CXCR5"       ,
"CXCL2->CXCR1"     ,    "CXCL2->CXCR2"  ,
"COL4A1->GP6"      ,    "COL4A1->ITGA2_ITGB1" ,
"COL4A1->ITGA3_ITGB1"  ,"COL4A1->ITGA9_ITGB1",  "COL4A1->SDC1"      ,  
 "COL4A3->CD44"        , "COL4A3->ITGA1_ITGB1" , "COL4A3->ITGA2_ITGB1" ,
 "COL4A3->ITGA3_ITGB1" , "COL4A3->SDC1"      ,   "COL4A4->CD44"        ,
"COL4A4->ITGA1_ITGB1" , "COL4A4->ITGA2_ITGB1" , "COL6A1->ITGA9_ITGB1" ,
"COL6A1->SDC1"       ,  "COL6A2->GP6"       ,   "COL6A3->ITGA1_ITGB1" ,
 "COL6A6->GP6"     ,     "COL6A6->ITGA1_ITGB1" , "COL9A2->CD44"    ,    
"COL9A2->ITGA1_ITGB1",  "COL9A2->SDC1"  )
migration_pairs <- c(
  "CCL19 -> CCR7",  # T细胞迁移
  "CXCL13 -> CXCR5",  # T细胞与B细胞相互作用
  "CXCL1 -> CXCR1",  # 可能与T细胞迁移相关
  "CCL3 -> CCR5",  # T细胞和NK细胞的招募
  "CCL8 -> CCR1",  # 与T细胞和NK细胞相关
  "ULBP1 -> KLRK1",  # 与NK细胞相关
  "ULBP2 -> KLRK1",  # 与NK细胞相关
  "ULBP2 -> NKG2D_HCST",  # 与NK细胞相关
  "CX3CL1 -> CX3CR1"  ,# 与T细胞和NK细胞迁移相关
  "CD80 -> CD28",  # T细胞激活
  "CD80 -> CTLA4",  # T细胞调控
  "CD40LG -> CD40",  # T细胞与抗原呈递细胞相互作用
  "IL2 -> IL2RB_IL2RG",  # T细胞和NK细胞激活
  "IL4 -> IL4R_IL13RA1",  # T细胞激活
  "IL10 -> IL10RA_IL10RB",  # 调控T细胞激活
  "LGALS9 -> CD44",  # T细胞激活
  "LGALS9 -> HAVCR2",  # T细胞激活
  "ULBP1 -> KLRK1",  # NK细胞激活
  "ULBP2 -> KLRK1",  # NK细胞激活
  "ULBP2 -> NKG2D_HCST",  # NK细胞激活
  "COL4A1 -> ITGA2_ITGB1",
  "COL4A1 -> ITGA3_ITGB1",
  "COL4A1 -> ITGA9_ITGB1",
  "COL4A1 -> SDC1",
  "COL4A3 -> ITGA1_ITGB1",
  "COL4A3 -> ITGA2_ITGB1",
  "COL4A3 -> ITGA3_ITGB1",
  "COL4A3 -> SDC1",
  "COL4A4 -> ITGA1_ITGB1",
  "COL4A4 -> ITGA2_ITGB1",
  "COL6A1 -> ITGA9_ITGB1",
  "COL6A1 -> SDC1",
  "COL6A3 -> ITGA1_ITGB1",
  "COL6A6 -> ITGA1_ITGB1",
  "COL9A2 -> ITGA1_ITGB1",
  "COL9A2 -> SDC1"
)
bind_net_count=bind_net_count[bind_net_count$path %>% stringr::str_replace("_"," -> ")%>%{ .%in%migration_pairs},]
bind_net_count$to=stringr::str_remove(bind_net_count$to," Plasma cells")
combos <- expand.grid(from = bind_net_count$from %>% unique, to = bind_net_count$to %>% unique,  stringsAsFactors = FALSE)
a=ggplot(bind_net_count[bind_net_count$from%in%c("TLS_CD74 B cells","ETS_CD74 B cells"),],aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_netVisual_bubble_from.pdf",
       a,height=50,width=14,limitsize = FALSE)
 
combos <- expand.grid(from = bind_net_count$from %>% unique, to = bind_net_count$to %>% unique,  stringsAsFactors = FALSE)
a=ggplot(bind_net_count[bind_net_count$to%in%c("TLS_CD74 B cells","ETS_CD74 B cells"),],
aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$to,"->",combos$from ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_netVisual_bubble_to.pdf",
       a,height=50,width=14,limitsize = FALSE)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)


extractEnrichedLR_list=list()
extractEnrichedLR_list_signaling=list()
for(sample in  names(cellchat_run)){
    for(signaling in cellchat_run[[sample]]@netP[[1]]){
        extractEnrichedLR_list[[paste(sample,signaling )]]=extractEnrichedLR(cellchat_run[[sample]], signaling =signaling, geneLR.return = FALSE)
    }
    extractEnrichedLR_list_signaling[[sample]]=cellchat_run[[sample]]@netP[[1]]
}

extractEnrichedLR_list_signaling %>% unlist %>% unique
draw_empty_plot <- function(title) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
}

pdf("output/fig3/fig3e_cell_chat_netVisual_aggregate.pdf",width=18,height=15)
par(mfrow = c(5,5), xpd=TRUE)
for(pathways.show in extractEnrichedLR_list_signaling%>% unlist %>% unique){
for (i in 1:length(cellchat_run)) {
    tryCatch ({  netVisual_aggregate(cellchat_run[[i]], signaling = pathways.show, layout = "circle", 
  edge.weight.max = 30, edge.width.max = 20, signaling.name = paste(pathways.show, names(cellchat_run)[i]),weight.scale =F)
  }, error = function(e) {
      draw_empty_plot(paste("Error in", pathways.show, names(cellchat_run)[i]))
  })
}
}
dev.off()
source("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/code_path/cell_chat.aplha.R")

split_sp=list()
for(sample in names(cellchat_run)){
    cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])
}
for(signaling in extractEnrichedLR_list_signaling%>% unlist %>% unique){


    split_sp=list()
    for(sample in names(cellchat_run)%>%sort){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])

print(paste(sample,signaling,sep="-"))
try({split_sp[[paste(signaling,sample,sep="-")]]=netVisual_aggregate(
                    object=cellchat_run[[sample]],
                    signaling = signaling,alpha.edge=1,
                    #signaling.name =paste(sample,signaling,sep="-"),
                    layout = "spatial",
                    sources.use = c("Yes_CD74 B cells", "No_CD74 B cells"),
 edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.05, vertex.label.cex = 3.5
                    )+ggplot2::labs(title=paste(sample,signaling,sep="-"))
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_sp/",signaling,".png"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
}



source("/home/wangdongbin/work/packages/标准库/分析/细胞通讯/cell_chat_sp.R")



object=net_agg;chech.name="FGF18_FGFR1";
                                    signaling.name = NULL; color.use = NULL;
subcell=NULL;
    thresh = 0.05; vertex.receiver = NULL; sources.use = NULL;
    targets.use = NULL; idents.use = NULL; top = 1; remove.isolate = FALSE;
    vertex.weight = 1; vertex.weight.max = NULL; vertex.size.max = NULL;
    weight.scale = TRUE; edge.weight.max = NULL; edge.width.max = 8;
    layout = c( "spatial"); pt.title = 12;
    title.space = 6; vertex.label.cex = 0.8; sample.use = NULL;
    alpha.image = 0.15; point.size = 1.5; group = NULL; cell.order = NULL;
    small.gap = 1; big.gap = 10; scale = FALSE; reduce = -1;
    show.legend = FALSE; legend.pos.x = 20; legend.pos.y = 20;max.scale=5;


migration_pairs <- "LGALS9->CD44
IL4->IL4R_IL13RA1
IL2->IL2RB_IL2RG
CXCL13->CXCR5
CXCL1->CXCR1
CD80->CTLA4
CD80->CD28
CCL19->CCR7" %>% stringr::str_split("\n")%>% unlist

split_sp=list()
for(sample in names(cellchat_run)){
    cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])
}
par(mfrow=c(1,1))

for(signaling in stringr::str_replace(migration_pairs,"->","_")){


    split_sp=list()
    for(sample in names(cellchat_run)%>%sort){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])

print(paste(sample,signaling,sep="-"))
try({split_sp[[paste(signaling,sample,sep="-")]]=netVisual_aggregate_check(
                    object=cellchat_run[[sample]],chech.name=signaling,
                    alpha.edge=1,
                    #signaling.name =paste(sample,signaling,sep="-"),
                    layout = "spatial",
                    sources.use = c("Yes_CD74 B cells", "No_CD74 B cells"),
 edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.05, vertex.label.cex = 3.5
                    )+ggplot2::labs(title=paste(sample,signaling,sep="-"))
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_sp2/",signaling,".pdf"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
}
for(signaling in stringr::str_replace(migration_pairs,"->","_")){


    split_sp=list()
    for(sample in names(cellchat_run)%>%sort){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])

print(paste(sample,signaling,sep="-"))
try({split_sp[[paste(signaling,sample,sep="-")]]=netVisual_aggregate_check(
                    object=cellchat_run[[sample]],chech.name=signaling,
                    alpha.edge=1,
                    #signaling.name =paste(sample,signaling,sep="-"),
                    layout = "spatial",
                    sources.use = c("Yes_CD74 B cells", "No_CD74 B cells"),
 edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.05, vertex.label.cex = 3.5
                    )+ggplot2::labs(title=paste(sample,signaling,sep="-"))
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_sp2/",signaling,".png"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
}
