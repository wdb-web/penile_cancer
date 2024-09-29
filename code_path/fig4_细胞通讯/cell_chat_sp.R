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

#cellchat=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cellchat.RDS")  
#names(cellchat_run)=meta.data$sample_id%>%unique
cellchat_run=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/2024_8_27_cellchat_k.min_5.RDS")  
#table(meta.data$sample_id,meta.data$TLS)


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
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/net.pdf",width=12, height=12)
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
cell_chat_count=cellchat_run_data_long%>%
    tidyr::pivot_wider(names_from=name,values_from=exp) %>%
    data.frame(check.names=F)
rownames(cell_chat_count)=cell_chat_count[,1]

pdf("output/fig4_细胞通讯/net.pdf",width=12, height=12)
netVisual_circle(cell_chat_count[c("DC"      ,         "ETS_CD74 B cells"  ,     "TLS_CD74 B cells",     "IGHM "       ,     "IGHV1 "     ,
"IGKV4 "   ,        "IGLV4 "    ,      "Myeloid"   ,       "NK & T"     ),c("DC"      ,         "ETS_CD74 B cells"  ,     "TLS_CD74 B cells",     "IGHM "       ,     "IGHV1 "     ,
"IGKV4 "   ,        "IGLV4 "    ,      "Myeloid"   ,       "NK & T"     )] %>% as.matrix,
                 weight.scale = T, label.edge= F)
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
                        filter(p<0.05)%>% dplyr::mutate(to=stringr::str_remove(name,"\\..*"),
                        path=stringr::str_remove(name,".*\\.")) %>%
                        group_by(from ,to,path ) %>%
                        summarise(count=n(),prob=mean(prob))

bind_net_count_to_wider=bind_net_count[,-5] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = count, values_fill = list(count = 0))


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
#bind_net_count$rate[bind_net_count$from=="ETS_CD74 B cells"]=bind_net_count$count[bind_net_count$from=="ETS_CD74 B cells"]/18
#bind_net_count$rate[bind_net_count$from=="TLS_CD74 B cells"]=bind_net_count$count[bind_net_count$from=="TLS_CD74 B cells"]/13
bind_net_count$rate=#[bind_net_count$from=="ETS_CD74 B cells"]
bind_net_count$count#[bind_net_count$from=="ETS_CD74 B cells"]/18

bind_net_count$rate=#[bind_net_count$from=="TLS_CD74 B cells"]
bind_net_count$count#[bind_net_count$from=="TLS_CD74 B cells"]/13


bind_net_count$from=stringr::str_replace(bind_net_count$from,"ETS","nTLS")
bind_net_count$to=stringr::str_replace(bind_net_count$to,"ETS","nTLS")

bind_net_count=bind_net_count %>% filter(from %in%c("TLS_CD74 B cells","nTLS_CD74 B cells") ,
                          to%in%c("DC" , "IGHM "  ,"IGHV1 " ,"IGKV4 " ,"IGLV4 " ,"Myeloid","NK & T"
                                  ) )


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
ligand_receptor_pairs <- c(
  "C3_CR2",
  "CLEC2C_KLRB1",
  "COL4A3_CD44",
  "COL4A3_ITGA1_ITGB1",
  "COL4A3_ITGA2_ITGB1",
  "FCER2A_CR2",
  "GZMA_F2R",
  "HLA-G_CD8A",
  "ICAM2_ITGAL_ITGB2",
  "ICAM2_ITGAM_ITGB2",
  "IL26_IL20RA_IL10RB",
  "NRG1_ERBB3",
  "RLN3_RXFP4",
  "TNFSF10_TNFRSF10B",
  "ULBP2_KLRK1",
  "ULBP2_NKG2D_HCST",
  "COL9A2_SDC1",
  "EBI3_IL27RA_IL6ST",
  "IL22_IL22RA1_IL10RB",
  "COL4A3_SDC1",
  "COL9A3_ITGA1_ITGB1",
  "COL9A3_ITGA2_ITGB1",
  "IFNB1_IFNAR1_IFNAR2",
  "LGALS9_P4HB",
  "CD274_PDCD1",
  "CXCL16_CXCR6",
  "IL33_IL1RL1_IL1RAP",
  "COL4A3_ITGA3_ITGB1",
  "COL9A3_CD44",
  "COL9A3_SDC1",
  "CTSG_F2R",
  "DLL4_NOTCH3",
  "JAM1_ITGAL_ITGB2",
  "NCAM1_FGFR1",
  "PF4_CXCR3",
  "RETN_CAP1",
  "RETN_TLR4",
  "TNFSF8_TNFRSF8",
  "CD34_SELP",
  "IAPP_CALCRL",
  "HBEGF_EGFR",
  "COL9A3_ITGA3_ITGB1",
  "COL9A2_ITGAV_ITGB8",
  "CTSG_F2RL3",
  "NTN1_UNC5C"
)
bind_net_count=bind_net_count[bind_net_count$path %in%ligand_receptor_pairs,]
bind_net_count$to=stringr::str_remove(bind_net_count$to," Plasma cells")
combos <- expand.grid(from = bind_net_count$from %>% unique, to = bind_net_count$to %>% unique,  stringsAsFactors = FALSE)

scale_max=function(data_scale,max=0.083,min=0.053){
data_scale[data_scale>max]=max
data_scale[data_scale<min]=min
return(data_scale)
}
data_scale=bind_net_count$prob#%>%scale_max
a=ggplot(bind_net_count[bind_net_count$from%in%c("TLS_CD74 B cells","nTLS_CD74 B cells"),],
         aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=scale_max(prob),size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/cell_chat_netVisual_bubble_from.pdf",
       a,height=15,width=6.5,limitsize = FALSE)


new_dotplot_li="ULBP2 -> NKG2D_HCST
ULBP2 -> KLRK1
JAM1 -> ITGAL_ITGB2
IL33 -> IL1RL1_IL1RAP
IL26 -> IL20RA_IL10RB
IL22 -> IL22RA1_IL10RB
ICAM2 -> ITGAL_ITGB2
FGF22 -> FGFR2
EFNA2 -> EPHA4
COL9A2 -> SDC1
COL4A3 -> ITGA1_ITGB1" %>% str_replace_all(" -> ","_") %>% stringr::str_split("\n") %>% unlist
#dd->bind_net_count
bind_net_count=bind_net_count[bind_net_count$path %in%new_dotplot_li,]
bind_net_count$to=stringr::str_remove(bind_net_count$to," Plasma cells")
combos <- expand.grid(from = bind_net_count$from %>% unique, to = bind_net_count$to %>% unique,  stringsAsFactors = FALSE)
a=ggplot(bind_net_count[bind_net_count$from%in%c("TLS_CD74 B cells","nTLS_CD74 B cells"),],
         aes(x= paste(from,"->",to) %>% factor( levels=paste(combos$from,"->",combos$to ) )
    ,y=paste(path %>% stringr::str_replace("_","->")),
color=prob,size=rate) )+geom_point()+theme_bw()+
scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+#ggplot2::scale_color_gradient(low="#fdeff0",high="#b3242a")
ggplot2::theme_bw()+#ggplot2::geom_text( colour = "black", fontface = "bold")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_size(range = c(1, 6))
ggsave("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_netVisual_bubble_from.pdf",
       a,height=4.5,width=6.5,limitsize = FALSE)

# 保存 所有通讯的结果的

library(ggalluvial)
colors20=c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

p=ggplot(bind_net_count,aes(axis1 = from, axis2 = path,
                            axis3 =to,y= count))+
  # scale_x_discrete(limits = c("Diet", "Microbiome", "Meta")) +
  geom_alluvium(aes(fill = bind_net_count$path),cex=0.4)+
  geom_stratum(width = 1/12, fill = "black", color = "grey")+
  #geom_stratum() +
  geom_stratum(width = 0.1) +
  theme_minimal()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
    scale_fill_manual(values = colors20) +
  theme(legend.position="none",axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),panel.grid=element_blank())#+

ggsave("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_alluvial.pdf",p,width = 15,height = 12)


cellchat=cellchat[-which(names(cellchat)%in%c("S12","S13","S16","S17","S4","S5","S9"))]

cellchat=cellchat_run
extractEnrichedLR_list=list()
extractEnrichedLR_list_signaling=list()
for(sample in  names(cellchat)){
    for(signaling in cellchat_run[[sample]]@netP[[1]]){
        extractEnrichedLR_list[[paste(sample,signaling )]]=
          extractEnrichedLR(object=cellchat_run[[sample]], signaling =signaling, geneLR.return = FALSE)
    }
    extractEnrichedLR_list_signaling[[sample]]=cellchat_run[[sample]]@netP[[1]]
}

extractEnrichedLR_list_signaling %>% unlist %>% unique
draw_empty_plot <- function(title) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
}

pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/fig4e_cell_chat_netVisual_aggregate.pdf",width=18,height=15)
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
library(circlize)
color.use = scPalette(length(bind_net_count[,1:2]%>% unlist %>% unique))
names(color.use) <- bind_net_count[,1:2]%>% unlist %>% unique
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/fig4b_cell_chat_netVisual_aggregate.pdf",width=18,height=15)

par(mfrow=c(2,2))
for(title.name in bind_net_count$path[bind_net_count$path%in%ligand_receptor_pairs]%>% unique){

CAF_bind_split=bind_net_count[bind_net_count$path==title.name,]
circos.clear()


# Prepare data
CAF_bind_split <- CAF_bind_split %>%
  select(from, to, n=rate)
net=CAF_bind_split
colnames(net)[1:2]=c("source","target")
#net$plot=T
cell.levels=table(bind_net_count[,1:2]%>% unlist)%>% names

cells.removed <- setdiff(cell.levels, as.character(union(net$source,
            net$target)))

            net.fake <- data.frame(cells.removed, cells.removed,
                1e-10 * sample(length(cells.removed), length(cells.removed)))
            colnames(net.fake) <- colnames(net)
            net <- rbind(net, net.fake)
            link.visible <- net[, 1:2]
            link.visible$plot <- F
            if (nrow(net) > nrow(net.fake)) {
                link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
            }
            scale = TRUE
num_sectors <- length(unique(c(as.character(net$source), as.character(net$target))))
net$source=stringr::str_remove(net$source," Plasma cells")
net$target=stringr::str_remove(net$target," Plasma cells")

# 创建一个与扇区数量相同长度的gap.degree参数
gap.degree <- rep(1, num_sectors)  # 将每个间隙设置为1度，可以根据需要调整
circos.par(gap.after=1,gap.degree=1)
# 创建弦图
  link.visible[is.na(link.visible)] <- FALSE  # 或者 TRUE，根据您的需要进行设置

chordDiagram(
    net[,c(1,2,3)],  order=bind_net_count[,1:2]%>% unlist %>% unique,
    scale = TRUE,  link.visible = net$n>0.000001#link.visible
  ,grid.col =color.use,
link.arr.type = "big.arrow",
    direction.type = c("diffHeight", "arrows"),
    directional = 1,
             annotationTrack = "grid"
)
circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
                niceFacing = TRUE,
adj = c(0.5, 0.5), cex = 0.8)  # 增大cex以提高可见性
}, bg.border = NA)
text(-0, 1.02, title.name, cex = 1)
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
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/output/fig4_细胞通讯/fig4cell_chat_sp/",signaling,".png"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
}


bind_net_count=bind_net %>% #[stringr::str_detect(bind_net$from ,"ETS_CD74 B cells")|stringr::str_detect(bind_net$from ,"TLS_CD74 B cells"),]
                        filter(p<0.05)%>% dplyr::mutate(to=stringr::str_remove(name,"\\..*"),
                        path=stringr::str_remove(name,".*\\.")) %>%
                        group_by(from ,to,path ) %>%
                        summarise(count=n(),prob=mean(prob))

bind_net_count_to_wider=bind_net_count[,-5] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = count, values_fill = list(count = 0))

No_rm_path=bind_net_count_to_wider$path[bind_net_count_to_wider$`TLS_CD74 B cells_NK & T` >=
bind_net_count_to_wider$`ETS_CD74 B cells_NK & T`  ] %>%
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
bind_net_count$from=stringr::str_replace(bind_net_count$from,"ETS","nTLS")
bind_net_count$to=stringr::str_replace(bind_net_count$to,"ETS","nTLS")

bind_net_count=bind_net_count %>% filter(from %in%c("TLS_CD74 B cells","nTLS_CD74 B cells") )


ligand_receptor_pairs <- c(
  "C3_CR2",
  "CLEC2C_KLRB1",
  "COL4A3_CD44",
  "COL4A3_ITGA1_ITGB1",
  "COL4A3_ITGA2_ITGB1",
  "FCER2A_CR2",
  "GZMA_F2R",
  "HLA-G_CD8A",
  "ICAM2_ITGAL_ITGB2",
  "ICAM2_ITGAM_ITGB2",
  "IL26_IL20RA_IL10RB",
  "NRG1_ERBB3",
  "RLN3_RXFP4",
  "TNFSF10_TNFRSF10B",
  "ULBP2_KLRK1",
  "ULBP2_NKG2D_HCST",
  "COL9A2_SDC1",
  "EBI3_IL27RA_IL6ST",
  "IL22_IL22RA1_IL10RB",
  "COL4A3_SDC1",
  "COL9A3_ITGA1_ITGB1",
  "COL9A3_ITGA2_ITGB1",
  "IFNB1_IFNAR1_IFNAR2",
  "LGALS9_P4HB",
  "CD274_PDCD1",
  "CXCL16_CXCR6",
  "IL33_IL1RL1_IL1RAP",
  "COL4A3_ITGA3_ITGB1",
  "COL9A3_CD44",
  "COL9A3_SDC1",
  "CTSG_F2R",
  "DLL4_NOTCH3",
  "JAM1_ITGAL_ITGB2",
  "NCAM1_FGFR1",
  "PF4_CXCR3",
  "RETN_CAP1",
  "RETN_TLR4",
  "TNFSF8_TNFRSF8",
  "CD34_SELP",
  "IAPP_CALCRL",
  "HBEGF_EGFR",
  "COL9A3_ITGA3_ITGB1",
  "COL9A2_ITGAV_ITGB8",
  "CTSG_F2RL3",
  "NTN1_UNC5C"
)

bind_net_count=bind_net_count[bind_net_count$path %in%ligand_receptor_pairs,]


library(circlize)
color.use = scPalette(length(bind_net_count[,1:2]%>% unlist %>% unique))
names(color.use) <- bind_net_count[,1:2]%>% unlist %>% unique
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/fig4b_cell_chat_和弦图_所有细胞.pdf",width=18,height=15)

par(mfrow=c(2,2))
for(title.name in bind_net_count$path[bind_net_count$path%in%ligand_receptor_pairs]%>% unique){

CAF_bind_split=bind_net_count[bind_net_count$path==title.name,]
circos.clear()


# Prepare data
CAF_bind_split <- CAF_bind_split %>%
  select(from, to, n=rate)
net=CAF_bind_split
colnames(net)[1:2]=c("source","target")
#net$plot=T
cell.levels=table(bind_net_count[,1:2]%>% unlist)%>% names

cells.removed <- setdiff(cell.levels, as.character(union(net$source,
            net$target)))

            net.fake <- data.frame(cells.removed, cells.removed,
                1e-10 * sample(length(cells.removed), length(cells.removed)))
            colnames(net.fake) <- colnames(net)
            net <- rbind(net, net.fake)
            link.visible <- net[, 1:2]
            link.visible$plot <- F
            if (nrow(net) > nrow(net.fake)) {
                link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
            }
            scale = TRUE
num_sectors <- length(unique(c(as.character(net$source), as.character(net$target))))
net$source=stringr::str_remove(net$source," Plasma cells")
net$target=stringr::str_remove(net$target," Plasma cells")

# 创建一个与扇区数量相同长度的gap.degree参数
gap.degree <- rep(1, num_sectors)  # 将每个间隙设置为1度，可以根据需要调整
circos.par(gap.after=1,gap.degree=1)
# 创建弦图
  link.visible[is.na(link.visible)] <- FALSE  # 或者 TRUE，根据您的需要进行设置

chordDiagram(
    net[,c(1,2,3)],  order=bind_net_count[,1:2]%>% unlist %>% unique,
    scale = TRUE,  link.visible = net$n>0.000001#link.visible
  ,grid.col =color.use,
link.arr.type = "big.arrow",
    direction.type = c("diffHeight", "arrows"),
    directional = 1,
             annotationTrack = "grid"
)
circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
                niceFacing = TRUE,
adj = c(0.5, 0.5), cex = 0.8)  # 增大cex以提高可见性
}, bg.border = NA)
text(-0, 1.02, title.name, cex = 1)
}
dev.off()


# 原位 

bind_net_count

celltype_color_dict =c('Cancer'= '#1f77b4',
 'B cell'= '#ff7f0e',
 'Fibroblast'= '#2ca02c',
 'Myeloid'='#8c564b',
 'Endothelial'= '#9467bd',
 'TLS_CD74 B cells'=  "#db4d57",
 'ETS_CD74 B cells'= '#e377c2',
 'DC'= '#7f7f7f',
 'NK & T'= '#bcbd22',
 "IGHV1 "=  "#689eac",
  "IGHM " = '#ADB467',
 "IGLV4 "  ='#CECF84' ,
  "IGKV4 "  = '#5DA05A')


split_sp_raw <- list()  # initialize an empty list to store plots
for(sample in names(cellchat_run)){
    plot_title <- paste(sample, "no form Treg to CTL", sep=" - ")
    split_sp_raw[[sample]] <- draw_empty_ggplot(plot_title)
}

source("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/细胞通讯画图定制版本.R")

line=c("FASL_FAS","CD86_CTLA4")
signaling="FASL_FAS"
sub_line=c( "B cell"       ,    "Myeloid"         ,    "TLS_CD74 B cells" ,"ETS_CD74 B cells", "DC"      ,        
"NK & T"         ,  "IGHV1 "    ,       "IGHM "      ,     "IGLV4 "     ,     
"IGKV4 "   )
for(signals in ligand_receptor_pairs ){
    system(paste( "echo ",signals ,">>run.sh"))
split_sp=split_sp_raw
for(sample in c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])
print(paste(sample,signals,sep="-"))
try({split_sp[[paste(signals,sample,sep="-")]]=netVisual_aggregate_check(
                    object=cellchat_run[[sample]],chech.name=signals,
                    alpha.edge=1,
                    #signaling.name =paste(sample,signaling,sep="-"),
                    layout = "spatial",
                    sources.use = c("TLS_CD74 B cells"),targets.use=sub_line,
 edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.09, vertex.label.cex = 5,subcell=sub_line
                    )+
 ggplot2::labs(title=paste(sample,signals,sep="-"))+scale_color_manual(values=celltype_color_dict)
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/cell_chta_imm/",signals,".pdf"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
system(paste( "echo ",signals ,"is ok >>run.sh"))
}

for(signals in ligand_receptor_pairs ){
    system(paste( "echo ",signals ,">>run.sh"))
split_sp=split_sp_raw
for(sample in c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])
print(paste(sample,signals,sep="-"))
try({split_sp[[paste(signals,sample,sep="-")]]=netVisual_aggregate_check(
                    object=cellchat_run[[sample]],chech.name=signals,
                    alpha.edge=1,
                    #signaling.name =paste(sample,signaling,sep="-"),
                    layout = "spatial",
                    sources.use = c("TLS_CD74 B cells"),targets.use=sub_line,
 edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.09, vertex.label.cex = 0,subcell=sub_line
                    )+
 ggplot2::labs(title=paste(sample,signals,sep="-"))+scale_color_manual(values=celltype_color_dict)
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/cell_chta_imm/",signals,".png"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE,dpi=300)
system(paste( "echo ",signals ,"is ok >>run.sh"))
}



for(signals in new_dotplot_li ){
    system(paste( "echo ",signals ,">>run.sh"))
split_sp=split_sp_raw
for(sample in c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])
print(paste(sample,signals,sep="-"))
try({split_sp[[paste(signals,sample,sep="-")]]=netVisual_aggregate_check(
                    object=cellchat_run[[sample]],chech.name=signals,
                    alpha.edge=1,
                    #signaling.name =paste(sample,signaling,sep="-"),
                    layout = "spatial",
                    sources.use = c("TLS_CD74 B cells"),targets.use=sub_line,
 edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.09, vertex.label.cex = 5,subcell=sub_line
                    )+
 ggplot2::labs(title=paste(sample,signals,sep="-"))+scale_color_manual(values=celltype_color_dict)
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/cell_chta_imm/",signals,".pdf"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
system(paste( "echo ",signals ,"is ok >>run.sh"))
}




library(ggplot2)
library(ggnewscale)
for(signals in new_dotplot_li ){
    system(paste( "echo ",signals ,">>run.sh"))
split_sp=split_sp_raw
a=spatialFeaturePlot(cellchat_run[[sample]], pairLR.use =  signals,
 point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
plot_=ggplot()+geom_point(data=a[[1]]$data[a[[1]]$data$feature.data>0.001,],aes(x=x,y=y,color=feature.data))+
  # 第一种颜色渐变
  scale_color_gradient(low =  "#FBDFE2", high ="#B83945", name = stringr::str_remove(signals,"_.*")) +
  # 添加新的颜色尺度
  ggnewscale::new_scale("color") +
  # 第二种颜色渐变
    geom_point(data=a[[2]]$data[a[[2]]$data$feature.data>0.001,],aes(x=x,y=y,color=feature.data))+
  scale_color_gradient(low = "#CFE7C4", high = "#4F845C", name = stringr::str_remove(signals,".*?_")) +
  theme_minimal()+ labs(title=signals)+ggplot2::coord_fixed()



ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位/",signals,".pdf"),
plot_, width =16*2, height = 20*2,limitsize = FALSE)
system(paste( "echo ",signals ,"is ok >>run.sh"))
}

spatialFeaturePlot(cellchat, pairLR.use = "IGF1_IGF1R", point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

for(signals in ligand_receptor_pairs ){
    system(paste( "echo ",signals ,">>run.sh"))
split_sp=split_sp_raw
for(sample in c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")){
    #cellchat_run[[sample]]=netAnalysis_computeCentrality(cellchat_run[[sample]])
print(paste(sample,signals,sep="-"))

try({
a=spatialFeaturePlot(cellchat_run[[sample]], pairLR.use =  signals,
 point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

split_sp[[paste(signals,sample,sep="-")]]=
ggplot()+geom_point(data=a[[1]]$data[a[[1]]$data$feature.data>0.001,],aes(x=x,y=y,color=feature.data))+
  # 第一种颜色渐变
  scale_color_gradient(low =  "#FBDFE2", high ="#B83945", name = stringr::str_remove(signals,"_.*")) +
  # 添加新的颜色尺度
  ggnewscale::new_scale("color") +
  # 第二种颜色渐变
    geom_point(data=a[[2]]$data[a[[2]]$data$feature.data>0.001,],aes(x=x,y=y,color=feature.data))+
  scale_color_gradient(low = "#CFE7C4", high = "#4F845C", name = stringr::str_remove(signals,".*?_")) +
  theme_minimal()+
 ggplot2::labs(title=paste(sample,signals,sep="-"))#+scale_color_manual(values=celltype_color_dict)
                    })
}
cowplot::ggsave2(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位/",signals,".pdf"),
cowplot::plot_grid(plotlist=split_sp), width = 20*2, height =16*2,limitsize = FALSE)
system(paste( "echo ",signals ,"is ok >>run.sh"))
}


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

B_meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
adata=read_h5ad("output/data/adata.h5ad")
Tcell=read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
meta.data<- Tcell$obs
x=t(Tcell$X)

a.meta.data=read_h5ad("output/data/dbscan.h5ad")$obs
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  
meta.data$TLS=ifelse(a.meta.data$dbscan_labels[match(rownames(meta.data),rownames(a.meta.data))]!=-1 ,"TLS","nTLS")  

X=t(adata$X)
adata_meta.data=adata$obs
adata_meta.data$cell_type=as.character(adata_meta.data$cell_type)

meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs

adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub=="CD74 B cells"]=paste0(
  ifelse(meta.data$TLS=="Yes","TLS","nTLS"),
"_",meta.data$cell_type_sub)[meta.data$cell_type_sub=="CD74 B cells"]
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub!="CD74 B cells"]=paste0(
meta.data$cell_type_sub)[meta.data$cell_type_sub!="CD74 B cells"] %>% stringr::str_remove("Plasma cells")
gene=X[unique(ligand_receptor_pairs ) %>%stringr::str_split("_") %>% unlist %>%unique%>%{.[.%in%rownames(X)]} ,]
gene_exp=data.frame(x=adata_meta.data$x_slide_mm[adata_meta.data$sample_id%in% c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")],
y=adata_meta.data$y_slide_mm[adata_meta.data$sample_id%in% c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")],
cell_type=adata_meta.data$cell_type[adata_meta.data$sample_id%in% c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19")],
t(gene)[adata_meta.data$sample_id%in% c("S2"  ,"S3" , "S6"  ,"S7",  "S8" ,"S10", "S11" ,"S14" ,"S15" ,"S18", "S19"),])



for(signals in ligand_receptor_pairs){
    
    
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
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位/",signals,".pdf"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
})
}











# 基因表达 
COL4A3_CD44
COL4A3_ITGA1_ITGB1
EBI3_IL27RA_IL6ST
signaling=c("EBI3_IL27RA_IL6ST","COL4A3_ITGA1_ITGB1","COL4A3_CD44")
signaling[1]->signals
for(signals in ligand_receptor_pairs%>% make.names ){
try({
  system(paste( "echo ",signals ,">>run.sh"))
split_st=stringr::str_split(signals,"_")[[1]]
ggplot_obj=ggplot()+
geom_point(data=gene_exp,aes_string(x="x",y="y"),color="#DCDCDC")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" ,],aes_string(x="x",y="y"),color="#A9A9A9")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" & gene_exp[,split_st[1]]>0,],
aes_string(x="x",y="y",color=split_st[1]),size=5)+
  # 第一种颜色渐变
scale_color_gradient(low =  "#FBDFE2", high ="#B83945", name = split_st[1]) +theme_classic()+ggplot2::coord_fixed()
  # 添加新的颜色尺度
  # 第二种颜色渐变
if(length(split_st)==2){
ggplot_obj=ggplot_obj+ ggnewscale::new_scale("color") + 
geom_point(data=gene_exp[gene_exp$cell_type %in% "NK & T"& gene_exp[,split_st[2]]>0,], 
aes_string(x="x",y="y",color=split_st[2]),size=5) +
  scale_color_gradient(low = "#CFE7C4", high = "#4F845C", name = split_st[2]) +
 ggplot2::labs(title=signals)   
  }
#+scale_color_manual(values=celltype_color_dict)
  if(length(split_st)>2){
     gene_exp$ bing_gene_exp=0
 gene_exp$bing_gene_exp[gene_exp$cell_type %in% "NK & T"]= 
  gene_exp[gene_exp$cell_type %in% "NK & T",split_st[2]]+gene_exp[gene_exp$cell_type %in% "NK & T",split_st[3]]

ggplot_obj=ggplot_obj+      ggnewscale::new_scale("color") +
geom_point(data=gene_exp[gene_exp$cell_type %in% "NK & T"&  gene_exp$bing_gene_exp>0,],aes_string(x="x",y="y",color="bing_gene_exp"),size=5)+
  scale_color_gradient(low = "#CFE7C4", high = "#4F845C", name =paste0(split_st[2],"_", split_st[3])) 
  }
#ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".pdf"),ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".png"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
})
}



# COL4A3_CD44
signals="COL4A3_CD44"
    system(paste( "echo ",signals ,">>run.sh"))
split_st=stringr::str_split(signals,"_")[[1]]
ggplot_obj=ggplot()+
geom_point(data=gene_exp,aes_string(x="x",y="y"),color="#DCDCDC")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" ,],aes_string(x="x",y="y"),color="#A9A9A9")

  # 添加新的颜色尺度
  # 第二种颜色渐变

ggplot_obj=ggplot_obj+ geom_point(data=gene_exp[gene_exp$cell_type %in% "NK & T"& gene_exp[,split_st[2]]>0,]%>%.[order(.[,split_st[2]]),], aes_string(x="x",y="y",color=split_st[2]),size=5) +
  scale_color_gradient(low = "#FFFFFF", high = '#C311FF', name = split_st[2]) +
 ggplot2::labs(title=paste(sample,signals,sep="-"))   +ggnewscale::new_scale("color") + 
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" & gene_exp[,split_st[1]]>0,]%>%.[order(.[,split_st[1]]),],
aes_string(x="x",y="y",color=split_st[1]),size=5)+
  # 第一种颜色渐变
scale_color_gradient(low =  "#FFFFFF", high ='#00FF1A', name = split_st[1]) +theme_classic()
  

ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".pdf"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".png"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)





# EBI3_IL27RA_IL6ST
signals="EBI3_IL27RA_IL6ST"
    system(paste( "echo ",signals ,">>run.sh"))
split_st=stringr::str_split(signals,"_")[[1]]
ggplot_obj=ggplot()+
geom_point(data=gene_exp,aes_string(x="x",y="y"),color="#DCDCDC")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" ,],aes_string(x="x",y="y"),color="#A9A9A9")

  # 添加新的颜色尺度
  # 第二种颜色渐变
#+scale_color_manual(values=celltype_color_dict)
 
     gene_exp$ bing_gene_exp=0
 gene_exp$bing_gene_exp[gene_exp$cell_type %in% "NK & T"]=  gene_exp[gene_exp$cell_type %in% "NK & T",split_st[2]]+gene_exp[gene_exp$cell_type %in% "NK & T",split_st[3]]

ggplot_obj=ggplot_obj+      
geom_point(data=gene_exp[gene_exp$cell_type %in% "NK & T"&  gene_exp$bing_gene_exp>0,]%>%.[order(.[,"bing_gene_exp"]),],aes_string(x="x",y="y",color="bing_gene_exp"),size=5)+
  scale_color_gradient(low = "#FFFFFF", high = '#003fff', name =paste0(split_st[2],"_", split_st[3])) +ggnewscale::new_scale("color") +
 ggplot2::labs(title=paste(signals)) +  geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" & gene_exp[,split_st[1]]>0,]%>%.[order(.[,split_st[1]]),],
aes_string(x="x",y="y",color=split_st[1]),size=5)+
  # 第一种颜色渐变
scale_color_gradient(low =  "#FFFFFF", high ='#ff7f0e', name = split_st[1]) +theme_classic()
  
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".pdf"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".png"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)




# COL4A3_ITGA1_ITGB1
signals="COL4A3_ITGA1_ITGB1"
system(paste( "echo ",signals ,">>run.sh"))
split_st=stringr::str_split(signals,"_")[[1]]
ggplot_obj=ggplot()+
geom_point(data=gene_exp,aes_string(x="x",y="y"),color="#DCDCDC")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" ,],aes_string(x="x",y="y"),color="#A9A9A9")+
geom_point(data=gene_exp[gene_exp$cell_type %in% "TLS_CD74 B cells" & gene_exp[,split_st[1]]>0,]%>%.[order(.[,split_st[1]]),],
aes_string(x="x",y="y",color=split_st[1]),size=5)+
  # 第一种颜色渐变
scale_color_gradient(low =  "#FFFFFF", high ="#00EAFF", name = split_st[1]) +theme_classic()
  # 添加新的颜色尺度

#+scale_color_manual(values=celltype_color_dict)
  if(length(split_st)>2){
     gene_exp$ bing_gene_exp=0
 gene_exp$bing_gene_exp[gene_exp$cell_type %in% "NK & T"]=  gene_exp[gene_exp$cell_type %in% "NK & T",split_st[2]]+gene_exp[gene_exp$cell_type %in% "NK & T",split_st[3]]

ggplot_obj=ggplot_obj+      ggnewscale::new_scale("color") +
geom_point(data=gene_exp[gene_exp$cell_type %in% "NK & T"&  gene_exp$bing_gene_exp>0,]%>%.[order(.[,"bing_gene_exp"]),],aes_string(x="x",y="y",color="bing_gene_exp"),size=5)+
  scale_color_gradient(low = "#FFFFFF", high = "#d62728", name =paste0(split_st[2],"_", split_st[3])) +
 ggplot2::labs(title=paste(signals))   
}
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".pdf"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)
ggsave(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/原位过滤/",signals,".png"),
ggplot_obj, width =16*2, height = 20*2,limitsize = FALSE)



library(cowplot)

# 使用 draw_colorbar 函数
draw_colorbar(ggplot2::scale_fill_gradient(low = "blue", high = "red"))
library(ComplexHeatmap)
# 创建自定义颜色函数
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# 定制 color bar
ht <- Heatmap(matrix(rnorm(100), 10), 
              col = col_fun, 
              heatmap_legend_param = list(
                title = "Expression",
                at = c(-2, 0, 2),
                labels = c("Low", "Medium", "High"),
                direction = "horizontal"
              ))
ComplexHeatmap::
# 只绘制 color bar
draw(ht, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")
library(ComplexHeatmap)
c("#FFFFFF", "#00FF1A")
c("#FFFFFF", "#C311FF")
c("#FFFFFF", "#ff7f0e")
c("#FFFFFF", "#003fff")
c("#FFFFFF", "#00EAFF")
c("#FFFFFF", "#d62728")
# 创建颜色函数
library(circlize)
col_fun1 <- colorRamp2(c( 0, 2), c( c("#FFFFFF", "#00FF1A")))
col_fun2 <- colorRamp2(c( 0, 2), c( c("#FFFFFF", "#C311FF")))
col_fun3<- colorRamp2(c( 0, 2), c( c("#FFFFFF", "#ff7f0e")))
col_fun4 <- colorRamp2(c( 0, 2), c( c("#FFFFFF", "#003fff")))
col_fun5 <- colorRamp2(c( 0, 2), c( c("#FFFFFF", "#00EAFF")))
col_fun6 <- colorRamp2(c( 0, 2), c( c("#FFFFFF", "#d62728")))


# 创建图例对象
lgd1 <- Legend(
  col_fun = col_fun1,               # 颜色映射函数
  title = "COL4A3_CD44 COL4A3 Expression",            # 图例标题
  at = c(0, 2),                # 标尺的刻度位置
  labels = c("Low", "High"),  # 刻度标签
  direction = "horizontal",        # 图例方向
  legend_width = unit(4, "cm")     # 图例宽度
)
# 创2建图例对象
lgd2 <- Legend(
  col_fun = col_fun2,               # 颜色映射函数
  title = "COL4A3_CD44 CD44 Expression",            # 图例标题
  at = c(0, 2),                # 标尺的刻度位置
  labels = c("Low", "High"),  # 刻度标签
  direction = "horizontal",        # 图例方向
  legend_width = unit(4, "cm")     # 图例宽度
)
# 创建图例对象
lgd3 <- Legend(
  col_fun = col_fun3,               # 颜色映射函数
  title = "EBI3_IL27RA_IL6ST EBI3 Expression",            # 图例标题
  at = c(0, 2),                # 标尺的刻度位置
  labels = c("Low", "High"),  # 刻度标签
  direction = "horizontal",        # 图例方向
  legend_width = unit(4, "cm")     # 图例宽度
)
# 创建图例对象
lgd4 <- Legend(
  col_fun = col_fun4,               # 颜色映射函数
  title = "EBI3_IL27RA_IL6ST IL27RA_IL6ST Expression",            # 图例标题
  at = c(0, 2),                # 标尺的刻度位置
  labels = c("Low", "High"),  # 刻度标签
  direction = "horizontal",        # 图例方向
  legend_width = unit(4, "cm")     # 图例宽度
)
# 创建图例对象
lgd5 <- Legend(
  col_fun = col_fun5,               # 颜色映射函数
  title = "COL4A3_ITGA1_ITGB1 COL4A3 Expression",            # 图例标题
  at = c(0, 2),                # 标尺的刻度位置
  labels = c("Low", "High"),  # 刻度标签
  direction = "horizontal",        # 图例方向
  legend_width = unit(4, "cm")     # 图例宽度
)
# 创建图例对象
lgd6 <- Legend(
  col_fun = col_fun6,               # 颜色映射函数
  title = "COL4A3_ITGA1_ITGB1 ITGA1_ITGB1 Expression",            # 图例标题
  at = c(0, 2),                # 标尺的刻度位置
  labels = c("Low", "High"),  # 刻度标签
  direction = "horizontal",        # 图例方向
  legend_width = unit(4, "cm")     # 图例宽度
)


# 创建单独的图例
lgds <- packLegend(lgd1, lgd2, lgd3, lgd4, lgd5, lgd6)

# 绘制图例到PDF文件
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig4_细胞通讯/color_bar.pdf")
draw(lgds)
dev.off()
