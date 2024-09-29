library(dplyr)
library(Seurat)
library(ggplot2)
library(CellChat)
library(tidyr)
library(circlize)
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
future::plan("multisession", workers =10)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
future::plan("multisession", workers =10)  # Reduce the number of workers
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")

cellchat_run=readRDS("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cellchat.RDS")
#table(meta.data$sample_id,meta.data$TLS)


a=list()

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
bind_net_count=bind_net[stringr::str_detect(bind_net$from ,"No_CD74 B cells")|stringr::str_detect(bind_net$from ,"Yes_CD74 B cells"),] %>%
                        filter(p<0.05)%>% dplyr::mutate(to=stringr::str_remove(name,"\\..*"),path=stringr::str_remove(name,".*\\.")) %>%
                        group_by(from ,to,path ) %>%
                        summarise(count=n(),prob=mean(prob))

bind_net_count_to_wider=bind_net_count[,-5] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = count, values_fill = list(count = 0))

No_rm_path=bind_net_count_to_wider$path[bind_net_count_to_wider$`Yes_CD74 B cells_NK & T`/15 >
bind_net_count_to_wider$`No_CD74 B cells_NK & T`/18  ] %>%
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
bind_net_count$rate[bind_net_count$from=="No_CD74 B cells"]=bind_net_count$count[bind_net_count$from=="No_CD74 B cells"]/18
bind_net_count$rate[bind_net_count$from=="Yes_CD74 B cells"]=bind_net_count$count[bind_net_count$from=="Yes_CD74 B cells"]/13

bind_net_count$from=stringr::str_replace(bind_net_count$from,"No","NOTLS")%>%
stringr::str_replace("Yes","TLS")%>%
stringr::str_remove_all(" Plasma cells")
bind_net_count$to=stringr::str_replace(bind_net_count$to,"No","NOTLS")%>%
stringr::str_replace("Yes","TLS")%>%
stringr::str_remove_all(" Plasma cells")
bind_net_count_to_wider=bind_net_count[,c(1,2,3,7)] %>% # 为了筛选 pos和 NK 更高的 样本
  pivot_wider(names_from = c(from, to), values_from = rate, values_fill = list(rate = 0))

bind_net_count=bind_net_count[bind_net_count$path %in% bind_net_count_to_wider$path[bind_net_count_to_wider$`TLS_CD74 B cells_NK & T`/bind_net_count_to_wider$`NOTLS_CD74 B cells_NK & T` >1.4],]

color.use = scPalette(length(bind_net_count[,1:2]%>% unlist %>% unique))
names(color.use) <- bind_net_count[,1:2]%>% unlist %>% unique
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
bind_net_count$path=stringr::str_replace(bind_net_count$path,"_"," -> ")
pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/fig3_circos.track.pdf",width=12, height=12)

par(mfrow=c(2,2))
for(title.name in bind_net_count$path[bind_net_count$path%in%migration_pairs]%>% unique){

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
