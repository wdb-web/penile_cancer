
library(Seurat)
library(AUCell)
library(SCENIC)
library(dplyr)
library(grid)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(RcisTarget)
library(dplyr)
library(igraph)
library(ggraph)


dotHeatmap <- function (enrichmentDf,
                        var.x="Topic", var.y="ID", 
                        var.col="FC", col.low="dodgerblue", col.mid="floralwhite", col.high="brown1", 
                        var.size="p.adjust", min.size=1, max.size=8,
                        ...){
  require(data.table)
  require(ggplot2)
  
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  p <- ggplot(data=enrichmentDf, mapping=aes_string(x=var.x, y=var.y)) + 
    geom_point(mapping=aes_string(size=var.size, color=var.col)) +
    scale_radius(range=c(min.size, max.size)) +
    scale_colour_gradientn(colors=colorPal(10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
          axis.text.x=element_text(angle=90, hjust=1),
          ...)
  return(p)
}





obs=read.csv("output/data/TF_loop/Treg.obs.csv")
reg=data.table::fread("output/TF/20000_gene/regulons.csv")
net=data.table::fread("output/TF/Treg_cell_sample_fliter.tsv")
auc=data.table::fread("output/TF/20000_gene/auc_mtx.csv")
cellType=c("high-differentiated","moderate-differentiated","low-differentiated")
# 筛选 TF
umap=read.csv("/home/wangdongbin/work/2024/2024_3_29_E003_colorectal/output/fig2/adata_Treg_umap.csv")

umap=umap[(umap$cell%in%auc$Cell),]
# AUC 分析
a=auc %>% data.frame(check.names = F)
rownames(a)=a[,1]
a=a[,-1]
#
# UMAP color 

rss <- data.frame(AUC=a[umap$cell,],umap,check.names=F)
rownames(obs)=obs$cell
obs=obs[obs$Grades!="",]
a=a[rownames(a)%in% rownames(obs),]
obs$Grades=obs$Grades|>stringr::str_replace("high-differentiated","moderate-and-high-differentiated")|>stringr::str_replace("moderate-differentiated","moderate-and-high-differentiated")

rss <- calcRSS(AUC=(a%>%t),
               cellAnnotation=obs$obs %>%
  dplyr::mutate(cell_annotation = dplyr::case_when(
    leiden == '0' ~ 'Memory CTL',
    leiden == '1' ~ 'Proliferating CTL',
    leiden == '2' ~ 'Stressed CTL',
    leiden == '3' ~ 'ANXA1+ CTL',
    TRUE ~ 'Unknown'  # Optional: catch-all for any values not specified
  ) )%>%.$cell_annotation)

rss=na.omit(rss)
# AUC画图
rssPlot <- plotRSS(rss)
plot=rssPlot$df
plot$Topic= factor(plot$Topic,rssPlot$rowOrder)
j=ggplot(plot,aes(y=Topic,x=cellType,size=RSS,color=Z))+geom_point()+
#scale_color_gradient2( low = '#330066',mid =  '#66CC66',high = '#FFCC33')+ 
theme_classic()+ scale_size(range = c(2, 6))+ theme(
    axis.title.x = element_text(),  # x轴标题加粗
    axis.title.y = element_text(),  # y轴标题加粗
    axis.text.x = element_text(face = "bold", size = 10,angle=45,hjust=1,vjust=1),  # x轴刻度标签加粗
    axis.text.y = element_text(face = "bold", size = 10)  # y轴刻度标签加粗
  ) + viridis::scale_color_viridis(option="D")

ggsave("/home/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/output/fig4/rssPlot_CTL_cell.pdf", j,width =3,height = 13)


# 网络图
net=net[net$TF %in% stringr::str_remove(plot$Topic,"\\(.*"),]

net_plot=net %>%group_by(TF) %>% top_n(wt=importance,n=10)
colnames(net_plot)= c("from","to","importance")

g="STAT5A，ZNF200，HSF1，ZNF23，ZNF580，E2F4，BACH1，ILF2，ZBTB25，RELB，NFKB1，REL，FOSL2，IRF2 ，FOXO1，REL，FOSL2，IRF2 ，FOXO1，ILF2
STAT5A，ZNF200，HSF1，ZNF23，ZNF580，E2F4，BACH1，ILF2，ZBTB25，RELB，NFKB1，REL，FOSL2，IRF2 ，FOXO1，REL，FOSL2，IRF2 ，FOXO1，ILF2"%>%stringr:: str_split("，")%>%unlist()
net_plot=net_plot[net_plot$from%in%g,]


actors <- data.frame(name=c(net_plot$from %>% unique,net_plot$to %>% unique),
lable=c(net_plot$from %>% unique,net_plot$to %>% unique),
                     type=c(rep("TF",net_plot$from %>% unique%>%length),
                     rep("Target",net_plot$to %>% unique %>%length)))

actors$type[actors$name%in%(names(table(actors$name))[table(actors$name)>1])]="TF"
actors=actors[!duplicated(actors),]
rownames(actors)=actors$name
g <- graph_from_data_frame(net_plot, directed=FALSE, vertices=actors)

color_dict=c("TF"="#E3AD68","Target"="#ACD45E")

a=ggraph(g, layout = 'graphopt') + 
  geom_edge_link(alpha = 1, aes(width  = (importance)),
                 label_dodge = unit(0.01, 'mm'),
                 arrow = arrow(type = "closed", length = unit(4, 'mm')),
                 end_cap = circle(0.01, 'mm'),color="darkgray"#, arrow = arrow(length = unit(4, 'mm'))
                 ) +ggraph::scale_edge_size(range = c(0.5, 1))+  #scale_size(range = c(0.0005, 0.001))+ 
  geom_node_point(aes(color=(type)),cex=8,alpha = 1)+  
  scale_color_manual(values = color_dict)+ # custom function for edge color
theme_void() +
  geom_node_text(aes( label = name), family = "serif",cex=1) 
ggsave("output/fig2/转录因子net.pdf",a,height=8,width=10)



rss <- calcRSS(AUC=(a%>%t),
               cellAnnotation=paste(obs[rownames(a),]$Grades,obs[rownames(a),]$cell_type_all_Treg) )

rss=na.omit(rss)
# AUC画图
rssPlot <- plotRSS(rss)
plot=rssPlot$df
plot$Topic= factor(plot$Topic,rssPlot$rowOrder)
j=ggplot(plot,aes(y=Topic,x=cellType,size=RSS,color=Z))+geom_point()+
#scale_color_gradient2( low = '#330066',mid =  '#66CC66',high = '#FFCC33')+ 
theme_classic()+ scale_size(range = c(2, 6))+ theme(
    axis.title.x = element_text(),  # x轴标题加粗
    axis.title.y = element_text(),  # y轴标题加粗
    axis.text.x = element_text(face = "bold", size = 10,angle=45,hjust=0.5,vjust=0.5),  # x轴刻度标签加粗
    axis.text.y = element_text(face = "bold", size = 10)  # y轴刻度标签加粗
  ) + viridis::scale_color_viridis(option="D")

#ggsave("output/fig2/rssPlot_Treg_cell.pdf", j,width =3,height = 13)
j=ggplot(plot,aes(y=Topic,x=cellType,size=RSS,color=Z))+geom_point()+
#scale_color_gradient2( low = '#330066',mid =  '#66CC66',high = '#FFCC33')+ 
theme_classic()+ scale_size(range = c(2, 6))+ theme(
    axis.title.x = element_text(),  # x轴标题加粗
    axis.title.y = element_text(),  # y轴标题加粗
    axis.text.x = element_text(face = "bold", size = 10,angle=45,hjust=0.5,vjust=0.5),  # x轴刻度标签加粗
    axis.text.y = element_text(face = "bold", size = 10)  # y轴刻度标签加粗
  ) + viridis::scale_color_viridis(option="D")

ggsave("output/fig5/rssPlot_Treg_cell.pdf", j,width =6,height = 20)



# 网络图
net=net[net$TF %in% stringr::str_remove(plot$Topic,"\\(.*"),]

net_plot=net %>%group_by(TF) %>% top_n(wt=importance,n=10)
colnames(net_plot)= c("from","to","importance")

g="STAT5A，ZNF200，HSF1，ZNF23，ZNF580，E2F4，BACH1，ILF2，ZBTB25，RELB，NFKB1，REL，FOSL2，IRF2 ，FOXO1，REL，FOSL2，IRF2 ，FOXO1，ILF2
STAT5A，ZNF200，HSF1，ZNF23，ZNF580，E2F4，BACH1，ILF2，ZBTB25，RELB，NFKB1，REL，FOSL2，IRF2 ，FOXO1，REL，FOSL2，IRF2 ，FOXO1，ILF2"%>%stringr:: str_split("，")%>%unlist()
net_plot=net_plot[net_plot$from%in%g,]


actors <- data.frame(name=c(net_plot$from %>% unique,net_plot$to %>% unique),
lable=c(net_plot$from %>% unique,net_plot$to %>% unique),
                     type=c(rep("TF",net_plot$from %>% unique%>%length),
                     rep("Target",net_plot$to %>% unique %>%length)))

actors$type[actors$name%in%(names(table(actors$name))[table(actors$name)>1])]="TF"
actors=actors[!duplicated(actors),]
rownames(actors)=actors$name
g <- graph_from_data_frame(net_plot, directed=FALSE, vertices=actors)

color_dict=c("TF"="#E3AD68","Target"="#ACD45E")

a=ggraph(g, layout ='kk' ) +# 'graphopt'
  geom_edge_link(alpha = 1, aes(width  = importance),
                 label_dodge = unit(0.01, 'mm'),
                 arrow = arrow(type = "closed", length = unit(4, 'mm')),
                 end_cap = circle(0.01, 'mm'),color="darkgray"#, arrow = arrow(length = unit(4, 'mm'))
                 ) +  scale_edge_width_continuous( range=c(1,2.5)) +
  geom_node_point(aes(color=(type)),cex=8,alpha = 1)+  
  scale_color_manual(values = color_dict)+ # custom function for edge color
theme_void() +
  geom_node_text(aes( label = name), family = "serif",cex=1) 
ggsave("output/fig2/转录因子net.pdf",a,height=8,width=10)
# AUC 可视化
reg=reg[-c(1:3),]


list_data=list()
gene=c(auc%>% colnames%>%stringr::str_remove("...$"))
#reg=reg[as.numeric(reg$q)<0.05,]
for(i in 1:nrow(reg)){
reg$V9%>% .[i] %>% stringr::str_split("\\)")%>% unlist()->k
k2=k%>% stringr::str_split_fixed(",|\\[",3)
k2[,3]=k2[,3] %>% as.numeric()
k2=na.omit(k2) %>% as.data.frame
k2[,2]=k2[,2]%>% stringr::str_remove_all("\\(|\\'| ")
k2[,3]=k2[,3] %>% as.numeric()
k2[,1]=reg$V1%>% .[i]
list_data[[i]]=k2
}
list_data=do.call(rbind,list_data)




#TF 网络图和AUC一起

#TF 网络图合并
reg=reg[reg$V1%in%c(stringr::str_remove(colnames(auc),"\\(.*")),]
reg=reg[reg$V7 %>% stringr::str_extract("(?<=q-value =) .*\\)")%>%stringr::str_remove("\\)") %>% 
as.numeric()%>% is.na() %>%{!.},]
reg$q=reg$V7 %>% stringr::str_extract("(?<=q-value =) .*\\)")%>%stringr::str_remove("\\)") %>% 
as.numeric()

list_data=list()
gene=c(auc%>% colnames%>%stringr::str_remove("...$"))
#reg=reg[as.numeric(reg$q)<0.05,]
for(i in 1:nrow(reg)){
reg$V9%>% .[i] %>% stringr::str_split("\\)")%>% unlist()->k
k2=k%>% stringr::str_split_fixed(",|\\[",3)
k2[,3]=k2[,3] %>% as.numeric()
k2=na.omit(k2) %>% as.data.frame
k2[,2]=k2[,2]%>% stringr::str_remove_all("\\(|\\'| ")
k2[,3]=k2[,3] %>% as.numeric()
k2[,1]=reg$V1%>% .[i]
list_data[[i]]=k2
}
list_data=do.call(rbind,list_data)
#TF 网络图和AUC一起
list_data=list_data[list_data$V1 %in% stringr::str_remove(plot$Topic,"\\(.*"),]



net_plot=list_data %>%group_by(V1) %>% top_n(wt=V3,n=20)
colnames(net_plot)= c("from","to","importance")
actors <- data.frame(name=c(net_plot$from %>% unique,net_plot$to %>% unique),
lable=c(net_plot$from %>% unique,net_plot$to %>% unique),
                     type=c(rep("TF",net_plot$from %>% unique%>%length),
                     rep("Target",net_plot$to %>% unique %>%length)))

actors$type[actors$name%in%(names(table(actors$name))[table(actors$name)>1])]="TF"
actors=actors[!duplicated(actors),]
rownames(actors)=actors$name
g <- graph_from_data_frame(net_plot, directed=FALSE, vertices=actors)

color_dict=c("TF"="#E3AD68","Target"="#ACD45E")

a=ggraph(g, layout = 'cactustree') + 
  geom_edge_link(alpha = 1, aes(width  = (importance)/10),
                 label_dodge = unit(0.1, 'mm'),
                 arrow = arrow(type = "closed", length = unit(4, 'mm')),
                 end_cap = circle(0.1, 'mm'),color="darkgray"#, arrow = arrow(length = unit(4, 'mm'))
                 ) + 
  geom_node_point(aes(color=(type)),cex=10,alpha = 1)+  
  scale_color_manual(values = color_dict)+ # custom function for edge color
theme_void()+  scale_size(range = c(0.01, 0.03)) +
  geom_node_text(aes( label = name), family = "serif") 

ggsave("output/fig2/转录因子net.pdf",a,height=25,width=30)
# 导出图数据到GraphML文件
write.graph(g, file = "output/fig2/graph_data.graphml", format = "graphml")


tfs <- c("ETV7(+)","EGR1(+) ","CEBPD(+)","BHLHE40(+)")
par(mfrow=c(2,2))
rss <- data.frame(AUC=a[umap$cell,],umap,check.names=F)
umap_sp=umap[,2:3]%>% as.matrix
rownames(umap_sp)=umap[,1]
k=AUCell::AUCell_plotTSNE(umap_sp, exprMat=a[umap$cell,tfs], plots = "AUC")
# AUC 分析
a=auc %>% data.frame(check.names = F)
rownames(a)=a[,1]
a=a[,-1]
#
# UMAP color 

rss <- data.frame(AUC=a[umap$cell,],umap,check.names=F)
rownames(obs)=obs$cell
obs=obs[obs$Grades!="",]
a=a[rownames(a)%in% rownames(obs),]
# Create the expression matrix specific to the cells and transcription factors
exprMat <- a[umap$cell, tfs]

# Plot the t-SNE with AUC values
