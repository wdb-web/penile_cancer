library(dplyr)
library(Seurat)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/home/dell/.miniconda/envs/py_scRNA/bin/python')
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
library(CellChat)
library(scatterpie)
library(CellChat)
library(anndata)
library(dplyr)
library(Seurat)
reticulate::use_python('/home/dell/.miniconda/envs/py_scRNA/bin/python')
future::plan("multisession", workers =10)  # Reduce the number of workers
options(future.globals.maxSize = 10000 * 1024^2) # 将最大值设为600 MiB
future::plan("multisession", workers =10)  # Reduce the number of workers
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")

meta.data=adata=read_h5ad("output/data/adata_B_TLS.h5ad")$obs
adata=read_h5ad("output/data/adata.h5ad")

X=t(adata$X)
adata_meta.data=adata$obs
adata_meta.data$cell_type=as.character(adata_meta.data$cell_type)
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))]=paste0(meta.data$TLS,"_",meta.data$cell_type_sub)
Seurat_data=CreateSeuratObject(X,meta.data =adata_meta.data )




CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB
) # use Secreted Signaling

conversion.factor = 0.18

# 创建进度条
options(future.globals.maxSize = 50 * 1024^3)
future::plan("multisession", workers =20)  # Reduce the number of workers


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
future::plan("multisession", workers =2)  # Reduce the number of workers


cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 422
#>



# Project the data
cellchat <- smoothData(cellchat, adj = PPI.human)
#cellchat <- projectData(cellchat, PPI.human)
#cellchat2 <- computeCommunProb(cellchat, type = "truncatedMean")
cellchat2 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
   distance.use = TRUE,raw.use = TRUE,
                              interaction.range = 250,
                            scale.distance =1000,k.min = 10, nboot = 100,
                              contact.dependent = TRUE, contact.range = 100)
gc()

cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
#saveRDS(cellchat2,paste0("/home/wangdongbin/项目/2024/hongkong_chinese_smi/hongkong_chinese_smi/output/cell_chat/ZFP36/",i,".RDS"))
return(cellchat2)
    }, future.seed = TRUE)


names(cellchat_run)=meta.data$sample_id%>%unique
saveRDS(cellchat_run,paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/cell_chat.RDS"))
cellchat_run=readRDS(paste0("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/cell_chat.RDS"))


pdf("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cellchat.pdf",width=20, height=10)
# Loop through each sample ID in 'cellchat_run'
for (i in names(cellchat_run)) {
  # Extract the CellChat object for the current sample
  cellchat_obj <- cellchat_run[[i]]

  # Perform desired analysis or visualization on the CellChat object
  # For example, you might want to visualize communication networks
 groupSize <- as.numeric(table(cellchat_obj@idents))
 par(mfrow = c(1,2), xpd=TRUE)
 netVisual_circle(cellchat_obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste("Number of interactions", i) )
 netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste("Interaction weights/strength", i) )
  # You can add more analyses or visualizations here
  # For example, netVisual_bubble or netAnalysis_contribution
}
# Close the PDF device
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
cell_chat_count=a%>% purrr::reduce(rbind) %>% group_by(form , name )%>% dplyr::summarise(exp=sum(value)) %>% tidyr::pivot_wider(names_from=name,values_from=exp) %>% data.frame(check.names=F)
rownames(cell_chat_count)=cell_chat_count[,1]
netVisual_circle(cell_chat_count[,-1] %>% as.matrix, weight.scale = T, label.edge= F)

# Sum the count matrices across all samples using purrr::Reduce
summed_counts <- purrr::reduce(a,`+`)
# Initialize summed_counts as a zero matrix with the same dimensions as the first element of list 'a'
summed_counts <- matrix(0, nrow = nrow(a[[9]]), ncol = ncol(a[[9]]))
rownames(summed_counts) <- rownames(a[[9]])
colnames(summed_counts) <- colnames(a[[9]])

# Loop through each dataframe in list 'a' and sum the counts
for (i in 1:length(a)) {
  # Ensure the rownames and colnames match
  if (all(rownames(a[[i]]) == rownames(summed_counts)) && all(colnames(a[[i]]) == colnames(summed_counts))) {
    summed_counts <- summed_counts + a[[i]]
  } else {
    stop("Rownames or colnames do not match in dataframe ", i)
  }
}

# Print the res
cellchat_merged <- mergeCellChat(cellchat_run, add.names = TRUE)
# 更新和处理合并后的CellChat对象

# 可视化合并后的数据，例如：
netVisual_bubble(cellchat_merged, remove.isolate = FALSE)
# Inspect the interaction_name_2 column
head(cellchat_merged@net$interaction_name_2)

# Sort the data manually
df.net <- cellchat_merged@net
df.net <- df.net[order(df.net$interaction_name_2), ]

# Ensure that the sorted data is correctly ordered
df.net_list=list()
for(i in names(cellchat_run)){
    cat(i)
        df.net <- subsetCommunication(cellchat_run[[i%>%as.character]], slot.name = "net", 
            sources.use = NULL, targets.use = NULL, 
            signaling = NULL, pairLR.use = NULL, thresh = 1)
        df.net$source.target <- paste(df.net$source, df.net$target, 
            sep = " -> ")
        df.net$sample=i
df.net_list[[i%>%as.character]]=df.net
}
# Now, try visualizing again
net=do.call(rbind,df.net_list)
net=net[net$target%in%c("No_CD74 B cells"    ,    "No_IGHM Plasma cells" ,  "No_IGHV1 Plasma cells" ,
 "No_IGKV4 Plasma cells" , "No_IGLV4 Plasma cells" , "Yes_CD74 B cells"   ,   
"Yes_IGHM Plasma cells"  ,"Yes_IGHV1 Plasma cells" ,"Yes_IGKV4 Plasma cells",
 "Yes_IGLV4 Plasma cells") &net$source%in%
 c("Cancer"         ,        "DC"    ,                 "Endothelial"         ,  
"Fibroblast"        ,     "Myeloid"       ,         "NK & T"  ),]

 c("Cancer"         ,        "DC"    ,                 "Endothelial"         ,  
"Fibroblast"        ,     "Myeloid"       ,         "NK & T"  ) %>% lapply(\(.x){
  paste0(c("No_CD74 B cells"    ,  "Yes_CD74 B cells"   ,   
 "No_IGHM Plasma cells" , "Yes_IGHM Plasma cells"  , 
 "No_IGHV1 Plasma cells" ,"Yes_IGHV1 Plasma cells" ,
 "No_IGKV4 Plasma cells" , "Yes_IGKV4 Plasma cells",
 "No_IGLV4 Plasma cells" ,  "Yes_IGLV4 Plasma cells"
 )," - ",.x )
})%>% unlist

net=net %>% group_by(source,target,source.target,interaction_name_2) %>% dplyr::summarise(prob=prob%>%mean,count=sum(pval<0.05))

        net$prob[net$prob == 0] <- NA
        net$prob.original <- net$prob
        net$prob <- -1/log(net$prob)
        angle = c(0, 45, 90)
        hjust = c(0, 1, 1)
        vjust = c(0, 1, 0.5)
        vjust.x = vjust[angle == 90]
        hjust.x = hjust[angle == 90]
net$source.target=factor(net$source.target, c("Cancer"         ,        "DC"    ,                 "Endothelial"         ,  
"Fibroblast"        ,     "Myeloid"       ,         "NK & T"  ) %>% lapply(\(.x){
  paste0(.x," -> ",c("No_CD74 B cells"    ,  "Yes_CD74 B cells"   ,   
 "No_IGHM Plasma cells" , "Yes_IGHM Plasma cells"  , 
 "No_IGHV1 Plasma cells" ,"Yes_IGHV1 Plasma cells" ,
 "No_IGKV4 Plasma cells" , "Yes_IGKV4 Plasma cells",
 "No_IGLV4 Plasma cells" ,  "Yes_IGLV4 Plasma cells"
 ) )
})%>% unlist)
    g <- ggplot(net, aes(x = source.target, y = interaction_name_2, 
        color = prob, size = count)) + geom_point(pch = 16) + 
        theme_linedraw() + theme(panel.grid.major = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, hjust = hjust.x, 
            vjust = vjust.x), axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
             
      g <- g + scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Spectral"))(99), 
            na.value = "white", limits = c(quantile(net$prob, 
                0, na.rm = T), quantile(net$prob, 1, na.rm = T)), 
            breaks = c(quantile(net$prob, 0, na.rm = T), quantile(net$prob, 
                1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob."))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_netVisual_bubble_to.pdf",g,width=15,height=40,limitsize = FALSE)


# Now, try visualizing again
net=do.call(rbind,df.net_list)
net=net[net$source%in%c("No_CD74 B cells"    ,    "No_IGHM Plasma cells" ,  "No_IGHV1 Plasma cells" ,
 "No_IGKV4 Plasma cells" , "No_IGLV4 Plasma cells" , "Yes_CD74 B cells"   ,   
"Yes_IGHM Plasma cells"  ,"Yes_IGHV1 Plasma cells" ,"Yes_IGKV4 Plasma cells",
 "Yes_IGLV4 Plasma cells") &net$target%in%
 c("Cancer"         ,        "DC"    ,                 "Endothelial"         ,  
"Fibroblast"        ,     "Myeloid"       ,         "NK & T"  ),]

 c("Cancer"         ,        "DC"    ,                 "Endothelial"         ,  
"Fibroblast"        ,     "Myeloid"       ,         "NK & T"  ) %>% lapply(\(.x){
  paste0(c("No_CD74 B cells"    ,  "Yes_CD74 B cells"   ,   
 "No_IGHM Plasma cells" , "Yes_IGHM Plasma cells"  , 
 "No_IGHV1 Plasma cells" ,"Yes_IGHV1 Plasma cells" ,
 "No_IGKV4 Plasma cells" , "Yes_IGKV4 Plasma cells",
 "No_IGLV4 Plasma cells" ,  "Yes_IGLV4 Plasma cells"
 )," - ",.x )
})%>% unlist

net=net %>% group_by(source,target,source.target,interaction_name_2) %>% dplyr::summarise(prob=prob%>%mean,count=sum(pval<0.05))

        net$prob[net$prob == 0] <- NA
        net$prob.original <- net$prob
        net$prob <- -1/log(net$prob)
        angle = c(0, 45, 90)
        hjust = c(0, 1, 1)
        vjust = c(0, 1, 0.5)
        vjust.x = vjust[angle == 90]
        hjust.x = hjust[angle == 90]
        net$source.target=factor(net$source.target, c("Cancer"         ,        "DC"    ,                 "Endothelial"         ,  
"Fibroblast"        ,     "Myeloid"       ,         "NK & T"  ) %>% lapply(\(.x){
  paste0(c("No_CD74 B cells"    ,  "Yes_CD74 B cells"   ,   
 "No_IGHM Plasma cells" , "Yes_IGHM Plasma cells"  , 
 "No_IGHV1 Plasma cells" ,"Yes_IGHV1 Plasma cells" ,
 "No_IGKV4 Plasma cells" , "Yes_IGKV4 Plasma cells",
 "No_IGLV4 Plasma cells" ,  "Yes_IGLV4 Plasma cells"
 )," -> ",.x )
})%>% unlist)
    g <- ggplot(net, aes(x = source.target, y = interaction_name_2, 
        color = prob, size = count)) + geom_point(pch = 16) + 
        theme_linedraw() + theme(panel.grid.major = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, hjust = hjust.x, 
            vjust = vjust.x), axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
             
      g <- g + scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Spectral"))(99), 
            na.value = "white", limits = c(quantile(net$prob, 
                0, na.rm = T), quantile(net$prob, 1, na.rm = T)), 
            breaks = c(quantile(net$prob, 0, na.rm = T), quantile(net$prob, 
                1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob."))
ggsave("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3/cell_chat_netVisual_bubble_from.pdf",g,width=15,height=40,limitsize = FALSE)

