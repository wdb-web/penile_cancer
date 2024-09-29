#/home/wangdongbin/.conda/envs/mistyr/bin/R
library(dplyr)
#library(SeuratDisk)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
library(anndata)
library("Matrix")
#library(Seurat)
library(future.apply)
library(ggplot2)
library(mistyR)
library(decoupleR)
library(janitor)
library(purrr)
library(rstatix)
library(ggpubr)
library(tibble)
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

# meta.data$cell_type_all 转变成T细胞亚类和CAF转变成 TLS_CD74.B.cells、ETS_CD74.B.cells
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub=="CD74 B cells"]=paste0(
  ifelse(meta.data$TLS=="Yes","TLS","ETS"),
"_",meta.data$cell_type_sub)[meta.data$cell_type_sub=="CD74 B cells"]
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub!="CD74 B cells"]=paste0(
meta.data$cell_type_sub)[meta.data$cell_type_sub!="CD74 B cells"] %>% stringr::str_remove("Plasma cells")

meta.data=adata_meta.data

# Expression data
# Here the dgRMatrix is converted to a dense matrix for vst compatibility reasons
# 保存xy 
geometry <- meta.data[,c('x_slide_mm', 'y_slide_mm')]
colnames(geometry) <- c("row", "col")
# Obtain genesets
# Use multivariate linear model to estimate activity
# MISTy 视图
# 合并数据
synthetic=data.frame(geometry,cell_type=meta.data$cell_type)
# 数据拆分
synthetic2=split( synthetic,meta.data$sample_id)
run_misty_bing <- function(name) {
#sample.expr <- synthetic2[[name]] %>% select(-c(row, col,cell_type))
  # 转变成为 标准字体 make.names会改变对应的 CAF名字，只能转化后才能维持原样
  sample_expr_cell_type <- synthetic2[[name]] %>% summarise(cell_type=cell_type %>%stringr::str_replace("\\+","_pos_") %>%stringr::str_replace("-","_neg_"))
  # 转化为哑变量
  sample.expr <- model.matrix(~ cell_type-1 , data = sample_expr_cell_type)
  # 坐标
  sample.pos <- synthetic2[[name]] %>% select(row, col)
  # 去掉方差为0的数据，这些数据会报错
  non_zero_variance <- sample.expr%>% data.frame %>% select_if(~var(.) != 0) # 去掉方差为0的数据，这些数据会报错
  # 哑变量 名字
  colnames(non_zero_variance)=colnames(non_zero_variance)%>%stringr::str_remove("cell_type") 
  # 运行 misty
  create_initial_view(non_zero_variance) %>%
  add_paraview(sample.pos, l = 30, family = "gaussian") %>%
  run_misty(paste0("results/", name))
  # 读取 misty
  misty.results <- collect_results(paste0("results/", name)) #read misty
  import=misty.results$importances.aggregated
  import$sample=name
  return(import)
}
run_read_misty_bing <- function(name) {
#sample.expr <- synthetic2[[name]] %>% select(-c(row, col,cell_type))
  misty.results <- collect_results(paste0("results/", name)) #read misty
  import=misty.results$importances.aggregated
  import$sample=name
  return(import)
}

result_folders_list <- names(synthetic2) %>% map(run_misty_bing)
result_folders_list <- names(synthetic2) %>% map(run_read_misty_bing)

library(dplyr)
#library(SeuratDisk)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
library(anndata)
library("Matrix")
#library(Seurat)
library(future.apply)
library(ggplot2)
library(mistyR)
library(decoupleR)
library(janitor)
library(purrr)
library(rstatix)
library(ggpubr)
library(OmnipathR)
library(tdiyr)
reticulate::use_python('/home/wangdongbin/.conda/envs/mistyr/bin/python')
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
setwd("/home/wangdongbin/data/2024/2024_7_13_E007_TLS")

meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs
adata=read_h5ad("output/data/adata.h5ad")

X=t(adata$X)
adata_meta.data=adata$obs
adata_meta.data$cell_type=as.character(adata_meta.data$cell_type)

# meta.data$cell_type_all 转变成T细胞亚类和CAF转变成 TLS_CD74.B.cells、ETS_CD74.B.cells
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub=="CD74 B cells"]=paste0(
  ifelse(meta.data$TLS=="Yes","TLS","ETS"),
"_",meta.data$cell_type_sub)[meta.data$cell_type_sub=="CD74 B cells"]
adata_meta.data$cell_type[match(rownames(meta.data),rownames(adata_meta.data))][meta.data$cell_type_sub!="CD74 B cells"]=paste0(
meta.data$cell_type_sub)[meta.data$cell_type_sub!="CD74 B cells"] %>% stringr::str_remove("Plasma cells")

meta.data=adata_meta.data

# Expression data
# Here the dgRMatrix is converted to a dense matrix for vst compatibility reasons
# 保存xy 
geometrys <- meta.data[,c('x_slide_mm', 'y_slide_mm')]
colnames(geometrys) <- c("row", "col")
# Obtain genesets
# Use multivariate linear model to estimate activity
# MISTy 视图
# 合并数据
synthetic=data.frame(geometrys,model.matrix(~ cell_type-1 , data = meta.data),X%>%t %>%as.matrix)
# 数据拆分
synthetic2=split( synthetic,meta.data$sample_id)
run_misty_gene_bing <- function(name) {
#sample.expr <- synthetic2[[name]] %>% select(-c(row, col,cell_type))
  # 转变成为 标准字体 make.names会改变对应的 CAF名字，只能转化后才能维持原样
  # 转化为哑变量
  sample.expr <- synthetic2[[name]] %>% select(!c(row, col))
  # 坐标
  sample.pos <- synthetic2[[name]] %>% select(row, col)
  # 去掉方差为0的数据，这些数据会报错
  non_zero_variance <- sample.expr%>% data.frame %>% select_if(~var(.) != 0) # 去掉方差为0的数据，这些数据会报错
  # 哑变量 名字
  colnames(non_zero_variance)=colnames(non_zero_variance)%>%stringr::str_remove("cell_type") 
  # 运行 misty
  create_initial_view(non_zero_variance) %>%
  add_paraview(sample.pos, l = 30, family = "gaussian") %>%
  run_misty(paste0("results/", name))
  # 读取 misty
  misty.results <- collect_results(paste0("results_gene/", name)) #read misty
  import=misty.results$importances.aggregated
  import$sample=name
  return(import)
}

run_read_misty_gene_bing <- function(name) {
#sample.expr <- synthetic2[[name]] %>% select(-c(row, col,cell_type))
  misty.results <- collect_results(paste0("results_gene/", name)) #read misty
  import=misty.results$importances.aggregated
  import$sample=name
  return(import)
}

result_folders_list <- names(synthetic2) %>% map(run_misty_gene_bing)
result_folders_list <- names(synthetic2) %>% map(run_read_misty_gene_bing)





# Matrix with important genes for each pathway
model <- get_progeny(organism = "human", top = 1000)
lig_rec <- import_intercell_network(interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra')),
                         transmitter_param = list(parent = 'ligand'),
                         receiver_param = list(parent = 'receptor'))
gene_names <- rownames(as.matrix(X)[(rowSums(as.matrix(
X%>%as.matrix > 0) / 
ncol(as.matrix(X)))) >= 0.05,]) 

for(sample in unique(meta.data$sample)){
# Use multivariate linear model to estimate activity
expression=X[,meta.data$sample==sample]%>%as.matrix
# Get unique ligands
ligands <- unique(lig_rec$source_genesymbol)
# Get expression of ligands in slide
slide_markers <- ligands[ligands %in% gene_names] 
ligand_expr <- t(as.matrix(expression[slide_markers,])) %>% clean_names()
geometry=geometrys[meta.data$sample==sample,]
#est_path_act <- run_mlm(X[,meta.data$sample==sample]%>%as.matrix, model[1:100,],.mor = NULL) 

model_matrix=model.matrix(~ cell_type-1 , data = meta.data[meta.data$sample_id==sample,]) %>% clean_names()
pathway_activity <- ligand_expr
pathway_act_view <- create_initial_view(as_tibble(pathway_activity) ) %>%
  add_paraview(geometry, l = 30, family = "constant")
ligand_view <- create_initial_view(as_tibble(model_matrix)  %>% clean_names()) %>%
  add_paraview(geometry, l = 30, family = "constant")
combined_views <- pathway_act_view %>% add_views(create_view("paraview.ligand.30", ligand_view[["paraview.30"]]$data, "para.ligand.30"))
run_misty(combined_views, "result/functional_ligand/" %>% paste0(sample)) 
#clean names
}         

# 基因表达情况
for(sample in unique(meta.data$sample)){
# Use multivariate linear model to estimate activity
expression=X[,meta.data$sample==sample]%>%as.matrix
# Get unique ligands
ligands <- unique(lig_rec$source_genesymbol)
gene_names=rownames(X)
# Get expression of ligands in slide
slide_markers <- ligands[ligands %in% gene_names] 
ligand_expr <- t(as.matrix(expression[slide_markers,])) %>% clean_names()
geometry=geometrys[meta.data$sample==sample,]
est_path_act <- ligand_expr %>% data.frame%>% tibble::rownames_to_column("condition") %>% pivot_longer(!condition)  %>% data.frame

est_path_act %>%  
# Put estimated pathway activities object into the correct format

est_path_act_wide <- est_path_act %>% 
  pivot_wider(id_cols = condition, names_from = source, values_from = score) %>%
  column_to_rownames("condition") 

# Clean names
colnames(est_path_act_wide)  <- est_path_act_wide %>% 
  clean_names(parsing_option = 0) %>% 
  colnames(.)

# Create a Seurat object
#seurat_vs[['progeny']] <- CreateAssayObject(counts = t(est_path_act_wide))

# Format for running MISTy later
pathway_activity <- t(as.matrix( t(est_path_act_wide)))
# Get ligands
# Get unique ligands
ligands <- unique(lig_rec$source_genesymbol)
gene_names <- rownames(expression[(rowSums(expression > 0) / ncol(expression)) >= 0.05,]) 
# Get expression of ligands in slide
slide_markers <- ligands[ligands %in% gene_names] 
ligand_expr <- t(as.matrix(expression[slide_markers,])) %>% clean_names()
model_matrix=model.matrix(~ cell_type-1 , data = meta.data[meta.data$sample_id==sample,]) %>% clean_names()
ligand_expr_cell_type=cbind(model_matrix,ligand_expr)
pathway_act_view <- create_initial_view(as_tibble(pathway_activity) ) %>%
  add_paraview(geometry, l = 30, family = "constant")
ligand_view <- create_initial_view(as_tibble(model_matrix)  %>% clean_names()) %>%
  add_paraview(geometry, l = 30, family = "constant")
combined_views <- pathway_act_view %>% add_views(create_view("paraview.ligand.30", ligand_view[["paraview.30"]]$data, "para.ligand.30"))
run_misty(combined_views, "result/functional_ligand/" %>% paste0(sample)) 
#clean names
}         


misty.results <- collect_results(list.files("result/functional_ligand/",full.names=T))

misty.results %>% plot_interaction_communities("intra",cutoff = 0.5)

summary(misty.results)

misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>% 
  summarize(mean.p = mean(value)) %>%
  arrange(mean.p)

misty.results %>% plot_view_contributions()

misty.results %>%
  plot_interaction_heatmap(view = "para.ligand.30", clean = FALSE, #trim = 0.3,
                           trim.measure = "gain.R2", cutoff=0)
#intra
misty.results %>%
  plot_interaction_heatmap(view = "intra", 
                           clean = TRUE, 
                           trim = 0,
                           trim.measure = "gain.R2",
                           cutoff = 0)
misty_results %>%
  plot_interaction_heatmap(view = "intra", clean = TRUE, #trim = 0.3,
                           trim.measure = "gain.R2", cutoff=0)
ggplot(geometry,aes(x=row,y=col))+
geom_point(aes(color=ligand_view$paraview.10$data$col6a3))+scale_colour_gradient2()


ggplot(geometry,aes(x=row,y=col))+
geom_point(aes(color=ligand_view$paraview.10$data$col6a3))+scale_colour_gradient2()

ggplot(geometry,aes(x=row,y=col))+
geom_point(aes(color=ligand_view$paraview.10$data$tgfb))+  scale_colour_gradient2()
#
result_folders_data=do.call(rbind,result_folders_list) %>% na.omit()

library(ggplot2)
library(dplyr)

# 计算Importance的平均值并绘图
a=result_folders_data %>%
  group_by(view, Predictor=Predictor%>% factor(levels=c("DC"          ,     "Endothelial"  ,    "Fibroblast"      ,
"IGHM."     ,       "IGHV1."    ,       "IGKV4.",           "IGLV4."  ,        
 "Myeloid"       ,   "NK...T"     ,"Cancer"   ,    "ETS_CD74.B.cells","TLS_CD74.B.cells"     )), 
 Target=Target%>% factor(c("DC"          ,     "Endothelial"  ,    "Fibroblast"      ,
"IGHM."     ,       "IGHV1."    ,       "IGKV4.",           "IGLV4."  ,        
 "Myeloid"       ,   "NK...T"     ,"Cancer"   ,    "ETS_CD74.B.cells","TLS_CD74.B.cells"     ))) %>%
  summarise(Importance = mean(Importance)) %>%
  ggplot(aes(x = Predictor, y = Target, color = Importance,size=abs(Importance))) +
  geom_point() +theme_bw() +
  facet_grid(~view) +scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 11,name="Spectral"))(99) %>%rev,
            na.value = "white")+
  #scale_color_gradientn(colors = rainbow(7)) +  # 彩虹渐变色
  theme(axis.text.x = element_text(angle = 45, hjust = 1))# X轴标签旋转45度

ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig3/fig3G_misty.pdf",a, width =12/3*2, height =6/3*2)

result_folders_data=do.call(rbind,result_folders_list) %>% na.omit()
stat.test <-   result_folders_data[result_folders_data$view == "intra" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ] %>%
  group_by(Target) %>%
  t_test( Importance~Predictor ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "Predictor")
intra = ggpaired(
  result_folders_data[result_folders_data$view == "intra" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ],
  x = "Predictor", y = "Importance", id = "sample", palette = "jco",
  color = "Predictor", line.color = "gray", line.size = 0.4, facet.by = "Target"
) + #ggplot2::facet_grid(~Target) +
  labs(title = "intra", x = "Predictor", y = "Importance") +
  scale_color_manual(values = c("#47B8B2", "#B8474D")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +stat_pvalue_manual(stat.test,label = "{p}")
  #stat_compare_means(method = "wilcox.test", paired = TRUE, label.y=2.1,label = "p.format", label.x = 1.5)
stat.test <-   result_folders_data[result_folders_data$view != "intra" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ] %>%
  group_by(Target) %>%
  t_test( Importance~Predictor ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "Predictor")
bxp + stat_pvalue_manual(stat.test)
# For para
para = ggpaired(
  result_folders_data[result_folders_data$view == "para.30" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ],
  x = "Predictor", y = "Importance", id = "sample", palette = "jco",
  color = "Predictor", line.color = "gray", line.size = 0.4, facet.by = "Target"
) + #ggplot2::facet_grid(~Target) +
  labs(title = "para.30", x = "Predictor", y = "Importance") +
  scale_color_manual(values = c("#47B8B2", "#B8474D")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_pvalue_manual(stat.test,label = "{p}")
  #stat_compare_means(method = "wilcox.test", paired = TRUE, label.y=2.1,label = "p.format", label.x = 1.5) 


# Box plot facetted by "dose"

# Saving the plots
cowplot::ggsave2(
  "/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig3/fig3H_misty_Target.pdf",
  cowplot::plot_grid(intra, para, ncol = 1),
  width = 10,height=20
)


intra=ggplot(result_folders_data[result_folders_data$view=="intra" &result_folders_data$Predictor%in%c("TLS_CD74.B.cells", "ETS_CD74.B.cells"),],
aes(x=Predictor,y=Importance,color=Predictor))+geom_boxplot()+
 ggplot2::facet_grid(~Target)+ggplot2::scale_color_manual(values=c("blue","red"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

para=ggplot(result_folders_data[result_folders_data$view=="para.30" &result_folders_data$Predictor%in%c("TLS_CD74.B.cells", "ETS_CD74.B.cells"),],
aes(x=Predictor,y=Importance,color=Predictor))+geom_boxplot()+
 ggplot2::facet_grid(~Target)+ggplot2::scale_color_manual(values=c("blue","red"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cowplot::ggsave2("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig3/fig3H_misty_Predictor.pdf",cowplot::plot_grid(intra,para,ncol=1))




library(hdf5r,lib.loc="/home/dell/.miniconda/envs/R_scRNA/lib/R/library")
library(dplyr)
#library(SeuratDisk)
library(ggpubr)

#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/home/dell/.miniconda/envs/py_scRNA/bin/python')
library(anndata)
library("Matrix")
#library(Seurat)
library(future.apply)
library(ggplot2)
library(decoupleR)
library(janitor)
library(purrr)
library(rstatix)

plan("multisession", workers =20)
options(future.globals.maxSize = 200 * 1024^3)  # 增加到2 GiB

adata=read_h5ad("/home/wangdongbin/data/2024/hongkong_chinese_smi/hongkong_chinese_smi/output/fig6/dis_adata_smi.h5ad")


dis_exp=adata$obs[c('group',"cell_subtype",'ETS_CD74.B.cells-CTL-recruitscore','TLS_CD74.B.cells-CTL-recruitscore')]

dis_exp=dis_exp[dis_exp$cell_subtype%in%c("TLS_CD74.B.cells" ,"ETS_CD74.B.cells"),]
dis_exp[is.nan(dis_exp[,3]),3]=0
dis_exp[is.nan(dis_exp[,4]),4]=0
test=dis_exp %>% tidyr::pivot_longer(!c( group, cell_subtype))%>%group_by( group ,cell_subtype) %>% summarise(value=mean(value))
test=rbind(test,data.frame(group=101, cell_subtype="ETS_CD74.B.cells",  value=0))
test$value=test$value
# For intra
stat.test <-   result_folders_data[result_folders_data$view == "intra" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ] %>%
  group_by(Target) %>%
  t_test( Importance~Predictor ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "Predictor")
intra = ggpaired(
  result_folders_data[result_folders_data$view == "intra" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ],
  x = "Predictor", y = "Importance", id = "sample", palette = "jco",
  color = "Predictor", line.color = "gray", line.size = 0.4, facet.by = "Target"
) + #ggplot2::facet_grid(~Target) +
  labs(title = "intra", x = "Predictor", y = "Importance") +
  scale_color_manual(values = c("#47B8B2", "#B8474D")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +stat_pvalue_manual(stat.test,label = "{p.adj}")
  #stat_compare_means(method = "wilcox.test", paired = TRUE, label.y=2.1,label = "p.format", label.x = 1.5)
stat.test <-   result_folders_data[result_folders_data$view != "intra" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ] %>%
  group_by(Target) %>%
  t_test( Importance~Predictor ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "Predictor")
bxp + stat_pvalue_manual(stat.test)
# For para
para = ggpaired(
  result_folders_data[result_folders_data$view == "para.30" & 
                        result_folders_data$Predictor %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells") & 
                        !result_folders_data$Target %in% c("TLS_CD74.B.cells", "ETS_CD74.B.cells"), ],
  x = "Predictor", y = "Importance", id = "sample", palette = "jco",
  color = "Predictor", line.color = "gray", line.size = 0.4, facet.by = "Target"
) + #ggplot2::facet_grid(~Target) +
  labs(title = "para.30", x = "Predictor", y = "Importance") +
  scale_color_manual(values = c("#47B8B2", "#B8474D")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+stat_pvalue_manual(stat.test,label = "{p.adj}")
  #stat_compare_means(method = "wilcox.test", paired = TRUE, label.y=2.1,label = "p.format", label.x = 1.5) 


# Box plot facetted by "dose"

# Saving the plots
cowplot::ggsave2(
  "/home/wangdongbin/data/2024/hongkong_chinese_smi/hongkong_chinese_smi/output/fig6/fig6G_misty.pdf",
  cowplot::plot_grid(intra, para, ncol = 1),
  width = 10,height=20
)


intra=ggplot(result_folders_data[result_folders_data$view=="intra" &result_folders_data$Predictor%in%c("TLS_CD74.B.cells", "ETS_CD74.B.cells"),],
aes(x=Predictor,y=Importance,color=Predictor))+geom_boxplot()+
 ggplot2::facet_grid(~Target)+ggplot2::scale_color_manual(values=c("blue","red"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

para=ggplot(result_folders_data[result_folders_data$view=="para.30" &result_folders_data$Predictor%in%c("TLS_CD74.B.cells", "ETS_CD74.B.cells"),],
aes(x=Predictor,y=Importance,color=Predictor))+geom_boxplot()+
 ggplot2::facet_grid(~Target)+ggplot2::scale_color_manual(values=c("blue","red"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cowplot::ggsave2("/home/wangdongbin/data/2024/hongkong_chinese_smi/hongkong_chinese_smi/output/fig6/fig6G_misty.pdf",cowplot::plot_grid(intra,para,ncol=1))


## misty.results improvements
misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>% 
  summarize(mean.p = mean(value)) %>%
  arrange(mean.p)


model <- get_progeny(organism = "human", top = 1000)
# Use multivariate linear model to estimate activity
est_path_act <- run_mlm(expression[,meta.data$group==sample], model,.mor = NULL) 
composition <- as_tibble(geometry[meta.data$group==sample,])
est_path_act_wide <- est_path_act %>% 
  pivot_wider(id_cols = condition, names_from = source, values_from = score) %>%
  column_to_rownames("condition") 

# Clean names
colnames(est_path_act_wide)  <- est_path_act_wide %>% 
  clean_names(parsing_option = 0) %>% 
  colnames(.)

# Create a Seurat object
seurat_vs[['progeny']] <- CreateAssayObject(counts = t(est_path_act_wide))

# Format for running MISTy later
pathway_activity <- t(est_path_act_wide)
# Get ligands
lig_rec <- OmnipathR::import_intercell_network(interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra')),
                         transmitter_param = list(parent = 'ligand'),
                         receiver_param = list(parent = 'receptor'))

# Get unique ligands
ligands <- unique(lig_rec$source_genesymbol)

# Get expression of ligands in slide
slide_markers <- ligands[ligands %in% gene_names] 
ligand_expr <- t(as.matrix(expression[slide_markers,])) %>% clean_names()

#clean names

pathway_act_view <- create_initial_view(as_tibble(pathway_activity) ) %>%
  add_paraview(geometry, l = 10, family = "constant")

ligand_view <- create_initial_view(as_tibble(ligand_expr)  %>% clean_names()) %>%
  add_paraview(geometry, l = 10, family = "constant")
combined_views <- pathway_act_view %>% add_views(create_view("paraview.ligand.10", ligand_view[["paraview.10"]]$data, "para.ligand.10"))
run_misty(combined_views, "result/functional_ligand")

misty_results <- collect_results("result/functional_ligand/")
misty_data <- mistyR::create_view  (counts = t(norm.data)[,meta.data$group==sample]%>% as.matrix, spatial = geometry[meta.data$group==sample,])
misty.intra <- create_initial_view( expression[,meta.data$group==sample]%>% as.matrix%>%t)
geometry <- adata$obs[,c('x_slide_mm', 'y_slide_mm')]
colnames(geometry) <- c("row", "col")
misty.views <- misty.intra %>% add_paraview(geometry[meta.data$group==sample,], l = 10)

neighbors <- nearest_neighbor_search(distances(as.matrix(geometry[meta.data$group==sample,])), k = 11)[-1, ]
# calculate the mean expression of the nearest neighbors for all markers
# for each cell in expr
expr=expression[,meta.data$group==sample]%>% as.matrix%>% as.data.frame(check.names=F)
nnexpr <- seq_len(nrow(expr)) %>%
  map_dfr(~ expr %>%
    slice(neighbors[, .x]) %>%
    colMeans())

nn.view <- create_view("nearest", nnexpr, "nn")
misty.views %>% run_misty()

nn.view