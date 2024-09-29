library(ggplot2)
library(dplyr)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/data/project/wangdongbin/.conda/envs/mistyr/bin/python')
library(anndata)
library(dplyr)
library("Matrix")
library(ggplot2)
setwd("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS")



meta.data=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$obs

adata=read_h5ad("output/data/adata.h5ad")
gene=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$X
gene_exp_pie=data.frame(x=meta.data$x_xlide_px,y=meta.data$y_xlide_px,CD74=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$X[,"CD74"]%>%unlist)
a=ggplot()+geom_point(data=adata$obs,aes(x=x_xlide_px, y=y_xlide_px),color="#DCDCDC")+
geom_point(data=meta.data[meta.data$TLS=="Yes",],aes(x=x_xlide_px, y= y_xlide_px),color="#A9A9A9")+ theme_classic()+
geom_point(data=gene_exp_pie[gene_exp_pie$CD74%>%order,],aes(x=x,y= y,color=CD74))+
scale_colour_gradient(low = "#FBDFE2", high = "#b1182d")+
#   scale_color_gradientn(values=c("#B83945","#E3E457"))+
ggplot2::coord_fixed()

ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_cell_exp.png",a, width = 15*2, height =20*2)
ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_cell_exp.pdf",a, width = 15*2, height =20*2)



gene_exp_pie=data.frame(x=meta.data$x_xlide_px,y=meta.data$y_xlide_px,HLA_DRA=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$X[,"HLA-DRA"]%>%unlist)
a=ggplot()+geom_point(data=adata$obs,aes(x=x_xlide_px, y=y_xlide_px),color="#DCDCDC")+
geom_point(data=meta.data[meta.data$TLS=="Yes",],aes(x=x_xlide_px, y= y_xlide_px),color="#A9A9A9")+ theme_classic()+
geom_point(data=gene_exp_pie[gene_exp_pie$HLA_DRA%>%order,],aes(x=x,y= y,color=HLA_DRA))+
scale_colour_gradient(low = "#FBDFE2", high = "#b1182d")+
#   scale_color_gradientn(values=c("#B83945","#E3E457"))+
ggplot2::coord_fixed()

ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_HLA_DRA_cell_exp.png",a, width = 15*2, height =20*2)
ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_HLA_DRA_cell_exp.pdf",a, width = 15*2, height =20*2)


gene_exp_pie=data.frame(x=meta.data$x_xlide_px,y=meta.data$y_xlide_px,CD74=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$X[,"CD74"]%>%unlist)
a=ggplot()+geom_point(data=adata$obs,aes(x=x_xlide_px, y=y_xlide_px),color="#DCDCDC")+
geom_point(data=meta.data[meta.data$TLS=="Yes",],aes(x=x_xlide_px, y= y_xlide_px),color="#A9A9A9")+ theme_classic()+
geom_point(data=gene_exp_pie[gene_exp_pie$CD74%>%order,],aes(x=x,y= y,color=CD74))+
scale_colour_gradient(low = "#FBDFE2", high = "#b1182d")+
#   scale_color_gradientn(values=c("#B83945","#E3E457"))+
ggplot2::coord_fixed()

ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_cell_exp.png",a, width = 15*2, height =20*2)
ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_cell_exp.pdf",a, width = 15*2, height =20*2)


gene_exp_pie=data.frame(x=meta.data$x_xlide_px,y=meta.data$y_xlide_px,B2M=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")$X[,"B2M"]%>%unlist)
a=ggplot()+geom_point(data=adata$obs,aes(x=x_xlide_px, y=y_xlide_px),color="#DCDCDC")+
geom_point(data=meta.data[meta.data$TLS=="Yes",],aes(x=x_xlide_px, y= y_xlide_px),color="#A9A9A9")+ theme_classic()+
geom_point(data=gene_exp_pie[gene_exp_pie$B2M%>%order,],aes(x=x,y= y,color=B2M))+
scale_colour_gradient(low = "#FBDFE2", high = "#04579B")+
#   scale_color_gradientn(values=c("#B83945","#E3E457"))+
ggplot2::coord_fixed()

ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_B2M_cell_exp.png",a, width = 15*2, height =20*2)
ggsave("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/fig4_B2M_cell_exp.pdf",a, width = 15*2, height =20*2)



# install.packages("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
# 定义一个颜色渐变
col_fun <- colorRamp2(c( 0, max(gene_exp_pie$B2M)), c("#DAE3FB", "#004B9B"))
# 创建图例
legend <- Legend(
  col_fun = col_fun,  # 使用颜色渐变
  title = "Expression Level",  # 图例标题
  labels_gp = gpar(fontsize = 10),  # 图例标签的字体大小
  title_gp = gpar(fontsize = 12, fontface = "bold"),  # 图例标题的字体大小和样式
  legend_height = unit(4, "cm"),  # 图例高度
  direction = "vertical"  # 图例方向（可以是 "horizontal"）
)
pdf("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig4/color_bar.pdf")
# 绘制图例
draw(legend)
dev.off()