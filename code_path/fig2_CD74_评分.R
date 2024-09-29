
library(dplyr)
library(Seurat)
#SeuratObject
#remotes::install_cran("Matrix", version = "1.5.3", repos = "https://cran.rstudio.com/")
reticulate::use_python('/data/project/wangdongbin/.conda/envs/mistyr/bin/python')
library(anndata)
library(dplyr)
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
library(scatterpie)
library(clusterProfiler)
library(enrichplot)
library(AUCell)
library(ggpubr)
library(reticulate)
#https://github.com/Moonerss/CIBERSORT
library(CIBERSORT)
library(GSVA)
#library(ImmuCellAI)
#library(xCell)
library(Seurat)

#Sys.setenv(PATH = paste("/data/project/wangdongbin/.miniconda/bin", Sys.getenv("PATH"), sep = ":"))
library(reticulate)
mixture_data=read.csv("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/raw_data/RNA_OS.csv") %>% na.omit
rownames(mixture_data)=mixture_data$X
mixture_matrix=mixture_data[,-c(1:4)] %>%t %>% as.matrix
#devtools::install_local("/home/wangdongbin/data/data/r-immunedeconv-2.1.2-r43hdfd78af_2.tar.bz2")
future::plan("multisession", workers =20)  # Reduce the number of workers
options(future.globals.maxSize = 100 * 1024^3) # 将最大值设为600 MiB
Tcell=read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")
meta.data<- Tcell$obs %>%
  dplyr::mutate(cell_annotation = dplyr::case_when(
    leiden == '0' ~ 'Memory CTL',
    leiden == '1' ~ 'Proliferating CTL',
    leiden == '2' ~ 'Stressed CTL',
    leiden == '3' ~ 'ANXA1+ CTL',
    TRUE ~ 'Unknown'  # Optional: catch-all for any values not specified
  ))
exp=Tcell$X %>% as.matrix
sig_matrix=data.frame(cell_type= meta.data[,'cell_annotation'],exp)
#　CIBERSORT
sig_matrix_mean=sig_matrix%>%group_by(cell_type) %>% summarise_if(is.numeric,mean) %>% data.frame(chech.name=FALSE)
rownames(sig_matrix_mean)=sig_matrix_mean[,1]
sig_matrix_mean=sig_matrix_mean[,-1] %>%t
sig_matrix_mean_matrix <- sig_matrix_mean %>% as.matrix

#sig_matrix, mixture_file, 
cibersort_results <- cibersort(sig_matrix=sig_matrix_mean_matrix, mixture_file=mixture_matrix, QN = FALSE)

# ssgsea 类型
x=t(Tcell$X) %>%as.matrix
Seurat_Object=CreateSeuratObject(x,meta.data=meta.data%>% data.frame)
Seurat_Object <- NormalizeData(Seurat_Object)

Idents(Seurat_Object)=meta.data$cell_annotation

marker=Seurat::FindAllMarkers(Seurat_Object)

marker_list=marker[marker$avg_log2FC>0.2 &marker$p_val_adj <0.01, ] %>% group_by(cluster) %>% 
    dplyr::arrange(avg_log2FC) %>% top_n(n=10,wt=avg_log2FC) %>%{split(.$gene ,.$cluster)}  #group_by(cluster) %>% group_split()
ssgsea_results <- gsva(mixture_matrix, marker_list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

# AUCell 类型

cells_rankings <- AUCell_buildRankings(mixture_matrix, nCores=5, plotStats=TRUE)


cells_AUC <- AUCell_calcAUC(marker_list, cells_rankings,nCores = 5, aucMaxRank=nrow(cells_rankings)*0.1)



writexl::write_xlsx(list(cibersort_results=cibersort_results %>% as.data.frame %>% tibble::rownames_to_column("sample"),
ssgsea_results=ssgsea_results %>% as.data.frame%>% tibble::rownames_to_column("sample"),
AUCell=AUCell::getAUC(cells_AUC) %>% as.matrix %>% as.data.frame%>% tibble::rownames_to_column("sample")),
"/home/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/output/CTL_score/score.xlsx")


library(ImmuCellAI)

# 假设你已经加载了ImmuCellAI包
# your_marker_list 是一个包含你新的marker基因的列表
cor_list=sig_matrix_mean_matrix[,-1] %>% apply(2 ,as.numeric) %>% t%>% na.omit  %>% cor  %>% data.frame
colnames(cor_list)=sig_matrix_mean_matrix[,1]
rownames(cor_list)=sig_matrix_mean_matrix[,1]
marker_exp_data=sig_matrix_mean_matrix[,-1] %>% apply(2 ,as.numeric) %>% t%>% na.omit
colnames(marker_exp_data)=sig_matrix_mean_matrix[,1]
result <- ImmuCellAI::ImmuCellAI_new(sample = mixture_matrixm,data_type, group_tag = marker_list)
ImmuCellAI::ImmuCellAI_new
function (sample, data_type, group_tag, response_tag, customer = 0, 
    sig_file = NULL, exp_file = NULL) 
{
    data("marker_exp")
    data("paper_marker")
    data("compensation_matrix")
    group_fre <<- c()
    ICB_response <<- c()
    layer1 = c(names(marker_list))

    paper_marker <<- marker_list
    compensation_matrix <<- cor_list[layer1, layer1]
    marker_exp <<- marker_exp[, layer1]
    layer1_abun <- Sample_abundance_calculation(sample, paper_marker, 
        marker_exp, data_type, customer)
    layer2 = c("CD4_naive", "CD8_naive", "Cytotoxic", "Exhausted", 
        "Tr1", "nTreg", "iTreg", "Th1", "Th2", "Th17", "Tfh", 
        "Central_memory", "Effector_memory", "MAIT")
    data("marker_exp")
    data("paper_marker")
    data("compensation_matrix")
    paper_marker <<- paper_marker[layer2]
    compensation_matrix <<- compensation_matrix[layer2, layer2]
    marker_exp <<- marker_exp[, layer2]
    layer2_abun <- Sample_abundance_calculation(sample, paper_marker, 
        marker_exp, data_type, customer)
    layer1_abun_normlized <- apply(layer1_abun[, -11], 1, function(x) (x/sum(x)))
    layer1_abun_normlized <- rbind(layer1_abun_normlized, layer1_abun[, 
        11])
    row.names(layer1_abun_normlized)[11] = colnames(layer1_abun)[11]
    if (group_tag) {
        group_content <<- sample[1, ]
    }
    layer2_abun_normlized <- apply(layer2_abun, 1, function(x) (x/sum(x)))
    CD4_sub <- c("CD4_naive", "Tr1", "nTreg", "iTreg", "Th1", 
        "Th2", "Th17", "Tfh")
    CD8_sub <- c("CD8_naive", "Cytotoxic", "Exhausted", "MAIT")
    CD4_sub_all <- c()
    for (i in seq(1, ncol(layer2_abun_normlized[CD4_sub, ]))) {
        CD4_sub_all <- cbind(CD4_sub_all, layer2_abun_normlized[CD4_sub, 
            i] * layer1_abun_normlized["CD4_T", i])
    }
    CD8_sub_all <- c()
    for (i in seq(1, ncol(layer2_abun_normlized[CD8_sub, ]))) {
        CD8_sub_all <- cbind(CD8_sub_all, layer2_abun_normlized[CD8_sub, 
            i] * layer1_abun_normlized["CD8_T", i])
    }
    Central_memory <- layer2_abun_normlized["Central_memory", 
        ] * (layer1_abun_normlized["CD8_T", ] + layer1_abun_normlized["CD4_T", 
        ])
    Effector_memory <- layer2_abun_normlized["Effector_memory", 
        ] * (layer1_abun_normlized["CD8_T", ] + layer1_abun_normlized["CD4_T", 
        ])
    all_norm = rbind(layer1_abun_normlized, CD4_sub_all, CD8_sub_all, 
        Central_memory, Effector_memory)
    result_layer(t(round(all_norm, 3)), group_tag, response_tag)
    return(list(Sample_abundance = T_FRE, Group_result = group_fre, 
        Response = ICB_response))
}