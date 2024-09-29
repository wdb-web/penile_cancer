import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
import tangram as tgw
import matplotlib.pyplot as plt

import anndata as ad


adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")

adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata_B_cell_type_dbscan_TLS.h5ad")
adata_st=adata_st[adata_st.obs['cell_type'].isin(['B cells','Plasma B'])]
#sc.pp.scale(adata_sc)
#sc.pp.scale(adata_st)
tg.pp_adatas(adata_st,adata_sc )
ad_map = tg.map_cells_to_space(adata_st,adata_sc ,

    #mode="cells",
    mode="clusters",
    cluster_label='cell_type_sub',  # .obs field w cell types
    density_prior='rna_count_based',
    num_epochs=500,
     device="cuda:0",
    #device='cpu',
)
ad_map.write_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/map_cells_to_space.h5ad")




ad_map2 = tg.map_cells_to_space(adata_st,adata_sc,
    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
     device="cuda:0",
    #device='cpu',
)
ad_map2.write_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/B_cell_map_cells_to_space_cell_to_sp.h5ad")



adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")

adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata.h5ad")

adata_st=adata_st[adata_st.obs['cell_type'].isin(['T cell'])]
 
adata_sc=adata_sc[adata_sc.obs['cell_type'].isin(['NK & T'])]
tg.pp_adatas(adata_st,adata_sc )

ad_map2 = tg.map_cells_to_space(adata_st,adata_sc,

    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
     device="cuda:0",
    #device='cpu',
)
ad_map2.write_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/T_cell_map_cells_to_space_cell_to_sp.h5ad")










#Cancer

adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")

adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata.h5ad")

adata_st=adata_st[adata_st.obs['cell_type'].isin(['Epithelial'])]
 
adata_sc=adata_sc[adata_sc.obs['cell_type'].isin(['Cancer'])]
tg.pp_adatas(adata_st,adata_sc )

ad_map2 = tg.map_cells_to_space(adata_st,adata_sc,

    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
     device="cuda:1",
   # device='cpu',
)
ad_map2.write_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Cancer_cell_map_cells_to_space_cell_to_sp.h5ad")



# 内皮



adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")

adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata.h5ad")

adata_st=adata_st[adata_st.obs['cell_type'].isin(['Endothelial'])]
 
adata_sc=adata_sc[adata_sc.obs['cell_type'].isin(['Endothelial'])]
tg.pp_adatas(adata_st,adata_sc )

ad_map2 = tg.map_cells_to_space(adata_st,adata_sc,

    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
     device="cuda:1",
    #device='cpu',
)
ad_map2.write_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Endothelial_cell_map_cells_to_space_cell_to_sp.h5ad")


# Fibroblast

import torch

adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")
adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata.h5ad")
adata_st=adata_st[adata_st.obs['cell_type'].isin(['Fibroblast'])]
adata_sc=adata_sc[adata_sc.obs['cell_type'].isin(['Fibroblast'])]

for sample in adata_st.obs['sample'].unique():
    adata_obj = adata_st[adata_st.obs['sample'] == sample]
    # 先尝试在 GPU 上运行
    try:
        # 预处理数据
        tg.pp_adatas(adata_obj, adata_sc)
        # 在 GPU 上执行 map_cells_to_space
        ad_map2 = tg.map_cells_to_space(
            adata_obj, adata_sc,
            mode="cells",
            density_prior='rna_count_based',
            num_epochs=500,
            device="cuda:1"
        )
    except RuntimeError as e:
        # 如果 CUDA 内存不足或其他错误发生，捕获异常并切换到 CPU
        print(f"CUDA 内存不足或其他错误，切换到 CPU 进行计算: {e}")
            # 在 CPU 上重新执行 map_cells_to_space
        ad_map2 = tg.map_cells_to_space(
                adata_obj, adata_sc,
                mode="cells",
                density_prior='rna_count_based',
                num_epochs=500,
                device="cpu"
            )
    # 保存结果
    ad_map2.write_h5ad(f"/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/{sample}_Fibroblast_cell_map_cells_to_space_cell_to_sp.h5ad")
    # 释放内存
    del ad_map2
    torch.cuda.empty_cache()
#ad_map2.write_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Fibroblast_cell_map_cells_to_space_cell_to_sp.h5ad")
# Myeloid

adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")

adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata.h5ad")

adata_st=adata_st[adata_st.obs['cell_type'].isin(['Myeloid','Mast'])]
 
adata_sc=adata_sc[adata_sc.obs['cell_type'].isin(['Myeloid','DC'])]

for sample in adata_st.obs['sample'].unique() :
    adata_obj=adata_st[adata_st.obs['sample']==sample]
    tg.pp_adatas(adata_obj,adata_sc )
    ad_map2 = tg.map_cells_to_space(adata_obj,adata_sc,
    mode="cells",
    density_prior='rna_count_based',
    num_epochs=500,
     device="cuda:1",
    #device='cpu',
    )
    ad_map2.write_h5ad(f"/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/{sample}_Myeloid_DC_cell_map_cells_to_space_cell_to_sp.h5ad")
    del ad_map2
    torch.cuda.empty_cache()


# 在调用 map_cells_to_space 之前或之后释放显存


# 定义文件夹路径和文件列表
data_folder = "/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/"
file_list = [f for f in os.listdir(data_folder) if f.endswith('_Myeloid_cell_map_cells_to_space_cell_to_sp.h5ad')]

# 读取并合并所有 .h5ad 文件
adatas = []
for file_name in file_list:
    adata = ad.read_h5ad(os.path.join(data_folder, file_name))
    adatas.append(adata)
    

adata_combined = tg.map_cells_to_space(adatas[0], adatas[0], mode="cells", density_prior='rna_count_based', num_epochs=0)

# 合并成一个大的 AnnData 对象
#adata_merged = ad.concat(adatas, axis=0)
adata_merged = ad.concat(adatas, axis=0, join='outer', label='sample', keys=[f.split('_')[0] for f in file_list])

# 保存合并后的 .h5ad 文件
output_file = "/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/Myeloid_cell_map_cells_to_space_merged.h5ad"
adata_merged.write_h5ad(output_file)

print(f"合并后的数据已保存到 {output_file}")



# Cancer
adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")

adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata.h5ad")

adata_st=adata_st[adata_st.obs['cell_type'].isin(['Epithelial'])]
 
adata_sc=adata_sc[adata_sc.obs['cell_type'].isin(['Cancer'])]
import numpy as np
import scanpy as sc

# 设置随机种子
np.random.seed(42)

# 确保在随机抽样之前所有细胞都被打乱
adata_sc = adata_sc[np.random.permutation(adata_sc.n_obs), :]

# 随机选择 100,000 个细胞
adata_sampled = adata_sc[np.random.choice(adata_sc.n_obs, 100000, replace=False), :]

# 输出结果
#adata_sampled=adata_sc

for sample in adata_st.obs['sample'].unique():
    adata_obj = adata_st[adata_st.obs['sample'] == sample]
    # 先尝试在 GPU 上运行
    try:
        # 预处理数据
        tg.pp_adatas(adata_obj, adata_sc)
        # 在 GPU 上执行 map_cells_to_space
        ad_map2 = tg.map_cells_to_space(
            adata_obj, adata_sc,
            mode="cells",
            density_prior='rna_count_based',
            num_epochs=500,
            device="cuda:1"
        )
    except RuntimeError as e:
        # 如果 CUDA 内存不足或其他错误发生，捕获异常并切换到 CPU
        print(f"CUDA 内存不足或其他错误，切换到 CPU 进行计算: {e}")
            # 在 CPU 上重新执行 map_cells_to_space
        ad_map2 = tg.map_cells_to_space(
                adata_obj, adata_sc,
                mode="cells",
                density_prior='rna_count_based',
                num_epochs=500,
                device="cpu"
            )
    # 保存结果
    ad_map2.write_h5ad(f"/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/{sample}_Cancer_cell_map_cells_to_space_cell_to_sp.h5ad")
    # 释放内存
    del ad_map2
# tg 插补


adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")

adata_sc = sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/adata_B_cell_type_dbscan_TLS.h5ad")
adata_st=adata_st[adata_st.obs['cell_type'].isin(['B cells','Plasma B'])]
#sc.pp.scale(adata_sc)
#sc.pp.scale(adata_st)

ad_map=sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_7_13_E007_TLS/output/data/be_B_cell_map_cells_to_space_cell_to_sp.h5ad")



ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_st)
genes = ['rgs6', 'satb2',  'cdh12']
genes=['mrgprx2', 'muc20', 'chrna2']
tg.plot_genes_sc(genes, adata_measured=ad_sp, adata_predicted=ad_ge, spot_size=50, scale_factor=0.1, perc=0.02, return_figure=False)