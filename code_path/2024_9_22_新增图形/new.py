import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
color_dict = {
"Epithelial" :"#E64B35",
"Endothelial": "#91D1C2CC",
"Macrophage": "#8491B4CC",
"Fibroblast" :"#7E6148CC",
"NK & T": "#F6D151",
"B cells": "#0F4C81",
"DC": "#60ACC0CC",
"Mast": "#D2568C"
}

adata.obs['cell_type'] = adata.obs['cell_type'].replace('Plastmass B', 'B cells')
adata.obs['cell_type'] = adata.obs['cell_type'].replace('T cell', 'NK & T')

adata=sc.read_h5ad("/home/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata_cell_type.h5ad")
sc.tl.tsne(adata,use_rep='X_pca_harmony')
sc.pl.tsne(adata,palette=color_dict,color="cell_type")
plt.subplots_adjust(left=0.1,right=0.8)  # 将left 的值调小可以向左移动内容
plt.savefig("/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig1/cell_type.png",dpi=600)

cell_types_genes = {
    'T lymphocytes': ['CD3D', 'CD3E', 'CD3G'],
    'B lymphocytes': ['CD79A', 'MS4A1', 'IGHG3'],
    'Myeloid cells': ['CD68', 'FCGR3A', 'LYZ'],
    'MAST cells': ['KIT', 'MS4A2', 'GATA2'],
    'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2'],
    'Endothelial cells': ['PECAM1', 'CLDN5',"RAMP2"]  , 
    'Epithelial cells': ['EPCAM', 'KRT19',  'KRT18']
}
color_dict = {
"Epithelial" :"#E64B35",
"Endothelial": "#91D1C2CC",
"Macrophage": "#8491B4CC",
"Fibroblast" :"#7E6148CC",
"NK & T": "#F6D151",
"B cells": "#0F4C81",
"DC": "#60ACC0CC",
"Mast": "#D2568C"
}


sorted_cell_types = list(cell_types_genes.keys())
order=['NK & T','B cells', 'Macrophage', 'Mast', 'Fibroblast','Endothelial', 'Epithelial']
# 创建一个函数来映射 cell_type 到排序的索引
def map_to_index(cell_type):
    return pd.Categorical(cell_type, categories=order, ordered=True)

# 应用这个函数到 adata.obs['cell_type'] 列
adata.obs['cell_type_ordered'] = adata.obs['cell_type'].apply(map_to_index)
adata.obs['cell_type_ordered'] = pd.Categorical(adata.obs['cell_type'], categories=order, ordered=True)

dp = sc.pl.matrixplot(adata, cell_types_genes,
                      standard_scale="var",groupby='cell_type_ordered', swap_axes=True,show=False
)
dp.add_totals().style(edge_color='black').show()

ax = plt.gca()
ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
plt.savefig(f'/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig1//marker基因_过滤.pdf')


import scanpy as sc 
import matplotlib.pyplot as plt
import pandas as pd
color_dict = {
"Cancer" :"#AC8BBF",
"Endothelial": "#91D1C2CC",
"Myeloid": "#8491B4CC",
"Fibroblast" :"#7E6148CC",
"NK & T": "#0F4C81",
"B cell": "#E64B35",
"DC": "#F6D151",
"SMC": "#D2568C"
}
adata=sc.read_h5ad("output/data/adata.h5ad")

fig,ax = plt.subplots(figsize=(6,4))
sc.pl.tsne(adata, color='cell_type',legend_loc='right margin',# add_outline=True, #
               legend_fontsize=8, frameon=False,#size=1,
               title='clustering of cells',ax=ax,show=False,palette=color_dict)
handles,labels = ax.get_legend_handles_labels()

plt.tight_layout()
plt.savefig(f'output/fig1/tsne_cell_type.png',dpi=1200)

sc.tl.tsne(adata,use_rep='X_pca_harmony', random_state=1234)
fig,ax = plt.subplots(figsize=(6,4))

sc.pl.tsne(adata, color=['cell_type'],legend_loc='right margin',# add_outline=True, #
               legend_fontsize=8, frameon=False,#size=1,
              vmax=1,show=False,palette=color_dict)

plt.tight_layout()
plt.savefig(f'output/fig1/tsne_cell_type.png',dpi=1200)


fig,ax = plt.subplots(figsize=(6,4))

sc.pl.umap(adata, color=['cell_type'],legend_loc='right margin',# add_outline=True, #
               legend_fontsize=8, frameon=False,#size=1,
              vmax=1,show=False,palette=color_dict)

plt.tight_layout()
plt.savefig(f'output/fig1/umap_cell_type.png',dpi=1200)
meta=adata.obs
plt.figure(figsize=(8, 10))
plt.scatter(meta['x_slide_mm'], meta['y_slide_mm'], c=meta["cell_type"].astype(str).map(color_dict), s=0.07, alpha=0.9, marker='o',edgecolors="none")
plt.savefig(f'output/fig1/scatter_cell_type.png',dpi=1200)
scatter = plt.scatter(meta['x_slide_mm'], meta['y_slide_mm'], 
                      c=meta["cell_type"].astype(str).map(color_dict), 
                      s=0.07, alpha=0.9, marker='o', edgecolors="none")

# 获取唯一的细胞类型
unique_types = meta['cell_type'].unique()

# 根据color_dict为每种细胞类型添加图例
handles = [plt.Line2D([0], [0], marker='o', color='w', label=cell_type, 
                      markerfacecolor=color_dict[cell_type], markersize=5) for cell_type in unique_types]

# 添加图例
plt.legend(handles=handles, title="Cell Types", loc="best", frameon=False)

# 显示图像
plt.savefig(f'output/fig1/scatter_cell_type.png',dpi=1200)
