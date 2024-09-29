import scanpy as sc
import os

from scipy.sparse import coo_matrix
from sklearn.neighbors import NearestNeighbors

import anndata as ad
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import scanpy as sc
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
adata=sc.read_h5ad("/data/project/tanghongzhen/data/project/E0001/adata.cellanno.h5ad")
os.chdir('/home/wangdongbin/work/2024/2024_7_13_E007_TLS')

# %%

sample_info = pd.read_excel('/data/project/tanghongzhen/data/project/E0001/sample_info.xlsx')
sample_info.index = sample_info['sample_id']
sample_info.index.name = None

patient_info = pd.read_excel('/data/project/tanghongzhen/data/project/E0001/patient_info.xlsx', index_col=None)
patient_info.index = patient_info['patient_id']
patient_info.index.name = None

# 整理样本信息并匹配到adata中
# patient_info按照样本名匹配到sample_info
sample_info[patient_info.columns] = patient_info.loc[sample_info['patient_id']].to_numpy()

# sample_info按照样本名匹配到adata
# adata.obs['sample'] = 'S' + adata.obs['sample'].astype(str)
adata.obs[sample_info.columns] = sample_info.loc[adata.obs['sample']].astype(str).to_numpy()
adata_B=adata[adata.obs['cell_type']=="B cell"].copy()
#sc.pp.highly_variable_genes(adata_B, n_top_genes=1000, batch_key='sample')
#sc.pp.pca(adata_B, use_highly_variable=True)
#sc.external.pp.harmony_integrate(adata_B, key=['sample'], max_iter_harmony=30)
#sc.pp.neighbors(adata_B, use_rep='X_pca_harmony', n_neighbors=30)
#sc.tl.umap(adata_B)
#sc.tl.leiden(adata_B, resolution=4)

#adata.write_h5ad("adata.cellanno_Bcell.h5ad")

#sc.pl.umap(adata_B,color="leiden")
#sc.pl.umap(adata_B,color="HPV")

n_neighbors = 201
sample = 'sample_id'
resolution = 0.1
# X，其中每一行代表一个细胞的x和y
# X=data.obs[['x_slide_mm', 'y_slide_mm']]
# 创建 NearestNeighbors 模型，并拟合特征矩阵 X
samples = adata.obs['sample'].unique()
raw_n_neighbors = n_neighbors
sparse_count_df = pd.DataFrame()
for sample in samples:
    sample_adata = adata[adata.obs['sample'] == sample]
    # if raw_n_neighbors == 'auto':
    #     n_neighbors = int(len(sample_adata) * 0.1)
    #     n_neighbors = max(n_neighbors, 30)
    #     n_neighbors = min(n_neighbors, 2000)
    #     n_neighbors = min(n_neighbors, len(sample_adata))
    print(n_neighbors)
    nn = NearestNeighbors(n_neighbors=n_neighbors)
    nn.fit(sample_adata.obs[['x_slide_mm', 'y_slide_mm']])
    distances, indices = nn.kneighbors(sample_adata.obs[['x_slide_mm', 'y_slide_mm']])
    indices = indices[:, :n_neighbors]  # 最近的n个细胞可能会因为距离相同而超过n个,所以取前n个
    indices = indices.flatten()

    # 获取所有最近邻细胞的细胞类型
    celltype_indices = sample_adata.obs['cell_type'].iloc[indices]
    sample_sparse_count_df = pd.DataFrame({'cell_type': celltype_indices})
    sample_sparse_count_df['cell'] = pd.Categorical(
        np.array([[i] * n_neighbors for i in sample_adata.obs_names]).flatten())
    sparse_count_df = pd.concat([sparse_count_df, sample_sparse_count_df], axis=0)

sparse_count_df = sparse_count_df.groupby(['cell_type', 'cell']).size().reset_index(name='Count')
count_df = sparse_count_df.pivot_table(index='cell', columns='cell_type', values='Count', aggfunc='sum',
                                       fill_value=0)
count_df=count_df[count_df.index.isin(adata_B.obs_names)]
sum(count_df[["B cell","NK & T","DC"]].T.sum() >120)


# 显示所有的出现次数

count_df[["B cell","NK & T","DC"]].T.sum()
plt.figure(figsize=(20, 6))
count_df[["B cell","NK & T","DC"]].T.sum().value_counts().sort_index().plot(kind='bar')
plt.title('Value Counts of Row Sums')
plt.xlabel('Sum Value')
plt.ylabel('Frequency')
plt.xticks(rotation=90, fontsize=4)
plt.savefig("Sum Value.png")

cellcount_adata = ad.AnnData(X=count_df)
cellcount_adata=cellcount_adata[cellcount_adata.obs_names.isin(adata_B.obs_names)]


# 假设 adata 是你的 AnnData 对象

all_adata=sc.read_h5ad("/data/project/tanghongzhen/data/project/E0001/adata.cellanno.h5ad")
adata.obs["x"]=all_adata.obs["x_xlide_px"]
adata.obs["y"]=all_adata.obs["y_xlide_px"]
adata_B=all_adata[all_adata.obs['cell_type']=="B cell"].copy()
adata_B.obs["center"]=count_df[["B cell","NK & T","DC"]].T.sum()
plt.figure(figsize=(1.5*15, 15))
B_cell=adata_B[adata_B.obs['center'].isin(list(range(100,200)))].obs
plt.scatter(B_cell['x_slide_mm'], B_cell['y_slide_mm'], c="red", s=0.04, alpha=0.5, marker='o')
plt.show()

leiden_counts = adata.obs['leiden'].value_counts()

# 只展示前30个Leiden群集的计数
top_30_leiden_counts = leiden_counts.head(30)
top_30_leiden_indices = top_30_leiden_counts.index

top_30_adata = adata#[adata.obs['leiden'].isin(top_30_leiden_indices)]
df=top_30_adata.obs

# 使用seaborn调色板自动生成40种颜色
palette = sns.color_palette("hsv", n_colors=30)
leiden_colors = {i: palette[i] for i in range(30)}
df = df[df['leiden'].astype(int).isin(leiden_colors.keys())]
# 为每个Leiden群集分配颜色
# 绘制散点图
# 计算 x 和 y 的范围
x_range = df['x'].max() - df['x'].min()
y_range = df['y'].max() - df['y'].min()
df['leiden'].value_counts()
#35 "14","15","16","17","18","19","20","22","26" list(range(0,5))
df_B=  df[df.index.isin(adata_B.obs_names) & df['leiden'].astype(int).isin(
   list(range(10,40)))
]

# 根据 x 和 y 的范围计算长宽比例
aspect_ratio = x_range / y_range
plt.figure(figsize=(aspect_ratio*15, 15))
colors = df_B['leiden'].astype(int).map(leiden_colors)
plt.scatter(df_B['x'], df_B['y'], c="red", s=0.04, alpha=0.5, marker='o')
plt.show()




plt.savefig("output/fig2/fig2_领域分析.png",dpi=600)
# 添加图例


adata_B=sc.read_h5ad("adata.cellanno_Bcell.h5ad")
color_dict = {
"Cancer" :"#E64B35",
"Endothelial": "#91D1C2CC",
"Myeloid": "#8491B4CC",
"Fibroblast" :"#7E6148CC",
"NK & T": "#F6D151",
"B cell": "#0F4C81",
"DC": "#60ACC0CC",
"SMC": "#D2568C"
}

plt.figure(figsize=(aspect_ratio*15, 15))
colors = all_adata.obs["cell_type"].astype(str).map(color_dict)
plt.scatter(all_adata.obs['x_slide_mm'], all_adata.obs['y_slide_mm'],
            c=colors, s=0.1, alpha=1, marker='o',edgecolors="none")
plt.savefig("cell_type.png",dpi=1200)



import matplotlib.pyplot as plt
import pandas as pd

# 定义感兴趣的基因列表及其对应的名称
genes_of_interest = {
"CD19":"CD19",
"CD22":"CD22",
    "MS4A1": "CD20",
    "CR2": "CD21",
    "CD68": "CD68",
    "FCER2": "CD23",
    "CD274": "PD-L1",
    "MKI67": "Ki67",
    "GLYCAM1": "PNAd",
    "CXCL13": "CXCL13",
    "PDPN": "Podoplanin",
    "BCL6": "Bcl-6"
}
# 定义感兴趣的基因列表及其对应的名称
genes_of_interest = {
"CD19":"CD19",
"CD22":"CD22",
    "MS4A1": "CD20",
    "CR2": "CD21",
    "CD68": "CD68",
    "FCER2": "CD23",
    "CD274": "PD-L1",
    "MKI67": "Ki67",
    "GLYCAM1": "PNAd",
    "CXCL13": "CXCL13",
    "PDPN": "Podoplanin",
    "BCL6": "Bcl-6"
}
# 设置绘图的宽高比
aspect_ratio = 1.5
color_dict = {i: plt.cm.viridis(i / 10.0) for i in range(10)}
import matplotlib.colors.Colormap as color
from matplotlib.colors import LinearSegmentedColormap
custom_cmap = LinearSegmentedColormap.from_list("mycmap", ["white", "#004D40"])

cmap = cm.get_cmap('inferno')
#inferno magma
custom_cmap = LinearSegmentedColormap.from_list("mycmap", ["white", "#004D40"])
adata_st=all_adata
adata_st.obsm['spatial'] = adata_st.obs[['x_slide_mm', 'y_slide_mm']]




cmap = cm.get_cmap('inferno')
#inferno magma
custom_cmap = LinearSegmentedColormap.from_list("mycmap", ["white", "#004D40"])
gene="CD3E"
cmap = plt.cm.get_cmap('inferno')
# inferno magma
gene_expression = all_adata[:, gene].to_df()
custom_cmap = LinearSegmentedColormap.from_list("mycmap", ["white", "#004D40"])
custom_cmap(gene_expression)
spatial_expression(adata_st, gene, cmap)
"CD19"
genes_of_interest = {
"CD19":"CD19",
"CD22":"CD22",
"MS4A1": "CD20",
"SDC1":"CD138",
#"IGHM":"IgM",
#"IGHD":"IgD",
#"IGHA":"IgA",
#"IGHE":"IgE",

"CR2":"CD21" ,
"FCER2": "CD23",
"MKI67":"Ki67",
"BCL6":"BCL6",
"AICDA":"AICDA"
}

gene="CD21";protein="CR2"
for gene, protein in genes_of_interest.items():
    # 获取感兴趣基因的表达数据
    gene_expression = all_adata[:, gene].to_df()

    # 将基因表达数据进行分级并映射到颜色
    norm = plt.Normalize(vmin=gene_expression.min(), vmax=gene_expression.max())
    expression_colors = cmap(norm(gene_expression))
    size = np.where(gene_expression>0,0.5,0.1)
    # 绘制基因表达的原位图像
    plt.figure(figsize=( 8, 10))
    #colors = all_adata.obs["cell_type"].astype(str).map(color_dict)
    #plt.scatter(all_adata.obs['x_slide_mm'], all_adata.obs['y_slide_mm'],
    #            c=colors, s=0.1, alpha=1, marker='o', edgecolors="none")
    plt.scatter(all_adata.obs['x_slide_mm'], all_adata.obs['y_slide_mm'],
                c=expression_colors, s=size   , alpha=1, marker='o', edgecolors="none")

    # 设置标题和保存图片
    plt.title(f"In situ expression of {protein}")
    plt.savefig(f"{protein}_expression.png", dpi=1200)
    plt.close()


plt.scatter(all_adata.obs['x_slide_mm'], all_adata.obs['y_slide_mm'],
                c=expression_colors, s=0.1, alpha=1, marker='o', edgecolors="none")
plt.show()

