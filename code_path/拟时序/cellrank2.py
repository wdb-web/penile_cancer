# coding='utf-8'

# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
#/data/project/tanghongzhen/software/anaconda3/envs/cellrank2/bin/python
import os
import math
import anndata as ad
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
import cellrank as cr
import warnings
import palantir
from cellrank.kernels import PseudotimeKernel
from matplotlib.backends.backend_pdf import PdfPages
from itertools import combinations
os.chdir('/home/wangdongbin/data/2024/2024_7_13_E007_TLS/')
os.makedirs('/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3', exist_ok=True)
outdir = '/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig3'
cr.settings.verbosity = 2
warnings.simplefilter("ignore", category=UserWarning)


# %%
def save_fig(name):
    plt.savefig(f'{outdir}/{name}.png', dpi=600)
    pdf = PdfPages(f'{outdir}/{name}.pdf')
    pdf.savefig()
    plt.close()
    pdf.close()


def degrees_to_vector(angle_degrees):
    angle_radians = math.radians(angle_degrees)
    dx = math.cos(angle_radians)
    dy = math.sin(angle_radians)
    return (dx, dy)


def distance_in_direction(point1, point2, angle_degrees):
    direction = degrees_to_vector(angle_degrees)
    # 计算两个点的坐标向量
    vector_point1 = [point1[0], point1[1]]
    vector_point2 = [point2[0], point2[1]]
    # 计算两个点的坐标向量在指定方向上的距离
    distance = 0
    for i in range(len(direction)):
        distance += (vector_point2[i] - vector_point1[i]) * direction[i]
    return distance


def farthest_point_at_angle(points, angle_degrees):
    # 计算中心点
    center = [sum([point[0] for point in points]) / len(points), sum([point[1] for point in points]) / len(points)]
    if angle_degrees is not None:  # 有角度则计算距离中心点最远的点
        max_distance = 0
        farthest_point = None
        for point in points:
            # 计算点与中心点的距离
            distance = distance_in_direction(center, point, angle_degrees)
            if distance > max_distance:
                max_distance = distance
                farthest_point = point
        return farthest_point
    else:  # 否则选择离中心点最近的点
        max_distance = np.inf
        farthest_point = None
        for point in points:
            point1 = center
            point2 = point
            # 计算点与中心点的距离
            distance = np.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)
            if distance < max_distance:
                max_distance = distance
                farthest_point = point
        return farthest_point


# %% 分析Patient4的肿瘤和淋巴样本 T细胞的拟时序
adata = sc.read_h5ad("/data/project/wangdongbin/data/2024//2024_7_13_E007_TLS/output/data/adata_B_cell_type_dbscan_TLS.h5ad")
adata2 = sc.read_h5ad('output/data/adata_B_TLS.h5ad')
adata.write_h5ad("/data/project/wangdongbin/data/2024//2024_7_13_E007_TLS/output/data/adata_B_cell_type_dbscan_TLS.h5ad")
# %% tcell拟时序
nkt_adata = ad.read_h5ad('output/data/adata_B_TLS.h5ad')
nkt_adata.obs['cell_subtype'] = nkt_adata.obs['cell_type_sub']
# 重新降维
# sc.tl.umap(nkt_adata)

palantir.preprocess.log_transform(nkt_adata)
dm_res = palantir.utils.run_diffusion_maps(nkt_adata, n_components=10)
ms_data = palantir.utils.determine_multiscale_space(nkt_adata)
#imputed_X = palantir.utils.run_magic_imputation(nkt_adata)

masks = palantir.presults.select_branch_cells(nkt_adata, q=.01, eps=.01)

pdf = PdfPages(f'{outdir}/Fig2-palantir-component.pdf')
palantir.plot.plot_diffusion_components(nkt_adata)
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.9)
# plt.show()
pdf.savefig()
plt.close()
pdf.close()

# 确定每种细胞的终点, 选择每个细胞类型距离坐标均值最近的细胞作为终点
cellids = []
# for cell_type, angle in zip(['Naive T', 'CTL', 'Prolif T - CD4', 'Prolif T - CD8', 'Treg', 'Exhausted T cell'], [270, 315, 90, 90, 270, 270]):

for cell_type, angle in zip(['CD74 B cells', 'IGKV4 Plasma cells', 'IGLV4 Plasma cells', 'IGHV1 Plasma cells', 'IGHM Plasma cells'], [270, 0, 90, 90, 270]):
    temp_nkt_adata = nkt_adata[nkt_adata.obs['cell_subtype'] == cell_type]
    x, y = farthest_point_at_angle(temp_nkt_adata.obsm['X_umap'], angle)
    obs_df = temp_nkt_adata.obs
    cellid = obs_df[(temp_nkt_adata.obsm['X_umap'][:, 0] == x) & (temp_nkt_adata.obsm['X_umap'][:, 1] == y)].index[0]
    cellids.append(cellid)
    # ax = sc.pl.umap(temp_nkt_adata, show=False)
    # ax.scatter(temp_nkt_adata.obsm['X_umap'][:, 0].mean(), temp_nkt_adata.obsm['X_umap'][:, 1].mean(), s=100,
    #            c='red')
    # ax.scatter(x, y, s=100, c='blue')
    # plt.show()
    # print(temp_nkt_adata.obs.loc[cellid, 'cell_subtype'])

terminal_states = nkt_adata.obs.loc[cellids, 'cell_subtype']

# terminal_states[start_cell] = 'proCAF'
palantir.plot.highlight_cells_on_umap(nkt_adata, terminal_states)
plt.tight_layout()
plt.savefig(f'{outdir}/Fig2-palantir-cell_subtype.png')

start_cell_CD74 = terminal_states[terminal_states == 'CD74 B cells'].index[0]
start_cell_cd8 = terminal_states[terminal_states == 'IGHV1 Plasma cells'].index[0]
terminal_states = terminal_states[~terminal_states.isin(['CD74 B cells', 'IGHV1 Plasma cells'])]

# %% start cell 设为 Prolif T - CD4
pr_res = palantir.core.run_palantir(nkt_adata, early_cell=start_cell_CD74, num_waypoints=500,
                                    terminal_states=terminal_states)
#pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]

palantir.plot.plot_palantir_results(nkt_adata, s=3)
plt.savefig(f'{outdir}/Fig2-plot_palantir_results.png',dpi=600)

sc.pl.umap(nkt_adata, color=['palantir_pseudotime', 'cell_subtype'], show=False, vmax=0.65)
fig = plt.gcf()
fig.set_size_inches(10, 4)
plt.tight_layout()
# plt.show()
plt.savefig(f'{outdir}/Fig2_start_B_palantir_cell_type.png',dpi=600)

masks = palantir.presults.select_branch_cells(nkt_adata, q=.01, eps=.01)
palantir.plot.plot_branch_selection(nkt_adata)
# plt.show()

palantir.plot.plot_trajectory(
    nkt_adata,
    "IGHV1 Plasma cells",
    cell_color="palantir_entropy",
    n_arrows=10,
    color="red",
    scanpy_kwargs=dict(cmap="viridis"),
    arrowprops=dict(arrowstyle="-|>,head_length=.5,head_width=.5"),
)
save_fig('Fig2-start-CD74-plot_trajectory')
more_genes = nkt_adata.var_names
gene_trends = palantir.presults.compute_gene_trends(
    nkt_adata,
    expression_key="MAGIC_imputed_data",
)

#trends = palantir.some_function_to_generate_trends(nkt_adata, "CD74 B cells")
#nkt_adata.varm['gene_trends_CD74 B cells'] = trends
gene_trends = palantir.presults.compute_gene_trends(
    nkt_adata,
    expression_key="MAGIC_imputed_data",
)
palantir.plot.plot_trend(nkt_adata, "CD74 B cells", "KLF1", color="nCount_RNA", position_layer="MAGIC_imputed_data")
plt.show()

imp_df = palantir.utils.run_magic_imputation( nkt_adata.to_df(), dm_res)

gene_trends = palantir.presults.compute_gene_trends( nkt_adata, more_genes)


communities = palantir.presults.cluster_gene_trends(nkt_adata,"CD74 B cells", more_genes)

# 假设有一个生成基因趋势的函数，如
#nkt_adata.varm['gene_trends_CD74 B cells'] = generate_gene_trends(nkt_adata, "CD74 B cells")
communities = palantir.presults.cluster_gene_trends(nkt_adata, "CD74 B cells", more_genes)
palantir.plot.plot_gene_trend_clusters(nkt_adata, "CD74 B cells")

pdf = PdfPages(f'{outdir}/Fig2-palantir-trajectory.pdf')
for cell_type in 'CD74 B cells', 'IGKV4 Plasma cells', 'IGLV4 Plasma cells', 'IGHV1 Plasma cells', 'IGHM Plasma cells':
    ax = palantir.plot.plot_trajectory(nkt_adata, cell_type, cell_color="palantir_pseudotime", n_arrows=10,
                                       color="violet", scanpy_kwargs=dict(cmap="viridis"),
                                       arrowprops=dict(arrowstyle="-|>,head_length=.5,head_width=.5"))
    # for mappable in ax.get_children():
    #     if hasattr(mappable, 'colorbar'):
    #         colorbar = mappable.colorbar
    #         colorbar.remove()
    #         break
    fig = ax.figure
    fig.set_size_inches(5, 4)
    plt.tight_layout()
    plt.show()
    pdf.savefig()
    plt.close()
pdf.close()

nkt_adata.write_h5ad('B_adata.CD74start.palantir.h5ad')

# CELLRANK2 分析驱动基因
pk = PseudotimeKernel(nkt_adata, time_key="palantir_pseudotime")
pk.compute_transition_matrix(threshold_scheme='hard', frac_to_keep=1)

# 分化方向
pk.plot_projection(basis="umap", color="palantir_pseudotime", legend_loc="on data", show=False, stream=True, s=5, alpha=1, cmap='viridis_r')
fig = plt.gcf()
ax = fig.axes[0]
left, bottom, width, height = ax.get_position().bounds
ax.set_position([left, bottom, 0.6, height])
ax = fig.axes[1]
left, bottom, width, height = ax.get_position().bounds
ax.set_position([0.85, bottom, width, height])
fig.set_size_inches(5, 4)
plt.tight_layout()
plt.savefig(f'{outdir}/Fig2_plot_projection_cell_type.png',dpi=600)

#
g = cr.estimators.GPCCA(pk)
#model = cr.models.GAMR(adata)
g.fit(cluster_key="cell_subtype", n_states=5)

g.set_terminal_states(states=['CD74 B cells', 'IGKV4 Plasma cells', 'IGLV4 Plasma cells', 'IGHV1 Plasma cells', 'IGHM Plasma cells'])
g.plot_macrostates(which="terminal", legend_loc="right", size=100, time_key='palantir_pseudotime', show=False,
                   figsize=(6, 4))
ax = plt.gca()
xticks = ax.get_xticks()
yticks = ax.get_yticks()
x_interval = (xticks[1] - xticks[0]) / 2
y_interval = (yticks[1] - yticks[0]) / 2
ax.set_xticks(np.arange(xticks[0], xticks[-1], x_interval), minor=True)
ax.set_yticks(np.arange(yticks[0], yticks[-1], y_interval), minor=True)
ax.grid(which='minor', linestyle=':', linewidth=1)
plt.tight_layout()
plt.show()
pdf = PdfPages('Fig2-cellrank-macrostates.pdf')
pdf.savefig()
plt.close()
pdf.close()

g.compute_fate_probabilities()

g.plot_fate_probabilities(same_plot=False, show=False)
plt.subplots_adjust(left=0.05, right=0.95)
pdf = PdfPages('Fig2-cellrank-fate-probabilities.pdf')
pdf.savefig()
plt.close()
pdf.close()

cr.pl.circular_projection(nkt_adata, keys=["cell_subtype"], legend_loc="right", show_edges=True,
                          label_rot='best', normalize_by_mean=True)
fig = plt.gcf()
fig.set_size_inches(4, 5)
ax = plt.gca()
ax.axis(True)
plt.title('Cell Type', y=1.05)
plt.subplots_adjust(0.05, 0.05, 0.95, 0.8333)
plt.show()

pdf = PdfPages('Fig2-cellrank-circular-projection.pdf')
pdf.savefig()
plt.close()
pdf.close()

driver_clusters = ['CD74 B cells', 'IGKV4 Plasma cells', 'IGLV4 Plasma cells', 'IGHV1 Plasma cells', 'IGHM Plasma cells']
delta_df = g.compute_lineage_drivers(
    lineages=["CD74 B cells_1"], cluster_key="cell_subtype", clusters=driver_clusters
)
delta_df.head(10)
nkt_adata.obs["fate_probabilities_CD74B1"] = g.fate_probabilities["CD74 B cells_1"].X.flatten()
sc.pl.embedding(
    nkt_adata,
    basis="umap",
    color=["fate_probabilities_CD74B1"] + list(delta_df.index[:8]),
    color_map="viridis",
    s=50,
    ncols=3,
    vmax="p96",
    show=False
)
pdf = PdfPages('Fig2-cellrank-fate-probabilities.pdf')
pdf.savefig()
plt.close()
pdf.close()

driver_df = g.compute_lineage_drivers()
driver_df.dropna(how='all', inplace=True)
nkt_adata.varm['terminal_lineage_drivers'].fillna(0, inplace=True)

pdf = PdfPages('Fig2-cellrank-lineage-drivers.pdf')
for cell_alpha, cell_beta in combinations(['CD74 B cells', 'IGKV4 Plasma cells', 'IGLV4 Plasma cells', 'IGHV1 Plasma cells', 'IGHM Plasma cells'], 2):
    # define set of genes to annotate
    alpha_genes = driver_df[f'{cell_alpha}_corr'].sort_values().index[-5:]
    beta_genes = driver_df[f'{cell_beta}_corr'].sort_values().index[-5:]
    genes_oi = {
        cell_alpha: alpha_genes,
        cell_beta: beta_genes,
    }
    # make sure all of these exist in AnnData
    assert [gene in nkt_adata.var_names for genes in genes_oi.values() for gene in genes], "Did not find all genes"
    # compute mean gene expression across all cells
    nkt_adata.var["mean expression"] = nkt_adata.X.A.mean(axis=0)
    # visualize in a scatter plot
    g.plot_lineage_drivers_correlation(
        lineage_x=cell_alpha,
        lineage_y=cell_beta,
        adjust_text=True,
        gene_sets=genes_oi,
        color="mean expression",
        legend_loc="none",
        figsize=(5, 5),
        dpi=150,
        fontsize=9,
        size=50,
    )
    pdf.savefig()
    plt.close()
pdf.close()


# cell rank 2
model = cr.models.GAMR(nkt_adata, n_knots=6, smoothing_penalty=10.0)

# compute putative drivers for the Beta trajectory
beta_drivers = g.compute_lineage_drivers(lineages="Beta")

# plot heatmap
cr.pl.heatmap(
    adata,
    model=model,  # use the model from before
    lineages="Beta",
    cluster_key="clusters",
    show_fate_probabilities=True,
    data_key="magic_imputed_data",
    genes=beta_drivers.head(40).index,
    time_key="palantir_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
)
