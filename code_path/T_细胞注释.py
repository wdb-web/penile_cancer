import scanpy as sc
import matplotlib.pyplot as plt


adata=sc.read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
#BCL6: 编码BCL6蛋白的基因。
#CXCR5: 编码CXCR5蛋白的基因。
#PDCD1: 编码PD-1蛋白的基因。
#CXCL13: 编码CXCL13蛋白的基因。
#IL21: 编码IL-21蛋白的基因。

sc.pl.umap(adata,color=["IL21","CXCL13","PDCD1","CXCR5","BCL6",'cell_type'],vmax=0.1)
plt.savefig("/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig5/fig5_Tcell_marker.png",dpi=600)

dp = sc.pl.dotplot(adata=adata, var_names=["IL21","CXCL13","PDCD1","CXCR5","BCL6","CD4","ICOS","CCR7"],
                   standard_scale='var',
                   groupby='leiden', 
                   return_fig=True)
#sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# dp.set_xticklabels(dp.get_xticks(), rotation=45)
dp.add_totals().show()
dp.get_axes()
ax = plt.gca()
ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
plt.savefig(f'/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig5/fig5_Tcell_marker.pdf')

adata=sc.read_h5ad("output/data/adata.h5ad")
dp = sc.pl.dotplot(adata=adata, var_names=["CD74"],
                   standard_scale='var',
                   groupby='leiden', 
                   return_fig=True)
#sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# dp.set_xticklabels(dp.get_xticks(), rotation=45)
dp.add_totals().show()
dp.get_axes()
ax = plt.gca()
ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
plt.savefig(f'/home/wangdongbin/work/2024/2024_7_13_E007_TLS/output/fig5/fig5_Tcell_marker.pdf')
adata=sc.read_h5ad("output/data/adata.h5ad")


dp = sc.pl.dotplot(adata=adata, var_names=["CD74"],
                   standard_scale='var',
                   groupby='cell_type', 
                   return_fig=True)
#sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
# dp.set_xticklabels(dp.get_xticks(), rotation=45)
dp.add_totals().show()
dp.get_axes()
ax = plt.gca()
ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
plt.savefig(f'output/fig2/fig2_CD74_marker.pdf')
