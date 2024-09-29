import scanpy as sc
import matplotlib.pyplot as plt
TLSdata=sc.read_h5ad("output/data/dbscan.h5ad")
T_data=sc.read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")
T_data.obs["TLS"]=TLSdata[TLSdata.obs_names.isin(T_data.obs_names)].obs['dbscan_labels']
color={ 'CTL': '#8c564b',
 'Naive T': '#e377c2',
 'DC': '#7f7f7f',
 'NK': '#bcbd22',
 'Treg': '#17becf',
 'Proliferative T': '#1f77b4',
 'γδ T': '#ff7f0e',
 'Th1': '#2ca02c',
 'Th17': 'blue'
 }
import matplotlib.pyplot as plt
T_data=T_data[~T_data.obs['cell_subtype'].isna()]

# Create a figure with 2 subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Adjust the figsize as needed
# Plot the first UMAP
sc.pl.umap(T_data[T_data.obs["TLS"]=="-1"], color='cell_subtype', #add_outline=True,
           legend_fontsize=8, frameon=False, ax=axs[0], show=False,
           title='Clustering of cells (nTLS)', palette=color)

# Plot the second UMAP
sc.pl.umap(T_data[T_data.obs["TLS"]!="-1"], color='cell_subtype',# add_outline=True,
           legend_fontsize=8, frameon=False, ax=axs[1], show=False,
           title='Clustering of cells (TLS)', palette=color)

# Save the combined figure as a PDF
plt.savefig('/home/wangdongbin/data/2024/2024_7_13_E007_TLS/output/fig6/T_umap_sample.pdf')

# Show the plot
plt.show()
