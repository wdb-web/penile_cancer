import scanpy as sc
import pandas as pd
import loompy as lp
import numpy as np

# cell read
cell=pd.read_csv("/home/wangdongbin/work/word/data/2024/2024_7_13_E007_TLS/output/data/Naive_T_cell_names.csv")
cell_names = cell['.']
T_cell= sc.read_h5ad("/home/wangdongbin/work/word/data/2024/2024_7_13_E007_TLS/output/data/be_T_cell_map_cells_to_space_cell_to_sp.h5ad")
T_cell_sp=T_cell.to_df().idxmax(axis=1).to_frame()

filtered_cells = T_cell_sp[T_cell_sp.iloc[:, 0].isin(cell_names)]
adata_st =sc.read_h5ad("/data/project/wangdongbin/data/2024/2024_2_3_scRNA_check_model/2024_2_3_scRNA_check_model/adata.subtype.h5ad")
adata_st=adata_st[adata_st.obs['cell_type'].isin(['T cell'])]

adata_st=adata_st[adata_st.obs_names.isin(filtered_cells.index)]

col_attrs = {"CellID": np.array(adata_st.obs_names)};
n_cells = adata_st.n_obs
row_attrs = {"Gene": np.array(adata_st.var_names),};
lp.create(f"/home/wangdongbin/work/word/data/2024/2024_7_13_E007_TLS/output/fig7/TF/adata_naive_T_Cell_SC.loom",adata_st.X.transpose(),
          row_attrs, {"CellID": np.array(adata_st.obs_names)});

Tcell=sc.read_h5ad("/data/project/tanghongzhen/data/project/E0001/nkt_adata.cell_subtype.h5ad")

Tcell=Tcell[Tcell.obs['cell_subtype'].isin(['Naive T'])]


col_attrs = {"CellID": np.array(Tcell.obs_names)};
n_cells = Tcell.n_obs
row_attrs = {"Gene": np.array(Tcell.var_names),};
lp.create(f"/home/wangdongbin/work/word/data/2024/2024_7_13_E007_TLS/output/fig7/TF/adata_naive_T_Cell_ST.loom",Tcell.X.transpose(),
          row_attrs, {"CellID": np.array(Tcell.obs_names)});

