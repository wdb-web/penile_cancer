import pandas as pd
import scanpy as sc
import decoupler as dc
import liana as li
from mudata import MuData
from anndata import AnnData

from liana.method import MistyData, genericMistyData, lrMistyData
from liana.method.sp import RandomForestModel, LinearModel, RobustLinearModel



adata_sc = sc.read_h5ad("output/data/adata_B_cell_type_dbscan_TLS.h5ad")
adata=sc.read_h5ad("output/data/adata.h5ad")
adata.obsm['spatial']=adata.obs[['x_slide_mm', 'y_slide_mm']].to_numpy(dtype='float32')
print(adata.uns['spatial'])

import numpy as np

# Convert the coordinates to a numpy array
xy_coords = adata.obs[['x_slide_mm', 'y_slide_mm']].to_numpy(dtype='float32')

# Insert the coordinates into the 'obsm' attribute under the key 'spatial'
adata.obsm['spatial'] = xy_coords
adata2 = AnnData(adata.X, obsm={"spatial": xy_coords},obs=adata.obs)
adata.uns['spatial'][library_id]
plot, _ = li.ut.query_bandwidth(coordinates=adata2.obsm['spatial'], start=0, end=500, interval_n=20)
plot, _ = li.mt.query_bandwidth(coordinates=adata2.obsm['spatial'], start=0, end=500, interval_n=20)
li.ut.spatial_neighbors(adata2, bandwidth=200, cutoff=0.1, kernel='gaussian', set_diag=True)
li.mt.bivariate(adata2,
                resource_name='consensus', # NOTE: uses HUMAN gene symbols!
                local_name='cosine', # Name of the function
                global_name="morans", # Name global function
                n_perms=100, # Number of permutations to calculate a p-value
                mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                add_categories=True, # Whether to add local categories to the results
                nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True
                )
li.method.bivariate(adata,
                resource_name='consensus', # NOTE: uses HUMAN gene symbols!
                local_name='cosine', # Name of the function
                global_name="morans", # Name global function
                n_perms=100, # Number of permutations to calculate a p-value
                mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                add_categories=True, # Whether to add local categories to the results
                nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True
                )
liana_results = li.liana_wrap(adata, method='connectome')
