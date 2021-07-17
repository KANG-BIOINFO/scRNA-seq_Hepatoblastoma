import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import loompy
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo',dpi=200)  # for beautified visualization

# load 
df_cells = pd.read_csv("/data/aronow/Kang/SingleCellHBPilot-HB17/integrate/tumor integration/ctype_tumors_v2_cluster6_separated.txt",
                        sep="\t", header=0, index_col=0, keep_default_na = False)
df_remove = pd.read_csv("/data/aronow/Kang/SingleCellHBPilot-HB17/integrate/tumor integration/ctype_remove.csv",
                        sep=",", header=0, index_col=0, keep_default_na = False)
remaining_cells = list(df_remove.loc[df_remove["remove"] == "F",:].index.values) # remove tumor cells from Tr-1, Tr-4, Tr-6, Tr-8 (with low sequencing depth) to avoid bias to the velocity.
df_cells = df_cells.loc[remaining_cells,:]        

df_umap = pd.read_csv("/data/aronow/Kang/SingleCellHBPilot-HB17/integrate/tumor integration/umap_tumors.txt",
                      sep="\t", header=0, index_col=0, keep_default_na = False)
df_umap = df_umap.loc[list(df_cells.index.values),:]

adata_merged = scv.read("hb_merged_loom_raw.h5ad")
adata_merged = adata_merged[list(df_cells.index.values),:] # extract tumor cells
for col in df_cells.columns:
    if (col not in adata_merged.obs.columns):
        adata_merged.obs[col] = df_cells[col]
adata_merged.obsm["X_umap"] = df_umap.values

# calculate 
scv.pp.filter_and_normalize(adata_merged)
scv.pp.moments(adata_merged)
adata_merged.write("hb_velocity_tumor_removed_v2_a.h5ad")
scv.tl.velocity(adata_merged, mode='stochastic')
scv.tl.velocity_graph(adata_merged)
scv.tl.velocity_pseudotime(adata_merged)
adata_merged.write("hb_velocity_tumor_removed_v2_b.h5ad")

# visualization
# colormap = ["#F8766D","#DE8B00","#B79F01","#7CAE00","#02BA38","#00C08B","#00BFC4","#00B4F0","#619CFF","#C475FF","#F564E3","#FF62AF"] # all clusters
colormap = ["#F8766D","#DE8B00","#B79F01","#7CAE00","#00C08B","#00BFC4","#00B4F0","#619CFF","#C475FF","#F564E3","#FF62AF"] # without c4
colormap = ["#F8766D","#B79F01","#7CAE00","#00C08B","#00B4F0","#619CFF","#C475FF","#F564E3","#FF62AF"] # remove c1,4,6

scv.pl.velocity_embedding(adata, color = "integrated_snn_res.0.3", save = "tumor_velocity_remove_v2_colorFormatted.png", palette = colormap, dpi = 96)
scv.pl.velocity_embedding_grid(adata, color = "integrated_snn_res.0.3", save = "tumor_velocity_grid_remove_v2_colorFormatted.png", palette = colormap, dpi = 96)
scv.pl.velocity_embedding_stream(adata, color = "integrated_snn_res.0.3", save = "tumor_velocity_stream_remove_v2_colorFormatted.png", palette = colormap, dpi = 96)

