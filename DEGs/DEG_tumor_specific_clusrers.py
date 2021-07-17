import scanpy as sc 
import numpy as np 
import pandas as pd 

adata_all = sc.read("../hb_integrated_normalized.h5ad")
df_cells = pd.read_csv("../ctype_integrated_harmony_TumorClustersAdded.txt",
                        sep = "\t", header=0, index_col=0, keep_default_na = False)
# df_cells = df_cells.loc[df_cells["Cell.class"].isin(["Tumor cell","Hepatocyte"])]
df_cells["tumor_cluster"] = ["Tr-"+str(i) for i in df_cells["tumor_cluster_03"]]
overlap_cells = np.intersect1d(adata_all.obs.index.values, df_cells.index.values)
df_cells = df_cells.loc[list(overlap_cells), :]
adata_all = adata_all[list(overlap_cells),:]
adata_all.obs = df_cells
# get either Cell class or Tumor cluster level
adata_all.obs["Cell.class_v2"] = [(adata_all.obs.iloc[i,9]) if ((adata_all.obs.iloc[i,9]) != "Tumor cell") else ((adata_all.obs.iloc[i,14])) for i in range(adata_all.shape[0])]
adata_tumor = adata_all[adata_all.obs["Cell.class"] == "Tumor cell",:]

# [‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’]
method = "t-test"

def format_DEGs(adata):
    keys = ["names","scores","logfoldchanges","pvals","pvals_adj","pts"]
    for i,key in enumerate(keys):
        a = pd.DataFrame(adata.uns["rank_genes_groups"][key]) # transfer to data frame
        b = pd.DataFrame(a.values.T.reshape(1,a.shape[0]*a.shape[1]).T) # reformat the data frame to one column
           
        if i == 0:
            b.columns = [key] # rename the column name
            b["cell group"] = sorted(list(a.columns)*a.shape[0]) # add cell group annotation
            b.set_index([key],inplace=True)
            b_merged = b
        else:
            if key in ["pts"]:
                pts_all = []
                for cell_group in np.unique(b_merged["cell group"]):
                    genes = b_merged.loc[b_merged["cell group"] == cell_group,:].index.values
                    pts_all = pts_all + list(a.loc[genes, cell_group])
                b_merged[key] = pts_all
            else:
                b_merged[key] = list(b[0])
        
    return b_merged


if __name__ == "__main__":
    # marker genes of Tr-2
    sc.tl.rank_genes_groups(adata_all, "Cell.class_v2", groups = ["Tr-2"], method = method, n_genes = adata_all.shape[1], pts = True)
    df_DEG = format_DEGs(adata_all)
    df_DEG.to_csv("tumor_clusters/DEG_Tr2_world.txt", sep="\t")

    sc.tl.rank_genes_groups(adata_tumor, "Cell.class_v2", groups = ["Tr-2"], method = method, n_genes = adata_tumor.shape[1], pts = True)
    df_DEG = format_DEGs(adata_tumor)
    df_DEG.to_csv("tumor_clusters/DEG_Tr2_tumor.txt", sep="\t")

    # marker genes of Tr-5
    sc.tl.rank_genes_groups(adata_all, "Cell.class_v2", groups = ["Tr-5"], method = method, n_genes = adata_all.shape[1], pts = True)
    df_DEG = format_DEGs(adata_all)
    df_DEG.to_csv("tumor_clusters/DEG_Tr5_world.txt", sep="\t")

    sc.tl.rank_genes_groups(adata_tumor, "Cell.class_v2", groups = ["Tr-5"], method = method, n_genes = adata_tumor.shape[1], pts = True)
    df_DEG = format_DEGs(adata_tumor)
    df_DEG.to_csv("tumor_clusters/DEG_Tr5_tumor.txt", sep="\t")

    # marker genes of Tr-3
    sc.tl.rank_genes_groups(adata_all, "Cell.class_v2", groups = ["Tr-3"], method = method, n_genes = adata_all.shape[1], pts = True)
    df_DEG = format_DEGs(adata_all)
    df_DEG.to_csv("tumor_clusters/DEG_Tr3_world.txt", sep="\t")

    sc.tl.rank_genes_groups(adata_tumor, "Cell.class_v2", groups = ["Tr-3"], method = method, n_genes = adata_tumor.shape[1], pts = True)
    df_DEG = format_DEGs(adata_tumor)
    df_DEG.to_csv("tumor_clusters/DEG_Tr3_tumor.txt", sep="\t")