import scanpy as sc 
import numpy as np 
import pandas as pd 

adata_all = sc.read("../hb_integrated_normalized.h5ad")
df_cells = pd.read_csv("../ctype_integrated_harmony_TumorClustersAdded.txt",
                        sep = "\t", header=0, index_col=0, keep_default_na = False)
df_cells = df_cells.loc[df_cells["Cell.class"].isin(["Tumor cell","Hepatocyte"])]
overlap_cells = np.intersect1d(adata_all.obs.index.values, df_cells.index.values)
df_cells = df_cells.loc[list(overlap_cells), :]
adata_all = adata_all[list(overlap_cells),:]
adata_all.obs = df_cells
adata_all.obs["Cell.class_v2"] = [(adata_all.obs.iloc[i,11] + "-" + adata_all.obs.iloc[i,9]) for i in range(adata_all.shape[0])]

# [‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’]
method = "t-test"

def format_DEGs(adata):
    keys = ["names","scores","logfoldchanges","pvals","pvals_adj","pts","pts_rest"]
    for i,key in enumerate(keys):
        a = pd.DataFrame(adata.uns["rank_genes_groups"][key]) # transfer to data frame
        b = pd.DataFrame(a.values.T.reshape(1,a.shape[0]*a.shape[1]).T) # reformat the data frame to one column
           
        if i == 0:
            b.columns = [key] # rename the column name
            b["cell group"] = sorted(list(a.columns)*a.shape[0]) # add cell group annotation
            b.set_index([key],inplace=True)
            b_merged = b
        else:
            if key in ["pts","pts_rest"]:
                pts_all = []
                for cell_group in np.unique(b_merged["cell group"]):
                    genes = b_merged.loc[b_merged["cell group"] == cell_group,:].index.values
                    pts_all = pts_all + list(a.loc[genes, cell_group])
                b_merged[key] = pts_all
            else:
                b_merged[key] = list(b[0])
        
    return b_merged

def getCellGroupDEG(adata, positive = "tumor-Tumor cell", negative = "background-Hepatocyte"):
    print(adata.shape)
    adata = adata[adata.obs["Cell.class_v2"].isin([positive, negative]), :]
    print(adata.shape)
    sc.tl.rank_genes_groups(adata, "Cell.class_v2", method = method, n_genes = adata.shape[1], pts = True)

    df_DEG = format_DEGs(adata)
    df_DEG = df_DEG.loc[df_DEG["cell group"] == positive,:]
    
    positive = positive.replace("/","_")
    negative = negative.replace("/","_")
    name = "./Fig4/" + positive + " vs. " + negative + ".txt"

    print("name")
    print(df_DEG.shape)
    df_DEG.to_csv(name, sep="\t")

if __name__ == "__main__":
    getCellGroupDEG(adata_all.copy(), positive= "tumor-Tumor cell", negative = "background-Hepatocyte")
    getCellGroupDEG(adata_all.copy(), positive= "pdx-Tumor cell", negative = "background-Hepatocyte")
    getCellGroupDEG(adata_all.copy(), positive= "pdx-Tumor cell", negative = "tumor-Tumor cell")