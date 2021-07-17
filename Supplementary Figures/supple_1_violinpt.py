import pandas as pd 
import numpy as np
import scanpy as sc
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

font = {'family' : 'Helvetica',
        'size'   : 22}
matplotlib.rc('font', **font)

adata = sc.read("hb_integrated_normalized.h5ad")
df_cells = pd.read_csv("ctype_integrated_harmony_TumorClustersAdded_cluster11Added.txt",sep="\t",header=0,index_col=0,keep_default_na=False)

adata = adata[list(df_cells.index.values),:]
adata.obs = df_cells

adata.obs["sample_group_v2"] = ["bkg" if i == "background" else i for i in adata.obs["sample_group"]]
adata.obs["Cell.class_v2"] = [(adata.obs["Cell.class"][i] + "-" + adata.obs["sample_group_v2"][i]) for i in range(adata.obs.shape[0])]

markers = ["MZB1","MS4A1","CD79A","CD79B",\
           "SCGB3A1","PROM1","CFTR","EPCAM","KRT7","KRT19","FXYD2","DEFB1","CD24","SPP1","SOX9","KRT18",\
           "PECAM1","VWF","CLEC4G","CLEC4M","FLT1",\
           "CRYAB","PRNP","CYGB","TAGLN","ACTA2","COL3A1","COL1A1","COL1A2","SPARC","SERPINE1",\
           "TF","CYP3A4","CYP2E1","APOB","ASGR1","PCK1","HP","ASS1","APOE","KRT8",\
           "LYZ","VCAN","IL18","CD68","CD163","S100A9","MAFB","VSIG4","CD5L","MARCO","CPVL","VCAM1",\
            "KLRB1","CD8A","CD3E",\
           "GPC3","DLK1","DUSP9","YAP1","CTNNB1","NKD1","COL2A1"]

sc.pl.stacked_violin(adata, markers, groupby='Cell.class_v2', rotation=90,figsize=(30,10), row_palette = sns.color_palette("pastel"), save = "Supp1_violin.png")
