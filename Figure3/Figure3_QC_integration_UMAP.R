library(Seurat)
library(ggplot2)
library(plyr)
library(harmony)

setwd("SingleCellHBPilot-HB17/integrate/")

# 1. read raw data
count = 1
for (dir in list.dirs("../processed/"))
{
  if (!dir %in% c("../processed/","../processed//hb30-background")) # remove hb30-background due to low quality
  {
    sample_name = strsplit(dir, "//")[[1]][2]
    hb.data = Read10X_h5(paste0(dir, "/filtered_feature_bc_matrix.h5"))

    hb = CreateSeuratObject(count = hb.data, min.cells = 0, min.features = 0)
    hb[["percent.mt"]] = PercentageFeatureSet(hb, pattern = "^MT-")

    # apply QC to each sample
    hb = hb[, (hb[["nCount_RNA"]]>800) & (hb[["nFeature_RNA"]] > 500) & (hb[["percent.mt"]] < 10)]
    hb = RenameCells(hb, new.names = paste0(sample_name, "_", colnames(hb)))
    hb$sample = sample_name

    if (count == 1) {hb.combined = hb}
    else {hb.combined = merge(hb.combined, y = hb)}

    count = count + 1
  }
}

# harmony integration
hb.combined = NormalizeData(hb.combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
hb.combined = RunHarmony(hb.combined, group.by.vars = "sample") # remove batch from each sample
hb.combined = RunUMAP(hb.combined, reduction = "harmony", dims = 1:30)
hb.combined = FindNeighbors(hb.combined, reduction = "harmony", dims = 1:30)
hb.combined <- FindClusters(hb.combined, resolution = 0.5) # use 3 resolutions
hb.combined <- FindClusters(hb.combined, resolution = 1)
hb.combined <- FindClusters(hb.combined, resolution = 2)

saveRDS(hb.combined, file = "./HB_integrate_reciPCA_7samples_harmony_perMito10.RDS")

# umap of cell classes (fig3a)
DimPlot(hb.combined, group.by = "Cell.class_concise", label = T, label.size = 8) +
  theme(legend.position = "none",
        legend.text = element_text(family = "Helvetica", size = 30),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

# umap of sample groups (fig3a)
DimPlot(hb.combined, group.by = "sample_group", cols = c("#3399FF","#FFCC99","#FF6666")) +
  theme(legend.position = "none",
        legend.text = element_text(family = "Helvetica", size = 30),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

# feature plots (figure 3b)
DefaultAssay(hb.combined) = "RNA"
FeaturePlot(hb.combined, features = c("GPC3","CYP3A4","OL6A3",
                                      "CD163","FLT1","PTPRC"), cols = c("#F5F5F5","red"), ncol = 3)
