library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(plyr)
library(EnhancedVolcano)
library(ggpubr)

hb.3groups = readRDS("HB_integrate_reciPCA_7samples_harmony_perMito10_3types.RDS")

# scatter plot of hepatocytes and tumor cells in 3 sample groups (Fig.4a)
Compare_genes = c("GPC3", "YAP1", "EPCAM", "DLK1", "KRT19","APOC3", "ALB", "HP","DUSP", "CTNNB1", "DLK1", "NKD1")
df_superbin = read.table("./superbinexpr_sampleGroup_CellClass_0302_concise.txt",
                         sep = "\t", header = T, row.names = 1) # superbinexor files are average expression table for each cell class in each sample group. They can downloaded from ToppCell.
df_superbin = df_superbin[,c("background.Hepatocyte.6648","tumor.Tumor.cell.34307","pdx.Tumor.cell.18124")]

ggplot(df_superbin, aes(x = `background.Hepatocyte.6648`,
                        y = `tumor.Tumor.cell.34307`)) + geom_point() +
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = "black", size = 2, fill = NA),
              axis.text = element_text(size = 20),
              axis.ticks.length = unit(0.25,"cm")) + xlab("") + ylab("") +
        geom_text_repel(data = df_superbin[Compare_genes,],
                        family = "Helvetica",
                        aes(label = rownames(df_superbin[Compare_genes,]), fontface = "bold"),
                        size = 5, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.35,"lines"),
                        colour = "#CC0000", segment.color = "grey", )

hb.lm = lm(`pdx.Tumor.cell.18124` ~ `tumor.Tumor.cell.34307`, data = df_superbin)
summary(hb.lm)$r.squared


# degs and volcano plot (Fig.4b)
hb.integrated$Cell.class_v2 = paste0(hb.integrated$sample_group, "-", hb.integrated$Cell.class)
# deg-A
deg_tumor_bg = read.table("./DEGs/Fig4/tumor-Tumor cell vs. background-Hepatocyte.txt",
                          sep = "\t", header=1)
deg_tumor_bg$names = as.character(deg_tumor_bg$names)
deg_tumor_bg = deg_tumor_bg[(deg_tumor_bg$pts > 0.05) | (deg_tumor_bg$pts_rest > 0.05),] # remove lowly expressed genes
deg_tumor_bg[deg_tumor_bg$logfoldchanges > 10,"logfoldchanges"] = 10
deg_tumor_bg[deg_tumor_bg$logfoldchanges < -10,"logfoldchanges"] = -10
EnhancedVolcano(deg_tumor_bg, lab = deg_tumor_bg$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("GPC3","DLK1","AFP", "DCDC2", "RPS27A", "RNF43","HMGA2", "NKD1", "NOTUM",
                              "HP","C3","CYP2B6","CRP","GHR","ACSL1","NAMPT","HPX","C9","NMT","PCK1","SDS"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "Hepatocyte (Bg) vs. Tumor Cells (Tumor)", subtitle = NULL, caption = NULL,transcriptLabSize=4.0,
                xlim = c(-10,10))

# deg-B
deg_pdx_bg = read.table("./DEGs/Fig4/pdx-Tumor cell vs. background-Hepatocyte.txt",
                        sep = "\t", header=1)
deg_pdx_bg$names = as.character(deg_pdx_bg$names)
deg_pdx_bg = deg_pdx_bg[(deg_pdx_bg$pts > 0.05) | (deg_pdx_bg$pts_rest > 0.05),] # remove lowly expressed genes
deg_pdx_bg[deg_pdx_bg$logfoldchanges > 10,"logfoldchanges"] = 10
deg_pdx_bg[deg_pdx_bg$logfoldchanges < -10,"logfoldchanges"] = -10
EnhancedVolcano(deg_pdx_bg, lab = deg_pdx_bg$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("WNK2","CTNND2","GPC3","DLK1","AFP", "DCDC2", "RPS27A", "RNF43","HMGA2", "NKD1", "NOTUM",
                              "HP","C3","CYP2B6","CRP","GHR","ACSL1","NAMPT","HPX","C9","NMT","PCK1","SDS"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "Hepatocyte (Bg) vs. Tumor Cells (PDX)", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-10,10))

# deg-C
deg_pdx_tumor = read.table("./DEGs/Fig4/pdx-Tumor cell vs. tumor-Tumor cell.txt",
                           sep = "\t", header=1)
deg_pdx_tumor$names = as.character(deg_pdx_tumor$names)
deg_pdx_tumor = deg_pdx_tumor[(deg_pdx_tumor$pts > 0.05) | (deg_pdx_tumor$pts_rest > 0.05),] # remove lowly expressed genes
deg_pdx_tumor[deg_pdx_tumor$logfoldchanges > 5,"logfoldchanges"] = 5
deg_pdx_tumor[deg_pdx_tumor$logfoldchanges < -5,"logfoldchanges"] = -5
EnhancedVolcano(deg_pdx_tumor, lab = deg_pdx_tumor$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("KCNQ1","APBB2","CASC9","CYP2E1","IGFBP1","CFH"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "Tumor Cells (Tumor) vs. Tumor Cells (PDX)", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-5,5))


# tumor genes (Fig.4c)
DefaultAssay(hb.3groups) = "RNA"
VlnPlot(hb.3groups,features = c("DLG1"), group.by = "Cell.class_v2",pt.size = 0,y.max = 5,cols = c("#4572C4","#FFCC99","#FF6666")) +
  coord_flip() + NoLegend() + theme(panel.border = element_rect(fill=NA,colour="black", size = 2),
                                    plot.title = element_blank(), axis.text.x = element_blank(),
                                    axis.text = element_blank(),axis.title = element_blank(),
                                    axis.ticks = element_blank())


# cell cycle scores (Fig.4d)
# load cell cycle genes
cellCycleGenes = read.table("/data/aronow/Kang/data/regev_lab_cell_cycle_genes.txt",sep="\t",header=1)
s_genes = as.character(cellCycleGenes$s_genes[1:43])
s_genes = intersect(s_genes, as.character(rownames(hb.integrated)))
g2m_genes = intersect(as.character(cellCycleGenes$g2m_genes), as.character(rownames(hb.integrated)))
all.genes = rownames(hb.integrated)

DefaultAssay(hb.integrated) = "RNA"
hb.integrated = ScaleData(hb.integrated, features = all.genes)

# calculate s scores and g2m scores by the sum of s genes or g2m genes
subset_hb = hb.integrated[s_genes,]
hb.integrated$s_score = apply(subset_hb@assays$RNA@scale.data, 2, sum)
subset_hb = hb.integrated[g2m_genes,]
hb.integrated$g2m_score = apply(subset_hb@assays$RNA@scale.data, 2, sum)

# draw box plot for 3 cell types
hb.integrated$Cell.class_v2 = paste0(hb.integrated$sample_group, "_", hb.integrated$Cell.class)
c = hb.integrated[,hb.integrated$Cell.class_v2 %in% c("background_Hepatocyte", "tumor_Tumor cell","pdx_Tumor cell")]
mycomparison = list(c("background_Hepatocyte","tumor_Tumor cell"),
                    c("background_Hepatocyte","pdx_Tumor cell"),
                    c("tumor_Tumor cell", "pdx_Tumor cell"))
c$Cell.class_v2 = as.character(c$Cell.class_v2)
c$Cell.class_v2 = factor(c$Cell.class_v2, levels = c("background_Hepatocyte",
                                                     "tumor_Tumor cell",
                                                     "pdx_Tumor cell"))

ggplot(c@meta.data, aes(x=`Cell.class_v2`, y=s_score, fill=`Cell.class_v2`)) +
  geom_boxplot(width = 0.5, lwd = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D","#02BA38","#619CFF"))+
  xlab("") + ylab("") + ylim(c(-20,75))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 3),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", size = 0.25),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "Helvetica", size = 26),
        legend.position = "none") +
  stat_compare_means(comparisons = mycomparison, method = "t.test", label.y = c(50,60,70), size = 6)
