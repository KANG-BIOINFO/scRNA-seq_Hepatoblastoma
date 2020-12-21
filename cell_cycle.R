library(Seurat)
library(plyr)
library(ggplot2)
library(SeuratWrappers)
library(RColorBrewer)
library(tibble)
library(ggpubr)

setwd("/data/aronow/Kang/SingleCellHBPilot-HB17/results/cellxgene/")
brks <- c(0,0.25,0.5,0.75,1)


### loading datasets 
trimmed_hepatoblastoma = readRDS("trimmed_hepatoblastoma_v2.RDS")

# load cell cycle genes
cellCycleGenes = read.table("regev_lab_cell_cycle_genes.txt",sep="\t",header=1)
s_genes = as.character(cellCycleGenes$s_genes[1:43])
s_genes = intersect(s_genes, as.character(rownames(trimmed_hepatoblastoma)))
g2m_genes = intersect(as.character(cellCycleGenes$g2m_genes), as.character(rownames(trimmed_hepatoblastoma)))
all.genes = rownames(trimmed_hepatoblastoma)
trimmed_hepatoblastoma = ScaleData(trimmed_hepatoblastoma, features = all.genes)

# calculate s scores and g2m scores by the sum of s genes or g2m genes
subset_hepatoblastoma = trimmed_hepatoblastoma[s_genes,]
trimmed_hepatoblastoma$s_score = apply(subset_hepatoblastoma@assays$RNA@scale.data, 2, sum)

subset_hepatoblastoma = trimmed_hepatoblastoma[g2m_genes,]
trimmed_hepatoblastoma$g2m_score = apply(subset_hepatoblastoma@assays$RNA@scale.data, 2, sum)


# draw box plot for 3 cell types 
c = trimmed_hepatoblastoma[,trimmed_hepatoblastoma$cell_type_v2 %in% c("Liver Background_Hepatocytes", 
                                                                      "Human Hepatoblastoma_Tumor cells",
                                                                      "PDX_Tumor cells")]
dim(c)
mycomparison = list(c("Liver Background_Hepatocytes","Human Hepatoblastoma_Tumor cells"), 
                    c("Liver Background_Hepatocytes","PDX_Tumor cells"),
                    c("Human Hepatoblastoma_Tumor cells", "PDX_Tumor cells"))
c$cell_type_v2 = as.character(c$cell_type_v2)
c$cell_type_v2 = factor(c$cell_type_v2, levels = c("Liver Background_Hepatocytes", 
                                                   "Human Hepatoblastoma_Tumor cells",
                                                   "PDX_Tumor cells"))
ggplot(c@meta.data, aes(x=cell_type_v2, y=g2m_score, fill=cell_type_v2)) + 
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


# box plot for tumor clusters
# s scores
d = trimmed_hepatoblastoma[,trimmed_hepatoblastoma$tumor_cluster_v2!=""]

ggplot(d@meta.data, aes(x=tumor_cluster_v2, y=s_score, fill=tumor_cluster_v2)) + 
        geom_boxplot(width = 0.5, lwd = 0.8, outlier.shape = NA) + 
        scale_fill_manual(values = c("#FFFFFF","#F8766D","#D89000","#A3A500","#02BF7D","#00BFC4","#02B0F6","#9490FF","#E76BF3","#FF62BC"))+
        xlab("") + ylab("") + ylim(c(-15,60)) +
      theme(panel.border = element_rect(fill = NA, colour = "black", size = 3),
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey", size = 0.25),
            axis.text.y = element_text(family = "Helvetica", size = 26),
            axis.text.x = element_text(family = "Helvetica", size = 16, face = "bold"),
            legend.position = "none") 

# g2m scores
ggplot(d@meta.data, aes(x=tumor_cluster_v2, y=g2m_score, fill=tumor_cluster_v2)) + 
        geom_boxplot(width = 0.5, lwd = 0.8, outlier.shape = NA) + 
        scale_fill_manual(values = c("#FFFFFF","#F8766D","#D89000","#A3A500","#02BF7D","#00BFC4","#02B0F6","#9490FF","#E76BF3","#FF62BC"))+
        xlab("") + ylab("") + ylim(c(-25,100))+
        theme(panel.border = element_rect(fill = NA, colour = "black", size = 3),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "grey", size = 0.25),
              axis.text.y = element_text(family = "Helvetica", size = 26),
              axis.text.x = element_text(family = "Helvetica", size = 16, face = "bold"),
              legend.position = "none") 













