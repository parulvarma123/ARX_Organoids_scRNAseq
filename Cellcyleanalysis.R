library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggridges)
library(dplyr)
library(harmony)
library(ggplot2)

hCS30D <- readRDS(file = "~/hCS30DIVCelltype10_25_2023.rds")
hCS30D <- UpdateSeuratObject(hCS30D)
DimPlot(hCS30D, reduction = "umap")

#Subset for Cycling Progenitors
CP <- subset(hCS30D, idents = c("Cycling Progenitors"))
CP <- NormalizeData(CP)
CP <- FindVariableFeatures(CP, selection.method = "vst")
CP <- ScaleData(CP, features = rownames(CP))
CP <- RunPCA(CP, features = VariableFeatures(CP))

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Assign cell cycle scores 
#First, we assign each cell a score, based on its expression of G2/M and S phase markers. 
#These marker sets should be anticorrelated in their expression levels, 
#and cells expressing neither are likely not cycling and in G1 phase.

CP <- CellCycleScoring(CP, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(CP[[]])
aa <- table(CP$Phase, CP$orig.ident)
write.table(aa, file ="~/Cell_cylce_phases_hCS30DIV_CyclingProg_030625.txt", sep = "\t", quote=F,row.names =T)

bb <- table(CP$Phase, CP$OrgIdent)
write.table(bb, file ="~/Cell_cylce_phases_hCS30DIV_CyclingProg_Ctrl_vs_Mut_030625.txt", sep = "\t", quote=F,row.names =T)

cc <- table(CP$orig.ident)