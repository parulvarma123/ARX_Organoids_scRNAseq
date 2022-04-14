library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

#Seurat V4 was used for differential gene expression analysis
#Used EnchancedVolcano to make Volcano plots

#Load the previously saved hCS30DIVCelltype.rds file
hCS30D <- readRDS(file = "~/hCS30DIVCelltype.rds")
DimPlot(hCS30D)

# Differential gene expression analysis between Control and ARX Mutant Radial Glial Cells
hCS30D$celltype.OrgIdent <- paste(Idents(hCS30D), hCS30D$OrgIdent, sep = "_")
hCS30D$celltype <- Idents(hCS30D)
Idents(hCS30D) <- "celltype.OrgIdent"
RadialGlia <- FindMarkers(hCS30D, ident.1 = "Radial Glia_Mut", ident.2 = "Radial Glia_Ctrl", verbose = FALSE)
head(RadialGlia, n = 15)
write.table(RadialGlia, file ="~/hCS30D_RadialGliaDGE.txt", sep = "\t", quote=F,row.names =T)

#save the file and one can rearrange and clean in excel and save as .csv
#Load the .csv file and make volcano plot

res1 <- read.csv(file ='~/hCS30D_RadialGliaDGE.csv', header = TRUE)

#Make rownames 
row.names(res1) <- res1$Gene

#Make different colors for logFC>1 and LogFC<1.0
keyvals <- ifelse(
  res1$avg_log2FC < -1.0, 'cornflowerblue',
  ifelse(res1$avg_log2FC > 1.0, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'p<0.05 & logFC>1'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'cornflowerblue'] <- 'p<0.05 & logFC<1'

#Volcanoplot
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "hCS30D_RadialGlia_DGE",
                pCutoff = 1.30103,
                FCcutoff = 1.0,
                pointSize = 1.5,
                labSize = 2.0,
                colCustom = keyvals,
                colAlpha = 1,
                cutoffLineType = "blank",
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                legendLabels = c("Not Sig.", "logFC", "Sig p-value", "p<0.05 & abs(logFC)>1"),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 3.0)
ggsave("~/hCS30D_RadialGlia_VolcanoPlot.jpeg", width = 9, height = 7, units = c("in"), dpi = 300, bg= "white")
dev.off()
