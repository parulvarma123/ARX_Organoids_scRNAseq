library(Seurat)
library(dplyr)
library(ggplot2)

#Heatmaps for all datasets

#Heatmap for hCS30DIV dataset
#Open hCS30DIV Cell Type File
hCS30D <- readRDS(file = "~/hCS30DIVCelltype.rds")
DimPlot(hCS30D)

#Subset all the controls for making heatmap
hCS30Dctrl1 <- subset(hCS30D, subset = OrgIdent == "Ctrl")

#Subsampling 2000 cells from control subsets for heatmap 
subsampledhCS30DIV <- hCS30Dctrl1[, sample(colnames(hCS30Dctrl1), size =2000, replace=F)]

#Preferred order for cell types
y=levels(subsampledhCS30DIV)
mylevels1 = levels(subsampledhCS30DIV)[c(1,3,7,5,6,2,8,4)]
Idents(subsampledhCS30DIV) <- factor(Idents(subsampledhCS30DIV), levels = mylevels1)

#Genes used in the heatmap
features <- c("VIM", "HES1", "HES5", "PAX6", "NES", "EMX1", "SOX2", "MKI67", "TOP2A", "DLL3", "DLL1", "NHLH2", "LHX9", "TBR1", "DLX5", "GAD2", "S100B", "SOX10")

#Making heatmap
heatmaphCS30D <- DoHeatmap(subsampledhCS30DIV,
                           features = features,
                           assay = 'RNA',
                           group.by = "ident", 
                           slot = "data", 
                           lines.width = 3,
                           disp.min = 1.0,
                           disp.max = 2.0,
                           group.colors = c("hotpink1", "deepskyblue","violetred4","gold1", "orangered","olivedrab4","purple3", "gray"),
                           group.bar.height = 0.05)
heatmaphCS30D + scale_fill_gradientn(colors = c("gray96",  "red"))
ggsave("~/Heatmaps/hCS30DHeatmap.tiff", width = 18, height = 8, units = c("in"), dpi = 300, bg= "white")
dev.off()

##############################
#Heatmap for hCS120DIV dataset

#Open hCS120DIV Cell Type File 
hCS120DIV <- readRDS(file = "~/hCS120DIVCelltype.rds")
DimPlot(hCS120DIV)

#Subset all the controls for making heatmap
hCS120DCtrl <- subset(hCS120DIV, subset = OrgIdent == "Ctrl")

#Subsample 2000 cells from control subset for making heatmap
subsampled120DIV <- hCS120DCtrl[, sample(colnames(hCS120DCtrl), size =2000, replace=F)]

#Preferred order of cell types
y2=levels(subsampled120DIV)
mylevels3 =  levels(subsampled120DIV)[c(8,2,6,7,1,3,5,9,4)]
Idents(subsampled120DIV) <- factor(Idents(subsampled120DIV), levels = mylevels3)

#Genes used in the heatmap
genes = c("VIM", "HES1", "HES5","SOX2","PAX6", "NES","EMX1", "TNC","FAM107A", "HOPX", "MKI67", "TOP2A", "EOMES","NEUROD1", "NEUROD6", "BCL11B", "TBR1", "SATB2", "CUX2", "RELN", "SCGN", "DLX5", "GAD2", "PAX2")

#Making heatmap
heatmaphCS120DIV <- DoHeatmap(subsampled120DIV,
                                features = genes,
                                assay = 'RNA',
                                group.by = "ident",
                                slot = "data", 
                                lines.width = 3,
                                disp.min = 0.1,
                                disp.max = 2.0,
                                group.colors = c("hotpink1", "navyblue", "deepskyblue", "violetred4","gold1", "orangered", "chartreuse", "green4","gray"),
                                group.bar.height = 0.05)

heatmaphCS120DIV + scale_fill_gradientn(colors = c("gray96",  "red"))
ggsave("~/Heatmaps/hCS120DIVHeatmap.tiff", width = 18, height = 8, units = c("in"), dpi = 300, bg= "white")
dev.off()

#############################
#Heatmap for hSS30DIV dataset

#Open hSS30DIV Cell Type File
hSS30DIV <- readRDS(file = "~/hSS30DIVCelltype.rds")

#Subset all the controls for making heatmap
hSS30DCtrl <- subset(hSS30DIV, subset = OrgIdent == "Ctrl")

#Subsample 2000 cells from control subset for making heatmap
subsampledhSS30DIV <- hSS30DCtrl[, sample(colnames(hSS30DCtrl), size =2000, replace=F)]

#Preferred order of cell types
z= levels(subsampledhSS30DIV)
mylevels2 =  levels(subsampledhSS30DIV)[c(2,3,1,4,7,5,6)]
Idents(subsampledhSS30DIV) <- factor(Idents(subsampledhSS30DIV), levels = mylevels2)

#Genes used in the heatmap
genes1 = c("VIM", "HES1", "HES5", "MKI67", "TOP2A","OTX2", "NKX2-1", "DLX1","DLX2", "DLX5", "LHX8", "GAD1", "GAD2",  "CALB1", "CALB2", "CXCR4","MAFB", "SST", "S100B", "SOX10")

#Making heatmap
heatmaphSS30DIV <- DoHeatmap(subsampledhSS30DIV,
                           features = genes1,
                           assay = 'RNA',
                           group.by = "ident", 
                           slot = "data", 
                           lines.width = 3,
                           disp.min = 0.1,
                           disp.max = 2.0,
                           group.colors = c("hotpink1", "deepskyblue", "violetred4","olivedrab3", "purple3","gray28", "gray"),
                           group.bar.height = 0.05)
heatmaphSS30DIV + scale_fill_gradientn(colors = c("gray96", "red"))
ggsave("~/Heatmaps/hSS30DIVHeatmap.tiff", width = 18, height = 8, units = c("in"), dpi = 300, bg= "white")
dev.off()

#############################
#Heatmap for hSS120DIV dataset

#Open hSS120DIV Cell Type File
hSS120DIV <- readRDS(file = "~/hSS120DCelltype.rds")
DimPlot(hSS120DIV)

#Subset all the controls for making heatmap
hSS120DCtrl <- subset(hSS120DIV, subset = OrgIdent == "Ctrl")

#Subsample 2000 cells from control subset for making heatmap
subsampledhSS120DIV <- hSS120DCtrl[, sample(colnames(hSS120DCtrl), size =2000, replace=F)]

#Preferred order of cell types
a= levels(subsampledhSS120DIV)
mylevels4 =  levels(subsampledhSS120DIV)[c(3,7,2,5,6,4,1)]
Idents(subsampledhSS120DIV) <- factor(Idents(subsampledhSS120DIV), levels = mylevels4)

#Genes used in the heatmap
genes4 = c("VIM", "HES1", "HES5", "MKI67", "TOP2A","OTX2", "NKX2-1", "DLX1","DLX2", "DLX5", "LHX8", "GAD1", "GAD2","CALB1", "CALB2", "CXCR4","MAFB", "SST")

#Making heatmap
heatmaphSS120DIV <- DoHeatmap(subsampledhSS120DIV,
                            features = genes4,
                            assay = 'RNA',
                            group.by = "ident", 
                            slot = "data", 
                            lines.width = 3,
                            disp.min = 0.1,
                            disp.max = 2.0,
                            group.colors = c("hotpink1", "deepskyblue", "violetred4","olivedrab3", "darkgreen","gray28", "gray"),
                            group.bar.height = 0.05)

heatmaphSS120DIV + scale_fill_gradientn(colors = c("gray90",  "red"))
ggsave("~/Heatmaps/hSS120DIVHeatmap.tiff", width = 18, height = 8, units = c("in"), dpi = 300, bg= "white")
dev.off()
