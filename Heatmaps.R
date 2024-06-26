library(Seurat)
library(dplyr)
library(ggplot2)

#Heatmaps for all datasets

###########################
#Heatmap for hCO30D,
#Have cleaned this by moving genes up and down 
#And also removing some genes which were not adding additional value
###########################

hCS30D <- readRDS(file = "/work/nac970/ARXWork2023/Scripts_ARX_Revision_2023/hCS30DIVCelltype10_25_2023.rds")
DimPlot(hCS30D)
hCS30Dctrl1 <- subset(hCS30D, subset = OrgIdent == "Ctrl")

#If we want to plot equal number of cells in each group
#maxcells  <- min(table(Idents(hCS30Dctrl1)))

subsampled <- hCS30Dctrl1[, sample(colnames(hCS30Dctrl1), size =2000, replace=F)]
y=levels(subsampled)

features <- c("FABP7", "CD99", "SFRP1", "HES5", "TTYH1", "CREB5", "C1orf61","DLK1", "PTPRZ1", "LINC01551", 
              "LIX1", "LHX2", "HMGA2", "B3GAT2", "IFI44L", "FOXG1", "SOX9", "EMX2", "FZD8", "HES1",
              "GLI3", "SLC1A3","PCLAF", "HMGB2", "TYMS","MAD2L1", "NFIA","PTTG1", "CENPF","CENPW","TACC3", "CDCA8", 
              "CCNA2", "PRC1", "HIST2H2AC", "NFIB",  "DMRTA2", "EMX1", "HIST1H1B", "HIST1H1D", "HIST1H1C", 
              "HIST1H2AJ", "HIST1H1E", "TOP2A", "CDK1", "NUSAP1", "HIST1H3B", "RRM2", "SPC25", "PBK", 
              "MKI67", "HIST1H3G", "NDC80", "UBE2C", "GTSE1", "AURKB", "HIST1H2AG", "RAD51AP1", "CDC20", 
              "PLK1", "CCNB1", "DLGAP5", "CCNB2", "CDKN3", "TPX2", "BIRC5", "SGO2", "ASPM", "HMMR", "TROAP",
              "CDCA3", "PSRC1",  "RSPO3", "RSPO2", "RSPO1", "TPBG", "WLS", "ZNF503", "PLS3",
              "TRH",  "SLIT2",  "HIST1H4D",  "AL031777.3",  "DLL3", "GADD45G", "RGS16", "DLL1", "ASCL1",
              "LMO1", "INSM1", "SMOC1", "SRRM4", "DLX2", "DLX1", "NEUROD4", "RASD1",
              "CBFA2T2", "MFNG", "NEUROG1", "KLHDC8A", "TBR1","FEZF2", "SAMD3", "PCSK2", "NRN1", "NSG2", "LHX5-AS1",
              "IGSF21","SLC17A6", "NRXN1", "TMEM163", "RIPOR2", "MYT1L", "CACNA2D1", "CHST8", "LHX5", "LHX2",
              "FOXP2", "PLPPR1", "CCSAP", "BCL11B","NEFM", "NEFL", "LHX9", "C1QL4", "ONECUT2", "EBF3", "SYT4", "EBF1",
              "MAB21L1", "PRPH", "SNCG", "ADCYAP1",  "NHLH2", "CACNA2D1","SLC17A6", "GREM2", "ONECUT1", "SSTR2")

#x=levels(hCS30Dctrl1)

mylevels =  levels(subsampled)[c(1,3,7,5,6,2,4)]
Idents(subsampled)<- factor(Idents(subsampled), levels = mylevels)
#Idents(hCS30Dctrl1) <- factor(Idents(hCS30Dctrl1), levels = mylevels)

heatmaphCS30D <- DoHeatmap(subsampled,
                           features = features,
                           group.by = "ident", 
                           assay = "RNA",
                           slot = "data", 
                           lines.width = 5,
                           disp.min = 1.0,
                           disp.max = 2.0,
                           group.colors = c("hotpink1", "deepskyblue","violetred4","gold1", "orangered", "olivedrab3", "lightgray"),
                           group.bar.height = 0.05)

heatmaphCS30D + scale_fill_gradientn(colors = c("gray98",  "red"))

ggsave("/work/nac970/ARXWork2023/hCS30D_Reanalysis_10_25_2023/hCS30DHeatmap_Cleaned_02_22_2024.tiff", width = 25, height = 40, units = c("in"), dpi = 300, bg= "white")
dev.off()

###########################
#Heatmap for hCO120DIV
#Cleaned
###########################
#Open hCS120DIV Cell Type File 

hCS120D <- readRDS(file = "/work/nac970/ARXWork2023/Scripts_ARX_Revision_2023/hCS120DIVCelltype10_25_2023.rds")
DimPlot(hCS120D)

#Subset all the controls for making heatmap
hCS120DCtrl <- subset(hCS120D, subset = OrgIdent == "Ctrl")

#Subsample 2000 cells from control subset for making heatmap
subsampled120DIV <- hCS120DCtrl[, sample(colnames(hCS120DCtrl), size =2000, replace=F)]

#Preferred order of cell types
y2=levels(subsampled120DIV)
mylevels3 =  levels(subsampled120DIV)[c(8,2,5,7,1,6,4,3)]
Idents(subsampled120DIV) <- factor(Idents(subsampled120DIV), levels = mylevels3)

#Genes used in the heatmap
genes = c("HES5", "HES6","SFRP1", "ASCL1", "DLL1", "SMOC1", "BCAN", "ETV1", "SCRG1", "AC004540.2",
          "RFX4", "LIMA1","MEIS2", "DOK5", "SOX2","TTYH1", "ID4","PON2",
          "NES",  "HOPX", "FAM107A", "CRYAB", "EDNRB", "S100B", "MGST1", "PMP2",
          "SLC1A3", "MT2A", "SPARCL1", "NTRK2", "TFPI", "MFGE8", "SERPINE2", "FAM89A", "S100A16", 
          "PLA2G16",  "IRX1", "RSPO3", "TPBG", "EFEMP1", "LINC01088", "CYBA", "ANXA2", "RSPO2",
          "IFITM3", "WLS", "HTRA1", "C5orf49", "S100A6", "NUPR1", "PLTP", "PIFO", "ATP1A2", 
          "LHX2",  "EMX2", "GINS2",  "IQGAP2", "GLI3", "SOX3", 
          "SYNE2", "TCIM", "COL11A1", "FGFBP3","LINC01896", "AC092958.1",  "PCLAF", "NUSAP1", "UBE2C", "BIRC5", "CENPF",
          "PBK", "CDK1", "MKI67", "PCLAF", "SPC25", "TPX2", "TOP2A", "NUF2", "CCNB2", "ASPM",
          "GTSE1","MAD2L1","ZWINT","PTTG1","SGO1","CDKN3","TAC3","NHLH1","EOMES","HES6","NEUROG2",
          "EMX1", "RCOR2", "INSM1", "IGFBPL1", "CRABP1", "NRN1", "MGP", "LHX2", "SFRP1", "GADD45G", "AC007952.4", "LINC01551",
          "ELAVL4", "GOLGA8A", "GPC2","MEF2C", "TBR1","FEZF2", "NEUROD6", "SLA", "NEUROD2", "GPR22", "PDE1A", "AC004158.1", 
           "FOXG1", "BCL11A", "LMO3", "SOX5",  "GRIA2", "SNCA", "LINC01551",
          "CAMKV", "ZBTB18", "CRYM","CELF4", "VSNL1", "GABRG2", "RELN", "RGS8", "UNC5D", "POU3F2", "NEFL", "CALY",
          "SCG2", "CALB2", "SCG5", "SCN2A", "SYT5", "CHD5", "CAMK2N2", "RAB3C", "ZCCHC12", "FXYD7",
          "CHGB", "ATP1A3", "SCN3B", "MAB21L1", "TAC1", "DTX4")

#Making heatmap
heatmaphCS120DIV <- DoHeatmap(subsampled120DIV,
                              features = genes,
                              assay = 'RNA',
                              group.by = "ident",
                              slot = "data", 
                              lines.width = 3,
                              disp.min = 0.1,
                              disp.max = 2.0,
                              group.colors = c("hotpink1", "navyblue", "deepskyblue", "violetred4","gold1", "orangered", "olivedrab3", "lightgray"),
                              group.bar.height = 0.05)

heatmaphCS120DIV + scale_fill_gradientn(colors = c("gray98",  "red"))
ggsave("/work/nac970/ARXWork2023/hCS120D_Reanalysis_10_25_2023/hCS120DIVHeatmap_Cleaned_All_Genes_02_22_2024.tiff", width = 25, height = 40, units = c("in"), dpi = 300, bg= "white")
dev.off()

#############################
#Heatmap for hGEO30DIV dataset
#Cleaned up 
####################################

#Open hSS30DIV Cell Type File
hSS30D <- readRDS(file = "/work/nac970/ARXWork2023/Scripts_ARX_Revision_2023/hSS30D_Celltype_11_02_2023.rds")
DimPlot(hSS30D)
#Subset all the controls for making heatmap
hSS30DCtrl <- subset(hSS30D, subset = OrgIdent == "Ctrl")

#Subsample 2000 cells from control subset for making heatmap
subsampledhSS30DIV <- hSS30DCtrl[, sample(colnames(hSS30DCtrl), size =2000, replace=F)]

#Preferred order of cell types
z= levels(subsampledhSS30DIV)
mylevels2 =  levels(subsampledhSS30DIV)[c(2,3,1,4,6,7,5)]
Idents(subsampledhSS30DIV) <- factor(Idents(subsampledhSS30DIV), levels = mylevels2)

#Genes used in the heatmap
genes1 = c("PTPRZ1","TTYH1", "C1orf61", "LINC01551","VEPH1","FABP7", "HES5","HES1", "SLC6A8", "MGST1", "AL139246.5", "GPC3", "FRZB", "ADGRV1", "SFTA3", "DLK1", 
           "RAB13","DLL1","RGS16", "GADD45G","DLL3","ASCL1","INSM1", "LMO1","GINS2", 
            "SNHG18", "OLIG1",  "LINC01833", "CD82",
           "FZD5", "NKX6-2", "NTRK2","NKX2-8","SLC1A3","DIO3", "NKX2-2", "NKX2-2", "CTGF",
           "IRX2", "AC099489.1", "C5orf38", "RFX4", "VCAM1", "PLP1", "LITAF", "PDLIM3", "PTX3",
           "FGFBP3", "ENPP2","PDE1A","LRP2","FGFR3","NR2F2", "TOX3", "NKAIN4", "NKX2-1", 
           "CCNB1IP1","RAB34",
           "QPRT", "MDFI", "EVA1B", "NUSAP1", "MKI67", "HIST1H1B", "PBK",
           "GTSE1","MAD2L1", "CDKN3","DLGAP5","TOP2A","CDK1","PCLAF","SPC25","HIST1H1D","SGO1",
           "UBE2C","PIMREG","CDCA3","CENPF","CENPN","C21orf58","BIRC5","NKX2-8","NKX6-2","NKX2-2",
           "CENPU",
           "CENPA","PDLIM3","NTRK2","PRC1","NEK2","ST18","CKAP2L","MT1X","INSM1","CDCA8","RASD1",
           "ROBO3","RIPOR2", "NTRK1" ,"MIR7-3HG", "LHX8",
           "DLX1","DLX2", "DLX6-AS1","ARX","DLX5","ARL4D","ACTL6B","FGF14","SEZ6L2","AKR1C1","DLX1","NXK2-1", "SRRM4",
           "RTN1","INA","TTC9B","NSG2","GNG3", "GAP43","PPP2R2B","VGF","CELF4","SCG5","DNER","RALYL","TMOD2","SCG2",
           "GRIA2","TAC1","ISL1","CALY","GRID2","SCN2A","TTC9B","AKR1C2",
           "CAMK2B","RASD1","NKX2-1",
           "DLX5", "DLX6","SST","GAD1","DLX6-AS1","CHD5","MYT1L")

#Making heatmap
heatmaphSS30DIV <- DoHeatmap(subsampledhSS30DIV,
                             features = genes1,
                             assay = 'RNA',
                             group.by = "ident", 
                             slot = "data", 
                             lines.width = 3,
                             disp.min = 0.1,
                             disp.max = 2.0,
                             group.colors = c("hotpink1", "deepskyblue", "violetred4","olivedrab3", "orangered1", "purple3","lightgray"),
                             group.bar.height = 0.05)
heatmaphSS30DIV + scale_fill_gradientn(colors = c("gray98", "red"))
ggsave("/work/nac970/ARXWork2023/hSS30D_Reanalysis_11_02_2023/hSS30DIV_Cleaned_Heatmap_All_Genes_02_22_2024.tiff", width = 25, height = 40, units = c("in"), dpi = 300, bg= "white")
dev.off()


####################################
#Heatmap for hGEO120DIV dataset
#All the genes, cleaned up heatmap 
####################################

#Open hSS120DIV Cell Type File
hSS120D <- readRDS(file = "/work/nac970/ARXWork2023/Scripts_ARX_Revision_2023/hSS120D_Celltype_11_02_2023.rds")
DimPlot(hSS120D)

#Subset all the controls for making heatmap
hSS120DCtrl <- subset(hSS120D, subset = OrgIdent == "Ctrl")

#Subsample 2000 cells from control subset for making heatmap
subsampledhSS120DIV <- hSS120DCtrl[, sample(colnames(hSS120DCtrl), size =2000, replace=F)]

#Preferred order of cell types
a= levels(subsampledhSS120DIV)
mylevels4 =  levels(subsampledhSS120DIV)[c(2,4,5,1,6,3)]
Idents(subsampledhSS120DIV) <- factor(Idents(subsampledhSS120DIV), levels = mylevels4)

#Genes used in the heatmap
genes4 = c("SLC1A3", "MGST1", "ATP1A2", "SPARCL1", "AQP4", "TTYH1", "IGFBP7", "PON2", "SPARC", "BCAN",
           "HES5", "PLPP3","METRN", "VCAM1","LYPD1", "PTPRZ1", "ADGRV1","PSAT1","FABP7","GJA1",
           "ZFP36L1","HES1","BHLHE41","IRX1","NUPR1","PMP2","EDNRB","TNC","SPON1","ANXA5","TIMP3","NTRK2",
           "CD99","ZFP36L2","SLC2A1","NFE2L2","PTGDS","MT2A","GNG5","CRYAB","IGFBP5","CHPF","VIM","NDUFA4L2",
           "CD9","MT1X","NKX6-2","ANGPTL4","SCD","MT3","NUSAP1","TOP2A","CENPF","UBE2C","BIRC5","CDK1",
           "PBK","MKI67","SMC4","TPX2","PTTG1", "PCLAF", "NUF2", "ASPM", "CCNB2", "MAD2L1", 
           "CDKN3", "CENPW", "TYMS", "SPC25", "LHX8","NKX2-1", "TMEFF2", "GABRB2", "GABRG2", "GABRA1", "ENC1",
           "SCG5", "RUNX1T1", "PPP1R1A", "CNTN1", "CHGA", "PBX3", "GBX2","DLX6-AS1","DLX5", "GNG3", "SEMA3E","GRIA2","LHX5-AS1",
           "TSHZ2","SEZ6L2","ATP1A3","CELF5","ANK3","ELAVL4",
           "SCN2A","GABRG2","CELF4","CADM1","CALY","PPP1R1A","NSF","FXYD7","NSG2","ZFHX3","INA","TTC9B", "RAB3B",
           "SLITRK4", "TSHZ2", "SCN2A", "CRABP1", "MEF2C", "ERBB4", "AKAP5", "ZEB2", "NXPH1",
           "PLS3", "PDZRN4", "SP9","EPHA5","PDZRN3","ACKR3","GAD2","BCL11B",
           "IFI44","DACT1","SLAIN1",
           "ANGPT2", "CRABP1", "HOMER3", "TAC3", "RAB3IP","DLX2","SOX6", "MTRNR2L12", 
           "ARL4D", "ARX", "CALB2", "SST", "CORT", "GAD1", "ST6GAL2", "NPY", "SATB1", "LHX6", "TAC1", "STXBP6", 
           "ZNF385D", "AC027031.2", "KCNC2", "NETO2", "BTBD8", "VSNL1")

#Making heatmap
heatmaphSS120DIV <- DoHeatmap(subsampledhSS120DIV,
                              features = genes4,
                              assay = 'RNA',
                              group.by = "ident", 
                              slot = "data", 
                              lines.width = 3,
                              disp.min = 0.1,
                              disp.max = 2.0,
                              group.colors = c("hotpink1", "deepskyblue", "violetred4","olivedrab3", "orange", "lightgray"),
                              group.bar.height = 0.05)

heatmaphSS120DIV + scale_fill_gradientn(colors = c("gray98",  "red"))
ggsave("/work/nac970/ARXWork2023/hSS120D_Reanalysis_11_02_2023/hSS120DIV_Heatmap_Cleaned_All_Genes_02_23_24.tiff", width = 25, height = 40, units = c("in"), dpi = 300, bg= "white")
dev.off()

library(Seurat)
library(dplyr)
library(ggplot2)

###########################
#Heatmap for hCS30DIV Radial Glia,
#Have cleaned this by moving genes up and down 
#And also removing some genes which were not adding additional value
###########################
hCS30D_RadialGlia <- readRDS(file = "/work/nac970/ARXWork2023/Scripts_ARX_Revision_2023/RadialGlia_Seurat_Celltype_10_27_2023.rds")
DimPlot(hCS30D_RadialGlia)
hCS30Dctrl1 <- subset(hCS30D_RadialGlia, subset = OrgIdent == "Ctrl")

#If we want to plot equal number of cells in each group
#maxcells  <- min(table(Idents(hCS30Dctrl1)))

subsampled <- hCS30Dctrl1[, sample(colnames(hCS30Dctrl1), size =2000, replace=F)]
y=levels(subsampled)

features <- c("MEG3","NNAT","C8orf59","TLE4","SMS","RAB8B","MEST","SFRP2","HES5","HES6",
              "GAP43","PTPRZ1","MINOS1","LARP7","PON2","FKBP2","CDH2","IFI44L","B3GAT2",
              "UCHL1","EIF4A1","PCSK1N","DMRTA2","EEF1A2","EMX1","SYT1","STMN2","ARX",
              "RTN1","C1orf61","NRXN1","NFIA","GLCE","GABARAP","EFNB2","NME2","SFRP1","SNHG6","PAX3","IRX2",
              "CRABP1","CRNDE","NTRK2","CTGF","NR2F2","WLS","ARL4A","IGDCC3","FGFR3","PTX3",
              "RFX4","PPP1R1A","TCF12","POU3F2","RASAL2","MAPK10","RRM2","PBK","CDCA5",
              "NCAPG","SPC25","ESCO2","TK1","PKMYT1","MYBL2","ASF1B","PCLAF","CDK1","KIF15","CDC45",
              "NUSAP1","FBLN1","TPM1","SHTN1","KCNQ1OT1","CDKN1C","BIRC5",
              "MKI67","CENPF","PTTG1","CDC20","ARL6IP1","CCNB2",
              "AL139246.5","DLK1","AC092958.1","RGS16","GADD45G",
              "DLL1","DLL3","NEUROG2","PRAG1","ASCL1","ELAVL4","ZC3H12C","TMEM158","NIN",
              "GLUL","TAGLN3","STK17A","TFDP2","TRAF4","GPC2","APOE","KRT18","MTRNR2L8","S100A13")

Idents(subsampled)<- factor(Idents(subsampled), levels = y)
#Idents(hCS30Dctrl1) <- factor(Idents(hCS30Dctrl1), levels = mylevels)

heatmaphCS30D_RadialGlia <- DoHeatmap(subsampled,
                           features = features,
                           group.by = "ident", 
                           assay = "RNA",
                           slot = "data", 
                           lines.width = 5,
                           disp.min = 1.0,
                           disp.max = 2.0,
                           group.colors = c("violet", "turquoise2","green4","orange2", "chartreuse", "grey65", "red", "black","navyblue"),
                           group.bar.height = 0.05)

heatmaphCS30D_RadialGlia + scale_fill_gradientn(colors = c("gray98",  "red"))

ggsave("/work/nac970/ARXWork2023/hCS30D_Reanalysis_10_25_2023/Radial_Glia_Subset/hCS30D_RadialGlia_Heatmap_Cleaned_03_04_2024.tiff", width = 25, height = 40, units = c("in"), dpi = 300, bg= "white")
dev.off()
                          


                        

