library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(Signac)

#Load the previously saved hSS120DIV Seurat file 
hSS120D <- readRDS(file = "~/hSS120DIVSeur.rds")

#Find the markers for each cluster using Venice package
markers = Signac::VeniceAllMarkers(hSS120D, only.pos = T)
write.table(markers, file ="~/hSS120DIVClusterMarkers.txt", sep = "\t", quote=F, row.names =T)
            
#Rename the clusters now after understanding the markers
new.cluster.ids <- c("Inhibitory Neurons", 
                     "Inhibitory Neurons", 
                     "Radial Glial Cells", 
                     "Radial Glial Cells", 
                     "Undefined", 
                     "Radial Glial Cells",
                     "Undefined",
                     "Inhibitory Neurons", 
                     "Undefined", 
                     "Inhibitory Neurons",
                     "Cycling Progenitors", 
                     "Ventral Intermediate Progenitor", 
                     "Undefined", 
                     "Oligodendrocyte Progenitor Cells")

names(new.cluster.ids) <- levels(hSS120D)
hSS120D <- RenameIdents(hSS120D, new.cluster.ids)
            
DimPlot(hSS120D, reduction = "umap")
            
#Create a column where the identity of the cluster will be added to either Ctrl or Mut
hSS120D$celltype.OrgIdent <- paste(Idents(hSS120D), hSS120D$OrgIdent, sep = "_")
hSS120D$celltype <- Idents(hSS120D)
            
#Save the file 
saveRDS(hSS120D, file = "~/hSS120DIVCelltype.rds")

#########################################################            

#Custom color for UMAPs
            
#UMAP for control dataset
            
#Subset Control dataset first
hSS120DCtrl <- subset(hSS120D, subset = OrgIdent == "Ctrl")
hSS120DCtrlUMAP <- DimPlot(hSS120DCtrl, reduction = "umap",
                           cols= c("Ventral Intermediate Progenitors" = "violetred4", 
                                   "Radial Glial Cells" = "hotpink1",
                                   "Cycling Progenitors" = "deepskyblue", 
                                   "Inhibitory Neurons"= "olivedrab3",
                                   "Undefined" = "light gray", 
                                   "Oligodendrocyte Progenitor Cells" = "orange"))
            
hSS120DCtrlUMAP
#Save the file
ggsave("~/hSS120D_UMAP_Ctrl.jpeg", dpi = 300)
dev.off()

#UMAP for Control vs ARX Mutant
hSS120D_UMAP_CtrlvsARX <- DimPlot(hSS120D, reduction = "umap",
                                  cols= c("Ventral Intermediate Progenitors" = "violetred4", 
                                   "Radial Glial Cells" = "hotpink1",
                                   "Cycling Progenitors" = "deepskyblue", 
                                   "Inhibitory Neurons"= "olivedrab3",
                                   "Undefined" = "light gray", 
                                   "Oligodendrocyte Progenitor Cells" = "orange"), split.by = "OrgIdent")
hSS120D_UMAP_CtrlvsARX
ggsave("~/hSS120DUMAPCtrlvsARX.jpeg", width = 9, height = 3, units = c("in"),  dpi = 300)
dev.off()

#UMAP EachLine
hSS120D_UMAP_EachLine <- DimPlot(hSS120D, reduction = "umap",
                                 cols= c("Ventral Intermediate Progenitors" = "violetred4", 
                                   "Radial Glial Cells" = "hotpink1",
                                   "Cycling Progenitors" = "deepskyblue", 
                                   "Inhibitory Neurons"= "olivedrab3",
                                   "Undefined" = "light gray", 
                                   "Oligodendrocyte Progenitor Cells" = "orange"), split.by = "orig.ident")
hSS120D_UMAP_EachLine
ggsave("~/hSS120DUMAPEachLine.jpeg", width = 23, height = 4, units = c("in"),  dpi = 300)
dev.off()
