library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(Signac)

#Load the previously saved hCS30DIV Seurat file 
hCS30D <- readRDS(file = '~/hCS30DIVSeur.rds')

#Find the markers for each cluster using Venice package
markers = Signac::VeniceAllMarkers(hCS30D, only.pos = T)
write.table(markers, file ='~/hCS30DIVClusterMarkers.txt', sep = "\t", quote=F, row.names =T)
            
#Rename the clusters now after understanding the markers
new.cluster.ids <- c("Radial Glia", 
                     "Inhibitory Neurons", 
                     "Cycling Progenitors", 
                     "Undefined", 
                     "Radial Glia", 
                     "Cycling Progenitors",
                     "Undefined", 
                     "Deep Layer Neurons", 
                     "Cajal-Retzius Neurons",
                     "Undefined", "Inhibitory Neurons",
                     "Intermediate Progenitors",
                     "Undefined, 
                     "Undefined",
                     "Undefined", 
                     "Undefined", 
                     "Inhibitory Neurons", 
                     "Cycling Progenitors",
                     "Cycling Progenitors",
                     "Cycling Progenitors",
                     "Undefined")
            
names(new.cluster.ids) <- levels(hCS30D)
            
hCS30D <- RenameIdents(hCS30D, new.cluster.ids)
DimPlot(hCS30D, reduction = "umap", label = TRUE)
            
#Create a column where the identity of the cluster will be added to either Ctrl or Mut
hCS30D$celltype.OrgIdent <- paste(Idents(hCS30D), hCS30D$OrgIdent, sep = "_")
hCS30D$celltype <- Idents(hCS30D)
            
#Save this file for all further analysis 
saveRDS(hCS30D, file = '~/hCS30DIVCelltype.rds')
            
##########################################################
            
#Custom color for UMAPs
            
#UMAP for control dataset
            
#Subset Control dataset first 
hCS30DCtrl <- subset(hCS30D, subset = OrgIdent == "Ctrl")
hCS30DUMAPCtrl <- DimPlot(hCS30Dctrl, reduction = "umap", 
                          cols = c("Radial Glia" = "hotpink1", 
                                   "Inhibitory Neurons" = "olivedrab3", 
                                   "Cycling Progenitors" = "deepskyblue", 
                                   "Deep Layer Neurons" ="gold1",
                                   "Cajal-Retzius Neurons" = "orangered",
                                   "Intermediate Progenitors" = "violetred4",
                                   "Undefined" = "light gray"))
hCS30DUMAPCtrl
            
#Save the file
ggsave('~hCS30DUMAPCtrl.jpeg', dpi = 300)
dev.off()

#UMAP for Control vs ARX Mutant
hCS30CtrlvsARXUMAP <- DimPlot(hCS30D, reduction = "umap", 
                              cols = c("Radial Glia" = "hotpink1", 
                                   "Inhibitory Neurons" = "olivedrab3", 
                                   "Cycling Progenitors" = "deepskyblue", 
                                   "Deep Layer Neurons" ="gold1",
                                   "Cajal-Retzius Neurons" = "orangered",
                                   "Intermediate Progenitors" = "violetred4",
                                   "Undefined" = "light gray"), split.by = "OrgIdent")

hCS30CtrlvsARXUMAP

#Save the file
ggsave('~/hCS30DControlvsARXUMAP.jpeg', width = 8, height = 3, units = c("in"), dpi = 300)
dev.off()
        
#hCS30D Control vs ARX UMAP for each organoid (C1, C2, C3, M1, M2, M3)
hCS30CtrlvsARXUMAPeachline <- DimPlot(hCS30D, reduction = "umap", 
                                      cols = c("Radial Glia" = "hotpink1", 
                                   "Inhibitory Neurons" = "olivedrab3", 
                                   "Cycling Progenitors" = "deepskyblue", 
                                   "Deep Layer Neurons" ="gold1",
                                   "Cajal-Retzius Neurons" = "orangered",
                                   "Intermediate Progenitors" = "violetred4",
                                   "Undefined" = "light gray"), split.by = "orig.ident")
        
hCS30CtrlvsARXUMAPeachline
        
#Save the file
ggsave('~/hCS30DEachLineUMAP.jpeg', width = 20, height = 4, units = c("in"), dpi = 300)
dev.off()

