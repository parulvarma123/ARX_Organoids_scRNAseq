library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(Signac)

#Load the previously saved hSS30DIV Seurat file 
hSS30D <- readRDS(file = "~/hSS30DIVSeur.rds")

#Find the markers for each cluster using Venice package
markers = Signac::VeniceAllMarkers(hSS30D, only.pos = T)
write.table(markers, file ="~/hSS30DIVClusterMarkers.txt", sep = "\t", quote=F, row.names =T)
            
#Rename the clusters now after understanding the markers
new.cluster.ids <- c("Ventral Intermediate Progenitor", 
                     "Radial Glial Cells", "Cycling Progenitors", 
                     "Inhibitory Neurons",
                     "Inhibitory Neurons", 
                     "Radial Glial Cells", 
                     "Radial Glial Cells", 
                     "Undefined", 
                     "Undefined", 
                     "Radial Glial Cells",
                     "Cajal-Retzius Neurons",
                     "Radial Glial Cells", 
                     "Undefined", 
                     "Inhibitory Neurons", 
                     "Cycling Progenitors", 
                     "Radial Glial Cells", 
                     "Undefined", 
                     "Microglia", 
                     "Undefined", 
                     "Radial Glial Cells",
                     "Undefined", 
                     "Undefined",
                     "Cycling Progenitors")

names(new.cluster.ids) <- levels(hSS30D)
hSS30D <- RenameIdents(hSS30D, new.cluster.ids)
            
#Create a column where the identity of the cluster will be added to either Ctrl or Mut
hSS30D$celltype.OrgIdent <- paste(Idents(hSS30D), hSS30D$OrgIdent, sep = "_")
hSS30D$celltype <- Idents(hSS30D)
            
#Save this file for all further analysis
saveRDS(hSS30D, file = "~/hSS30DCelltype.rds")
            
#########################################################
            
#Custom color for UMAPs
            
#UMAP for control dataset
            
#Subset Control dataset first
hSS30DCtrl <- subset(hSS30D, subset = OrgIdent == "Ctrl")
hSS30D_UMAP_Ctrl <- DimPlot(hSS30DCtrl, reduction = "umap",
                            cols = c("Ventral Intermediate Progenitors" = "violetred4", 
                                     "Radial Glial Cells" = "hotpink1",
                                     "Cycling Progenitors" = "deepskyblue", 
                                     "Undefined"= "light gray",
                                     "Inhibitory Neurons" = "olivedrab3",
                                     "Cajal-Retzius Neurons" = "orangered",
                                     "Microglia" = "purple3"))
hSS30D_UMAP_Ctrl
            
#Save the file
ggsave("~/hSS30D_UMAP_Ctrl.jpeg", dpi = 300)
dev.off()

#UMAP for Control vs ARX Mutant
hSS30D_UMAP_CtrlvsARX <- DimPlot(hSS30D, reduction = "umap",
                                 cols= c("Ventral Intermediate Progenitors" = "violetred4", 
                                     "Radial Glial Cells" = "hotpink1",
                                     "Cycling Progenitors" = "deepskyblue", 
                                     "Undefined"= "light gray",
                                     "Inhibitory Neurons" = "olivedrab3",
                                     "Cajal-Retzius Neurons" = "orangered",
                                     "Microglia" = "purple3"), split.by = "OrgIdent")
hSS30D_UMAP_CtrlvsARX

#Save the file
ggsave("~/hSS30D_UMAP_CtrlvsARX.jpeg", width = 9, height = 3, units = c("in"),  dpi = 300)
dev.off()
        
#UMAP for hSS30D Control vs ARX each organoid (C1, C2, C3, M1, M2, M3)
hSS30DControlvsARXUMAP_EachLine <- DimPlot(hSS30D, reduction = "umap",
                                           cols= c("Ventral Intermediate Progenitors" = "violetred4", 
                                     "Radial Glial Cells" = "hotpink1",
                                     "Cycling Progenitors" = "deepskyblue", 
                                     "Undefined"= "light gray",
                                     "Inhibitory Neurons" = "olivedrab3",
                                     "Cajal-Retzius Neurons" = "orangered",
                                     "Microglia" = "purple3"), split.by = “orig.ident”)
        
hSS30DControlvsARXUMAP_EachLine
        
#Save the file 
ggsave("~/hSS120DUMAPEachLine.jpeg", width = 23, height = 4, units = c("in"),  dpi = 300)
dev.off()
