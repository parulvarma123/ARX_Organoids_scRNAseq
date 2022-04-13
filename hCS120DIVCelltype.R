library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(Signac)

#Load the previously saved hCS120DIV Seurat file 
hCS120D <- readRDS(file = "~/hCS120DIVSeur.rds")

#Find the markers for each cluster using Venice package
markers = Signac::VeniceAllMarkers(hCS120D, only.pos = T)
write.table(markers, file ="~/hCS120DIVClusterMarkers.txt", sep = "\t", quote=F, row.names =T)
            
#Rename the clusters now after understanding the markers
new.cluster.ids <- c("Deep Layer Neurons", 
                     "oRG", 
                     "oRG", 
                     "Upper Layer Neurons", 
                     "Other", 
                     "Deep Layer Neurons", 
                     "Interneuron Progenitor", 
                     "oRG", 
                     "Other", 
                     "Cycling Progenitors", 
                     "Other", 
                     "Intermediate Progenitors",
                     "Radial Glia", 
                     "Interneurons",
                     "Interneurons", 
                     "Other", 
                     "Other", 
                     "Other", 
                     "Other")

names(new.cluster.ids) <- levels(hCS120D)
hCS120D <- RenameIdents(hCS120D, new.cluster.ids)
            
DimPlot(hCS120D, reduction = "umap", label = TRUE)
            
#Create a column where the identity of the cluster will be added to either Ctrl or Mut
hCS120D$celltype.OrgIdent <- paste(Idents(hCS120D), hCS120D$OrgIdent, sep = "_")
hCS120D$celltype <- Idents(hCS120D)
            
#Save this file for all further analysis
saveRDS(hCS120D, file = '~/hCS120DCelltype.rds')
            
#########################################################
            
#Custom color for UMAPs
            
#UMAP for control dataset
            
#Subset Control dataset first
hCS120Dctrl <- subset(hCS120D, subset = OrgIdent == "Ctrl")
hCS120DUMAPCtrl <- DimPlot(hCS120Dctrl, reduction = "umap", 
                           cols = c("Radial Glia" = "hotpink1", 
                                    "Interneuron Progenitor" = "chartreuse", 
                                    "Cycling Progenitors" = "deepskyblue",
                                    "oRG"="navyblue", 
                                    "Deep Layer Neurons" ="gold1",
                                    "Interneurons" = "green4",
                                    "Upper Layer Neurons"= "orangered",
                                    "Intermediate Progenitors" = "violetred4",
                                    "Other"="gray"))
hCS120DUMAPCtrl
            
#Save the file 
ggsave("~/hCS120DUMAPCtrl.jpeg", dpi = 300)
dev.off()

#UMAP for Control vs ARX Mutant
hCS120DUMAPCtrlvsARX <- DimPlot(hCS120D, reduction = "umap", 
                                cols = c("Radial Glia" = "hotpink1",
                                         "Interneuron Progenitor" = "chartreuse"
                                         "Cycling Progenitors" = "deepskyblue",
                                         "oRG"="navyblue", 
                                         "Deep Layer Neurons" ="gold1",
                                         "Interneurons" = "green4",
                                         "Upper Layer Neurons"= "orangered",
                                         "Intermediate Progenitors" = "violetred4",
                                         "Other"="gray"), split.by = "OrgIdent")
hCS120DUMAPCtrlvsARX

#Save the file 
ggsave("~/hCS120DUMAPCtrlvsARX.jpeg",width = 8, height = 4, units = c("in"), dpi = 300)
dev.off()
        
#UMAP for hCS120D Control vs ARX each organoid(C1, C2, C3, M1, M2, M3)
hCS120DUMAPEachLine <- DimPlot(hCS120D, reduction = "umap", 
                               cols = c("Radial Glia" = "hotpink1", 
                                        "Interneuron Progenitor" = "chartreuse", 
                                        "Cycling Progenitors" = "deepskyblue",
                                        "oRG"="navyblue",
                                        "Deep Layer Neurons" ="gold1",
                                        "Interneurons" = "green4",
                                        "Upper Layer Neurons"= "orangered",
                                        "Intermediate Progenitors" = "violetred4",
                                        "Other"="gray"), split.by = "orig.ident")

hCS120DUMAPEachLine
        
#Save the file 
ggsave("~/hCS120DUMAPEachLine.jpeg",width = 20, height = 4, units = c("in"), dpi = 300)
dev.off()