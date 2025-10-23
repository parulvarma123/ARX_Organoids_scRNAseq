library(Seurat)
library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggridges)
library(dplyr)
library(udunits2)

#Loading hCS30D dataset 

hCS30D <- readRDS(file = "~/hCS30DIVCelltype10_25_2023.rds")
hCS30D <- UpdateSeuratObject(hCS30D)
DimPlot(hCS30D, reduction = "umap")

seur_subset <- subset(seur, idents = c("Radial Glia", "Cycling Progenitors", "Intermediate Progenitors", "Deep Layer Neurons"))
seur_subset <- NormalizeData(seur_subset)
seur_subset <- FindVariableFeatures(seur_subset)
seur_subset <- ScaleData(seur_subset)
seur_subset <- FindNeighbors(seur_subset, dims = 1:30)
seur_subset <- FindClusters(seur_subset, resolution = 0.5)
seur_subset <- RunUMAP(seur_subset, dims = 1:30, n.neighbors = 50)

DimPlot(seur_subset, label = TRUE)

#Create the CDS object 

fd = data.frame("gene_short_name" = rownames(GetAssayData(seur_subset)))
rownames(fd) = rownames(GetAssayData(seur_subset))
cds <- new_cell_data_set(GetAssayData(seur_subset,slot = "counts"),
                         cell_metadata = seur_subset@meta.data,
                         gene_metadata = fd)

#subset so control cells equals mutant cells 
table(colData(cds)$OrgIdent) # look at which one has fewer cells and downsample the other to that number 
ctrlcells = rownames(colData(cds)[colData(cds)$OrgIdent=="Ctrl",])
mutcells = sample(rownames(colData(cds)[colData(cds)$OrgIdent =="Mut",]),8755 ,replace = F)
cds = cds[,c(ctrlcells, mutcells)]

#Run Monocle
cds <- preprocess_cds(cds, num_dim = 30) #use the same # of top PCs as used for clustering
cds <- reduce_dimension(cds,umap.fast_sgd=TRUE)
cds <- cluster_cells(cds, partition_qval = 0.5) #force few partitions
cds <- learn_graph(cds)

levels = c("Radial Glia", "Cycling Progenitors", "Deep Layer Neurons", "Intermediate Progenitors")
cols= c("hotpink1", "deepskyblue", "gold1", "violetred4")
cols = cols[match(levels(factor(colData(seur_new)$celltype)),levels)]
p1 <- plot_cells(seur_new, color_cells_by = "celltype", 
                 cell_size = 0.5,
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 labels_per_group = 0) + scale_color_manual(values=cols) + theme_void()
p1
ggsave("~/hCS30D_Subset_Cells_Plot1_Celltype_Label.jpeg",width = 25, height = 15, units = c("in"), dpi = 300)
dev.off()


p2 = plot_cells(seur_new, color_cells_by = "seurat_clusters", show_trajectory_graph = F,
                label_groups_by_cluster = F, group_label_size = 4,  cell_size = 0.5) +theme_void()
p2
ggsave("~/hCS30D_Subset_Cells_Plot2_SeuratClusters.jpeg",width = 25, height = 15, units = c("in"), dpi = 300)
dev.off()

p3=plot_cells(seur_new, color_cells_by = "OrgIdent", show_trajectory_graph = F, label_cell_groups = F,cell_size = 0.5) + theme_void()
p3
p3+scale_color_manual(values = c("green3", "orange"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave("~/hCS30D_Subset_Cells_Plot3_CtrlvsARX.jpeg",width = 25, height = 15, units = c("in"), dpi = 300)
dev.off()

#Learn graph for new trajectory
seur_new <- learn_graph(seur_new)
seur_new <- order_cells(seur_new)
#Now choose nodes in different partitions
p4 = plot_cells(seur_new,show_trajectory_graph = T,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5,
                cell_size = 0.5) + theme_void()
p4
ggsave("~/hCS30D_Subset_Cells_Plot4_Psuedotime_StartingNode.jpeg",width = 25, height = 15, units = c("in"), dpi = 300)
dev.off()

#All plots together 
png("~/hCS30D_Subset_Cells_AllPlots_1.png",res=125,width=2000,height=400)
plot_grid(p1,p2,p3,p4,ncol=4,rel_widths = c(0.3,0.3,0.2,0.2))
dev.off()

#Converting the cds colData into a dataframe to read the pseudotime values 
#to be plotted on graphs and for doing statistical analysis 

df= as.data.frame(colData(seur_new))
df$pseudotime = pseudotime(seur_new)
df= df[is.finite(df$pseudotime),]

#Trajectory for all cells together
png("~/hCS30D_Psuedotime_Allcells.png", res= 125, height=450, width = 900)
ggplot(df,aes(x=pseudotime, y=OrgIdent))+
  geom_density_ridges(aes(fill=OrgIdent))+ ggtitle("All Cells")+ 
  scale_fill_manual(values = c("green3", "orange"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()

#Trajectory for different cell types 
#Change df$celltype based on which cell type is being analyzed.
#Example of cycling progenitors

png("~/Monocle_CyclingProgenitors.png", res= 125, height=168, width = 400)
ggplot(df[df$celltype == "Cycling Progenitors",],aes(x=pseudotime, y=OrgIdent))+
  geom_density_ridges(aes(fill=OrgIdent))+ ggtitle("Cycling Progenitors")
theme_classic()
dev.off()