library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(Signac)

#For our own reference 
#C1: hCS120DBB9
#C2: hCS120DKB12
#C3: hCS120DDT10
#M1: hCS120DEF10
#M2: hCS120DAB9
#M3: hCS120D301012

#Seurat V3.2 was used for all the initial analysis
#OrgIdent is Ctrl or Mut
#orig.ident is C1, C2, C3, M1, M2, M3

#Create Seurat Object for Control 1
hCS120DBB9data_dir <- '~/ControlN/outs/filtered_feature_bc_matrix'
hCS120DBB9data.data <- Read10X(data.dir = hCS120DBB9data_dir)
hCS120DBB9 <- CreateSeuratObject(counts = hCS120DBB9data.data, project = "C1", min.cells = 10, min.features = 500)
hCS120DBB9$OrgIdent <- "Ctrl"

#Create Seurat object for Control 2
hCS120DKB12data_dir <- '~/hCS120DKB12/outs/filtered_feature_bc_matrix'
hCS120DKB12data.data <- Read10X(data.dir = hCS120DKB12data_dir)
hCS120DKB12 <- CreateSeuratObject(counts = hCS120DKB12data.data, project = "C2", min.cells = 10, min.features = 500)
hCS120DKB12$OrgIdent <- "Ctrl"

#Create Seurat Object for Control 3 
hCS120DDT10data_dir <- '~/hCS120DDT10New/outs/filtered_feature_bc_matrix'
hCS120DDT10data.data <- Read10X(data.dir = hCS120DDT10data_dir)
hCS120DDT10 <- CreateSeuratObject(counts = hCS120DDT10data.data, project = "C3", min.cells = 10, min.features = 500)
hCS120DDT10$OrgIdent <- "Ctrl"

#Create Seurat Object for ARX Mutant 1 
hCS120DEF10data_dir <- '~/PatientN/outs/filtered_feature_bc_matrix'
hCS120DEF10data.data <- Read10X(data.dir = hCS120DEF10data_dir)
hCS120DEF10 <- CreateSeuratObject(counts = hCS120DEF10data.data, project = "M1", min.cells = 10, min.features = 500)
hCS120DEF10$OrgIdent <- "Mut"

#Create Seurat object for ARX Mutant 2
hCS120DAB9data_dir <- '~/hCS120DAB9/outs/filtered_feature_bc_matrix'
hCS120DAB9data.data <- Read10X(data.dir = hCS120DAB9data_dir)
hCS120DAB9 <- CreateSeuratObject(counts = hCS120DAB9data.data, project = "M2", min.cells = 10, min.features = 500)
hCS120DAB9$OrgIdent <- "Mut"

#Create Seurat Object for ARX Mutant 3
hCS120D301012data_dir <- '~/hCS120D301012New/outs/filtered_feature_bc_matrix'
hCS120D301012data.data <- Read10X(data.dir = hCS120D301012data_dir)
hCS120D301012 <- CreateSeuratObject(counts = hCS120D301012data.data, project = "M3", min.cells = 10, min.features = 500)
hCS120D301012$OrgIdent <- "Mut"

#Merge into one dataset 
hCS120Dallreplicates <- merge(hCS120DBB9,  y= c(hCS120DKB12, hCS120DDT10, hCS120DEF10,hCS120DAB9, hCS120D301012), add.cell.ids = c("C1", "C2", "C3", "M1", "M2", "M3"), project = "hCS120DAllReplicates")

# QC and removing mitochondrial genes and ribosomal genes as well 
mito.genes <- grep(pattern = "^MT", x = rownames(x=hCS120Dallreplicates), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(hCS120Dallreplicates, slot = 'counts')[mito.genes, ])/Matrix::colSums(GetAssayData(hCS120Dallreplicates, slot="counts"))
hCS120Dallreplicates[['percent.mito']] <- percent.mito

# Stashing away ribo genes 
ribo.genes <- grep("^RP[S,L]", rownames(hCS120Dallreplicates), value = TRUE)
percent.ribo <- Matrix::colSums(GetAssayData(hCS120Dallreplicates, slot ='counts')[ribo.genes, ])/Matrix::colSums(GetAssayData(hCS120Dallreplicates, slot="counts"))
hCS120Dallreplicates[['percent.ribo']] <- percent.ribo

#Combined QC Plots
VlnPlot(hCS120Dallreplicates, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol= 4)
plot1 <- FeatureScatter(hCS120Dallreplicates, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(hCS120Dallreplicates, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(hCS120Dallreplicates, feature1 = "nCount_RNA", feature2 = "percent.ribo")
CombinePlots(plots=list(plot1, plot2, plot3))

# Subsetting it with 7000 features based on QC
hCS120Dallreplicates <- subset(hCS120Dallreplicates, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mito < 0.2)

# Normalization using scaling factor 10^6
hCS120Dallreplicates <- NormalizeData(hCS120Dallreplicates, normalization.method = "LogNormalize", scale.factor = 1000000)

# Get variable genes 
hCS120Dallreplicates <- FindVariableFeatures(hCS120Dallreplicates, selection.method = 'mean.var.plot')

top10 <- head(VariableFeatures(hCS120Dallreplicates), 10)
plot1 <- VariableFeaturePlot(hCS120Dallreplicates)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(hCS120Dallreplicates)

#Scale data and regressing out unwanted sources of variation such as mitochondrial genes and ribosomal genes
hCS120Dallreplicates <- ScaleData(hCS120Dallreplicates, vars.to.regress = c("percent.mito", "percent.ribo"))

hCS120Dallreplicates <- RunPCA(hCS120Dallreplicates, features = VariableFeatures(hCS120Dallreplicates), verbose = F)
DimPlot(hCS120Dallreplicates, reduction= "pca")

# Running harmony for batch correction
hCS120Dallreplicates <- RunHarmony(hCS120Dallreplicates, group.by.vars = "orig.ident")

#Running UMAP with harmony embeddings
hCS120Dallreplicates<-  RunUMAP(hCS120Dallreplicates, reduction  = "harmony", dims = 1:30)

# Resolution 0.5 seems to be better 
hCS120Dallreplicates <- FindNeighbors(hCS120Dallreplicates, reduction = "harmony", dims = 1:30)
hCS120Dallreplicates <- FindClusters(hCS120Dallreplicates, resolution= 0.5)
DimPlot(hCS120Dallreplicates, reduction = "umap", label = TRUE)

#Save the Seurat file
saveRDS(hCS120Dallreplicates, file = "~/hCS120DIVSeur.rds")
