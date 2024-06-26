library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(Signac)

#For our own reference 
#C1: hSS30DBB9
#C2: hSS30DKB12
#C3: hSS30DDT10
#M1: hSS30DEF10
#M2: hSS30DAB9
#M3: hSS30D301012

#Seurat V3.2 was used for all the initial analysis
#OrgIdent is Ctrl or Mut
#orig.ident is C1, C2, C3, M1, M2, M3

#Create Seurat Object for Control 1
hSS30DBB9data_dir <- '~/hSSControl30New/outs/filtered_feature_bc_matrix'
hSS30DBB9data.data <- Read10X(data.dir = hSS30DBB9data_dir)
hSS30DBB9 <- CreateSeuratObject(counts = hSS30DBB9data.data, project = "C1", min.cells = 10, min.features = 500)
hSS30DBB9$OrgIdent <- "Ctrl"

#Create Seurat object for Contro 2
hSS30DKB12data_dir <- '~/hSSKB1230D/outs/filtered_feature_bc_matrix'
hSS30DKB12data.data <- Read10X(data.dir = hSS30DKB12data_dir)
hSS30DKB12 <- CreateSeuratObject(counts = hSS30DKB12data.data, project = "C2", min.cells = 10, min.features = 500)
hSS30DKB12$OrgIdent <- "Ctrl"

#Create Seurat Object for Control 3
hSS30DDT10data_dir <- '~/hSSDT1030D/outs/filtered_feature_bc_matrix'
hSS30DDT10data.data <- Read10X(data.dir = hSS30DDT10data_dir)
hSS30DDT10 <- CreateSeuratObject(counts = hSS30DDT10data.data, project = "C3", min.cells = 10, min.features = 500)
hSS30DDT10$OrgIdent <- "Ctrl"

#Create Seurat Object for ARX Mutant 1 
hSS30DEF10data_dir <- '~/hssEF1030D/outs/filtered_feature_bc_matrix'
hSS30DEF10data.data <- Read10X(data.dir = hSS30DEF10data_dir)
hSS30DEF10 <- CreateSeuratObject(counts = hSS30DEF10data.data, project = "M1", min.cells = 10, min.features = 500)
hSS30DEF10$OrgIdent <- "Mut"

#Create Seurat object for ARX Mutant 2
hSS30DAB9data_dir <- '~/hSS30DAB9NN/outs/filtered_feature_bc_matrix'
hSS30DAB9data.data <- Read10X(data.dir = hSS30DAB9data_dir)
hSS30DAB9 <- CreateSeuratObject(counts = hSS30DAB9data.data, project = "M2", min.cells = 10, min.features = 500)
hSS30DAB9$OrgIdent <- "Mut"

#Create Seurat Object for ARX Mutant 3
hSS30D301012data_dir <- '~/hSS30101230D/outs/filtered_feature_bc_matrix'
hSS30D301012data.data <- Read10X(data.dir = hSS30D301012data_dir)
hSS30D301012 <- CreateSeuratObject(counts = hSS30D301012data.data, project = "M3", min.cells = 10, min.features = 500)
hSS30D301012$OrgIdent <- "Mut"

#Merge into one dataset 
hSS30Dallreplicates <- merge(hSS30DBB9,  y= c(hSS30DKB12, hSS30DDT10, hSS30DEF10,hSS30DAB9, hSS30D301012), add.cell.ids = c("C1", "C2", "C3", "M1", "M2", "M3"), project = "hSS30AllReplicates")

# QC and removing mitochondrial genes and ribosomal genes as well 
mito.genes <- grep(pattern = "^MT", x = rownames(x=hSS30Dallreplicates), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(hSS30Dallreplicates, slot = 'counts')[mito.genes, ])/Matrix::colSums(GetAssayData(hSS30Dallreplicates, slot="counts"))
hSS30Dallreplicates[['percent.mito']] <- percent.mito

# Stashing away ribo genes 
ribo.genes <- grep("^RP[S,L]", rownames(hSS30Dallreplicates), value = TRUE)
percent.ribo <- Matrix::colSums(GetAssayData(hSS30Dallreplicates, slot ='counts')[ribo.genes, ])/Matrix::colSums(GetAssayData(hSS30Dallreplicates, slot="counts"))
hSS30Dallreplicates[['percent.ribo']] <- percent.ribo

#Combined QC Plots
VlnPlot(hSS30Dallreplicates, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol= 4)
plot1 <- FeatureScatter(hSS30Dallreplicates, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(hSS30Dallreplicates, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(hSS30Dallreplicates, feature1 = "nCount_RNA", feature2 = "percent.ribo")
CombinePlots(plots=list(plot1, plot2, plot3)) 

#Subsetting it with 7500 features based on QC
hSS30Dallreplicates <- subset(hSS30Dallreplicates, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mito < 0.2)

#Normalization using scaling factor 10^6
hSS30Dallreplicates <- NormalizeData(hSS30Dallreplicates, normalization.method = "LogNormalize", scale.factor = 1000000)

#Get variable genes 
hSS30Dallreplicates <- FindVariableFeatures(hSS30Dallreplicates, selection.method = 'mean.var.plot')

top10 <- head(VariableFeatures(hSS30Dallreplicates), 10)
plot1 <- VariableFeaturePlot(hSS30Dallreplicates)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(hSS30Dallreplicates)

#Scale data and regressing out unwanted sources of variation such as mitochondrial genes and ribosomal genes
hSS30Dallreplicates <- ScaleData(hSS30Dallreplicates, vars.to.regress = c("percent.mito", "percent.ribo"))

hSS30Dallreplicates <- RunPCA(hSS30Dallreplicates, features = VariableFeatures(hSS30Dallreplicates), verbose = F)
DimPlot(hSS30Dallreplicates, reduction= "pca")

#Running harmony for batch correction
hSS30Dallreplicates <- RunHarmony(hSS30Dallreplicates, group.by.vars = "orig.ident")

#Running UMAP with harmony embeddings
hSS30Dallreplicates<-  RunUMAP(hSS30Dallreplicates, reduction  = "harmony", dims = 1:30)

#Resolution 0.5 seems to be better 
hSS30Dallreplicates <- FindNeighbors(hSS30Dallreplicates, reduction = "harmony", dims = 1:30)
hSS30Dallreplicates <- FindClusters(hSS30Dallreplicates, resolution= 0.5)
DimPlot(hSS30Dallreplicates, reduction = "umap", label = TRUE)

#Save the Seurat file 
saveRDS(hSS30Dallreplicates, file= "~/hSS30DIVSeur.rds")
