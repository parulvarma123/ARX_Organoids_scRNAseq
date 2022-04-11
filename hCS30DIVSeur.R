library(Seurat)
library(dplyr)
library(harmony)
library(stringr)

#Preliminary Seurat analysis for Control vs ARX Mutant 30DIV human cortical spheroids (30DIV hCS)

#Create Seurat Object for Control1 
hCS30DBB9data_dir <- '~/hCSControl30d/outs/filtered_feature_bc_matrix'
hCS30DBB9data.data <- Read10X(data.dir = hCS30DBB9data_dir)
hCS30DBB9 <- CreateSeuratObject(counts = hCS30DBB9data.data, project = "C1", min.cells = 10, min.features = 500)
hCS30DBB9$OrgIdent <- "Ctrl"

#Create Seurat object for Control2
hCS30DKB12data_dir <- '~/hCSKB1230D/outs/filtered_feature_bc_matrix'
hCS30DKB12data.data <- Read10X(data.dir = hCS30DKB12data_dir)
hCS30DKB12 <- CreateSeuratObject(counts = hCS30DKB12data.data, project = "C2", min.cells = 10, min.features = 500)
hCS30DKB12$OrgIdent <- "Ctrl"

#Create Seurat Object for Control3
hCS30DDT10data_dir <- '~/hCSDT1030D/outs/filtered_feature_bc_matrix'
hCS30DDT10data.data <- Read10X(data.dir = hCS30DDT10data_dir)
hCS30DDT10 <- CreateSeuratObject(counts = hCS30DDT10data.data, project = "C3", min.cells = 10, min.features = 500)
hCS30DDT10$OrgIdent <- "Ctrl"

#Create Seurat Object for ARX Mutant1
hCS30DEF10data_dir <- '~/30DhcsEF10/outs/filtered_feature_bc_matrix'
hCS30DEF10data.data <- Read10X(data.dir = hCS30DEF10data_dir)
hCS30DEF10 <- CreateSeuratObject(counts = hCS30DEF10data.data, project = "M1", min.cells = 10, min.features = 500)
hCS30DEF10$OrgIdent <- "Mut"

#Create Seurat object for ARX Mutant2
hCS30DAB9data_dir <- '~/hCSAB930D/outs/filtered_feature_bc_matrix'
hCS30DAB9data.data <- Read10X(data.dir = hCS30DAB9data_dir)
hCS30DAB9 <- CreateSeuratObject(counts = hCS30DAB9data.data, project = "M2", min.cells = 10, min.features = 500)
hCS30DAB9$OrgIdent <- "Mut"

#Create Seurat Object for ARX Mutant3
hCS30D301012data_dir <- '~/30Dhcs30101-2/outs/filtered_feature_bc_matrix'
hCS30D301012data.data <- Read10X(data.dir = hCS30D301012data_dir)
hCS30D301012 <- CreateSeuratObject(counts = hCS30D301012data.data, project = "M3", min.cells = 10, min.features = 500)
hCS30D301012$OrgIdent <- "Mut"

#Merge into one dataset 
hCS30Dallreplicates <- merge(hCS30DBB9,  y= c(hCS30DKB12, hCS30DDT10, hCS30DEF10,hCS30DAB9, hCS30D301012), add.cell.ids = c("C1", "C2", "C3", "M1", "M2", "M3"), project = "hCS30AllReplicates")

#QC and removing mitochondrial genes and ribosomal genes as well 
mito.genes <- grep(pattern = "^MT", x = rownames(x=hCS30Dallreplicates), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(hCS30Dallreplicates, slot = 'counts')[mito.genes, ])/Matrix::colSums(GetAssayData(hCS30Dallreplicates, slot="counts"))
hCS30Dallreplicates[['percent.mito']] <- percent.mito

#Stashing away ribo genes 
ribo.genes <- grep("^RP[S,L]", rownames(hCS30Dallreplicates), value = TRUE)
percent.ribo <- Matrix::colSums(GetAssayData(hCS30Dallreplicates, slot ='counts')[ribo.genes, ])/Matrix::colSums(GetAssayData(hCS30Dallreplicates, slot="counts"))
hCS30Dallreplicates[['percent.ribo']] <- percent.ribo

#Combined QC plots
VlnPlot(hCS30Dallreplicates, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol= 4)
plot1 <- FeatureScatter(hCS30Dallreplicates, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(hCS30Dallreplicates, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(hCS30Dallreplicates, feature1 = "nCount_RNA", feature2 = "percent.ribo")
CombinePlots(plots=list(plot1, plot2, plot3))

#Subetting it with 6000 features based on QC
hCS30Dallreplicates <- subset(hCS30Dallreplicates, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < 0.2)

#Normalization using scaling factor 10^6
hCS30Dallreplicates <- NormalizeData(hCS30Dallreplicates, normalization.method = "LogNormalize", scale.factor = 1000000)

# Get variable genes 
hCS30Dallreplicates <- FindVariableFeatures(hCS30Dallreplicates, selection.method = 'mean.var.plot')

top10 <- head(VariableFeatures(hCS30Dallreplicates), 10)
plot1 <- VariableFeaturePlot(hCS30Dallreplicates)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(hCS30Dallreplicates)

#Scale data and regressing out unwanted sources of variation such as mitochondrial genes and ribosomal genes
hCS30Dallreplicates <- ScaleData(hCS30Dallreplicates, vars.to.regress = c("percent.mito", "percent.ribo"))

hCS30Dallreplicates <- RunPCA(hCS30Dallreplicates, features = VariableFeatures(hCS30Dallreplicates), verbose = F)
DimPlot(hCS30Dallreplicates, reduction= "pca")

#Running harmony for batch correction
hCS30Dallreplicates <- RunHarmony(hCS30Dallreplicates, group.by.vars = "orig.ident")

#Running UMAP with harmony embeddings
hCS30Dallreplicates<-  RunUMAP(hCS30Dallreplicates, reduction  = "harmony", dims = 1:30)

#Resolution 0.5 seems to be better 
hCS30Dallreplicates <- FindNeighbors(hCS30Dallreplicates, reduction = "harmony", dims = 1:30)
hCS30Dallreplicates <- FindClusters(hCS30Dallreplicates, resolution= 0.5)
DimPlot(hCS30Dallreplicates, reduction = "umap", label = TRUE)

#Save the Seurat file
saveRDS(file = "~/hCS30DSeur.rds")

