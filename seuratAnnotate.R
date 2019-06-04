#library("devtools")
#install_local("~/projects/SingleR/SingleR-master")
library("SingleR")
#install.packages("Seurat")
#install.packages("ggplot2")
#install.packages("sctransform")
library("Seurat")
library("ggplot2")
library("sctransform")

args = commandArgs(trailingOnly = TRUE)

# This workflow describes a Seurat clustering with automatic annotation of cell types. It uses Seurat v3 and sctransform, which allows for a better standardization of the workflow.
# First, specify some variables...
data.dir <- args [1]   #"/path/to/filtered_gene_bc_matrices/hg19/" 
species <- "Human" # or "Mouse"

data 	<- Read10X(data.dir = data.dir)
obj 	<- CreateSeuratObject(counts = data)
# Merging different Seurat objects is possible with the merge function

# store mitochondrial percentage in object meta data
obj 	<- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
pdf("Initial_QC.pdf")
p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p1)
dev.off()

# Subset object
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# run sctransform (other vars to regress?)
obj <- SCTransform(obj, vars.to.regress = c("percent.mt", "nUMI"), verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, verbose = FALSE, resolution = c(0.4, 0.8, 1.2))

pdf("InitialDimPlot.pdf")
p2 <- DimPlot(obj, label = TRUE) + NoLegend()
print(p2)
dev.off()

# Add SingleR cell types to Seurat object
singler <- CreateBigSingleRObject(obj@data, 
	annot = NULL,
	xy = obj@meta.data$xy,
	clusters = obj@meta.data$activ.ident,
	project.name = "tmp",
	min.genes = 200,
	technology = "10X", 
	species = species, 
	citation = "",
	ref.list = list(), 
	normalize.gene.length = F, 
	variable.genes = "de",
	fine.tune = T, 
	do.signatures = F, 
	do.main.types = T, 
	reduce.file.size = T,
	numCores = 16)
# Add annotations to original object
seurat.obj@meta.data$singlerMain <- singler$...



save(singler,file="./singler.RData")

