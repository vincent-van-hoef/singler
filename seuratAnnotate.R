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
data.dirs <- c(args[1],args[2], args[3], args[4], args[5])[length(args)] #"/path/to/filtered_gene_bc_matrices/hg19/" 
names(data.dirs) <- c("samp1", "samp2", "samp3", "samp4", "samp5")[length(args)]
species <- "Human" # or "Mouse"

data 	<- Read10X(data.dir = data.dirs)
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
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# run sctransform (other vars to regress?)
obj <- SCTransform(obj, vars.to.regress = c("percent.mt"), verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, verbose = FALSE, resolution = c(0.4, 0.8, 1.2))

pdf("InitialDimPlot.pdf")
p2 <- DimPlot(obj, label = TRUE) + NoLegend()
print(p2)
dev.off()

# Single Function for big objects
CreateBigSingleRObject = function(counts,annot=NULL,project.name,xy,clusters,N=10000,
                                  min.genes=200,technology='10X',
                                  species='Human',citation='',
                                  ref.list=list(),normalize.gene.length=F,
                                  variable.genes='de',fine.tune=T,
                                  reduce.file.size=T,do.signatures=F,
                                  do.main.types=T,
                                  temp.dir=getwd(), numCores = SingleR.numCores) {
  
  n = ncol(counts)
  s = seq(1,n,by=N)
  dir.create(paste0(temp.dir,'/singler.temp/'), showWarnings = FALSE)
  for (i in s) {
    print(i)
    A = seq(i,min(i+N-1,n))
    singler = CreateSinglerObject(counts[,A], annot = annot[A], project.name=project.name, 
                                  min.genes = min.genes,  technology = technology, 
                                  species = species, citation = citation,
                                  do.signatures = do.signatures, clusters = NULL,
                                  numCores = numCores)
    
    save(singler,file=paste0(temp.dir,'/singler.temp/',project.name,'.',i,'.RData'))
  }
  
  singler.objects.file <- list.files(paste0(temp.dir,'/singler.temp/'), 
                                     pattern='RData',full.names=T)
  
  singler.objects = list()
  for (i in 1:length(singler.objects.file)) {
    load(singler.objects.file[[i]])
    singler.objects[[i]] = singler
  }
  
  singler = SingleR.Combine(singler.objects,order = colnames(counts), 
                            clusters=clusters,xy=xy)
  
  singler
}


# Add SingleR cell types to Seurat object
singler <- CreateBigSingleRObject(GetAssayData(obj), 
	annot = NULL,
	xy = NULL,
	clusters = Idents(obj),
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
seurat.obj@meta.data$singlerMain <- singler$singler[[1]]$SingleR.single.main


save(singler,file="./singler.RData")
save(seurat.obj,file="./seurat_anot.RData")
