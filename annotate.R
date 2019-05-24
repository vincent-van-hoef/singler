#library("devtools")
#install_local("~/projects/SingleR/SingleR-master")
library("SingleR")

args <- commandArgs(trailingOnly = TRUE)

seurat.obj <- readRDS(args[1])
species     <- args[2]

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
                                  species = species, citation = citation, ref.list = ref.list,
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

if(species == "Human"){
singler <- CreateBigSingleRObject(seurat.obj@data, 
	annot = NULL,
	xy = NULL,
	project.name = "tmp",
	min.genes = 0,
	technology = "10X", 
	species = "Human", 
	citation = "",
	ref.list = list(blueprint_encode, hpca), 
	normalize.gene.length = F, 
	variable.genes = "de",
	fine.tune = T, 
	do.signatures = F, 
	do.main.types = T, 
	reduce.file.size = T,
	clusters = seurat.obj@meta.data$RawClusterNames,
	numCores = 16)
} else if (species == "Mouse") {
singler <- CreateBigSingleRObject(seurat.obj@data, 
	annot = NULL,
	xy = NULL,
	project.name = "tmp",
	min.genes = 0,
	technology = "10X", 
	species = "Mouse", 
	citation = "",
	ref.list = list(immgen, mouse.rnaseq), 
	normalize.gene.length = F, 
	variable.genes = "de",
	fine.tune = T, 
	do.signatures = F, 
	do.main.types = T, 
	reduce.file.size = T,
	clusters = seurat.obj@meta.data$RawClusterNames,
	numCores = 16)
} else {
stop("Please provide a valid species")
}

# Add annotations to original object
seurat.obj@meta.data$singlerMain <- singler$...



save(singler,file="./singler.RData")

