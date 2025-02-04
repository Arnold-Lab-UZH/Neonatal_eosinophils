########## This code clusters and annotates cells from cord and peripheral blood ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.Packages_and_functions.R")

########## peripheral blood experiment 1 
##### read in object 
obj <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/peripheral_blood.rds")

##### pre-processing and clustering 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

##### fastMNN integration to account for batch effects 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
DimPlot(obj, group.by = "mnn.clusters",raster=FALSE,reduction = "umap.mnn", label = TRUE)

obj <- JoinLayers(obj)

### DEGs of clusters 
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### marker gene expression
FeaturePlot(obj, features = "CLC") # Eosinophils cluster 0,1,2,3
FeaturePlot(obj, features = "CPA3") # Mast cells cluster 6
FeaturePlot(obj, features = "S100A8") #Neutrophils cluster 5
FeaturePlot(obj, features = "CX3CR1") #Monocytes cluster 4
FeaturePlot(obj, features = "GZMB") #Monocytes cluster 11

VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")
#cluster 4, 9, 10, 12 very low nFeatures, cluster 7 high mito, 

# rename
current.cluster.ids <- c(0:12)
new.cluster.ids <- c("Eosinophils","Eosinophils","Eosinophils","Eosinophils",
                     "lowQ","Neutrophils","Mast","other","Monocytes","lowQ","lowQ",
                     "T","other")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("Eosinophils","Neutrophils","T","Monocytes", "Mast"))

# save object 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/peripheral_blood_anno.rds")

########## cord blood and one peripheral blood control 
##### load R object 
obj <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Cord_blood.rds")

##### clustering 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

##### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.3, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
DimPlot(obj, group.by = "mnn.clusters",raster=FALSE,reduction = "umap.mnn", label = TRUE)

obj <- JoinLayers(obj)

### DEGs of clusters 
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### marker gene expression
FeaturePlot(obj, features = "CLC") # Eosinophils cluster 0
FeaturePlot(obj, features = "ELANE") #Granulocyte-monocyte progenitor (GMP) cluster 2 
FeaturePlot(obj, features = "HBB") # Red blood cells cluster 3, 5
FeaturePlot(obj, features = "EPX") #EPs (eosinophil progenitors) parts of cluster 7
FeaturePlot(obj, features = "S100A8") #Neutrophils cluster 4

VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")

#rename
current.cluster.ids <- c(0:8)
new.cluster.ids <- c("Eosinophils","unknown","GMPs","RBCs",
                     "Neutrophils","RBCs","lowQ","EPs","cycling_unknown")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

# subcluster EPs to remove the contamination 
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "EPs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

# rename clusters
current.cluster.ids <- c("cycling_unknown","Eosinophils","EPs_0","EPs_1","GMPs","lowQ", "Neutrophils",
                         "RBCs","unknown")
new.cluster.ids <- c("cycling_unknown","Eosinophils","lowQ","EPs","GMPs","lowQ", "Neutrophils",
                     "RBCs","unknown")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE, reduction = "umap.mnn")

Idents(subCl) <- "annotation"
obj <- subset(subCl, idents = c("Eosinophils","Neutrophils","EPs","GMPs"))
DimPlot(obj, group.by = "annotation", label = TRUE, reduction = "umap.mnn")

# save object 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/Cord_blood_anno.rds")

