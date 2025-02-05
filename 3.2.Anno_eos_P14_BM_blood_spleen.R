########## This code clusters and annotates cells from blood, BM and spleen eosinophil enriched datasets ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.3.Functions_annotation.R")

##### read in object 
obj <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_blood_bone_marrow_spleen.rds")

##### remove low quality cells 
VlnPlot(obj, features = "percent.mt")
obj <- subset(obj, subset = percent.mt < 30)

##### pre-processing and clustering 
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj)
print(ElbowPlot(obj)) 
obj <- FindNeighbors(object = obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 5, algorithm = 1)
obj <- RunUMAP(obj, dims = 1:10, seed.use = 5)
DimPlot(obj, label = TRUE, group.by = "condition")
DimPlot(obj, label = TRUE, group.by = "seurat_clusters", label.size = 10)

##### Annotation 
### SingleR for automatic annotation 
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj, assay = "RNA"), ref = mouse.se, labels = mouse.se$label.main)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, obj))
#cluster 13: ILCs/NK/T
#cluster 12: DCs
#cluster 1,6,9: Eosinophils
#cluster 0,2,3,4,8,11: Neutrophils
#cluster 5,7: B cells 
#cluster 10: Monocytes 

### check quality of the clusters 
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")

### DEGs of clusters 
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

### marker gene expression
FeaturePlot(obj, features = "S100a8") #Neutrophils cluster 1,2,3,4,8,11
FeaturePlot(obj, features = "Siglecf") #Eosinophils cluster 1,6,9
FeaturePlot(obj, features = "Epx") #Eosinophils progenitor cluster 9
FeaturePlot(obj, features = "Ighm") #B cells cluster 5,7 
FeaturePlot(obj, features = "Cx3cr1") #Monocytes cluster 10
FeaturePlot(obj, features = "Ccr9") #DCs cluster 12
FeaturePlot(obj, features = "Cd3e") #T/ILC/NK cluster 13

### rename clusters
current.cluster.ids <- c(0:13)
new.cluster.ids <- c("Neutrophils", "Eosinophils","Neutrophils","Neutrophils","Neutrophils","B","Eosinophils", 
                     "B","Neutrophils","Eosinophils","Monocytes","Neutrophils","DCs","T_ILC")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE)

##### save objects 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_blood_bone_marrow_spleen_annotated.rds")

