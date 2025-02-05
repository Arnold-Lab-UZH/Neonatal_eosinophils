########## This code clusters and annotates cells from SI and colon eosinophil enriched datasets ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.3.Functions_annotation.R")

##### read in object 
obj <- readRDS(file = "/scratch/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_colon_SI.rds")

##### pre-processing and clustering 
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj)
print(ElbowPlot(obj)) 
obj <- FindNeighbors(object = obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.2, random.seed = 5, algorithm = 2)
obj <- RunUMAP(obj, dims = 1:15, seed.use = 5)
DimPlot(obj, label = TRUE)

##### Annotation 
### SingleR for automatic annotation 
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj, assay = "RNA"), ref = mouse.se, labels = mouse.se$label.main)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, obj))

### check quality of the clusters 
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")
# cluster 7 has very high percent.mt and low nFeatures --> low Quality 

### DEGs of clusters 
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

### marker gene expression
FeaturePlot(obj, features = "Siglecf") #Eosinophils cluster 0,1,2
FeaturePlot(obj, features = "Epx") #Eosinophils cluster 0,1,2
FeaturePlot(obj, features = "Col1a1") #Fibroblasts cluster 3
FeaturePlot(obj, features = "Cd209a") #DCs cluster 4
FeaturePlot(obj, features = "Notch3") #Stromal cluster 5
FeaturePlot(obj, features = "Tff3") #Epithelial cluster 6
FeaturePlot(obj, features = "C1qc") #Macrophages cluster 8
FeaturePlot(obj, features = "S100a8") #Neutrophils cluster 9
FeaturePlot(obj, features = "Pecam1") #Endothelial cluster 10 

### rename clusters
current.cluster.ids <- c(0:10)
new.cluster.ids <- c("Eosinophils", "Eosinophils","Eosinophils","Fibroblasts","DCs","Stromal","Epithelial", "lowQ","Macrophages","Neutrophils","Endothelial")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE)

##### save objects 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_colon_SI_anno.rds")
