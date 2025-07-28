########## This code clusters and annotates cells from CD45 enriched colon neo P16 after eosinophil depletion at day 5, 7, and 10 (Cre+) and Cre- as a control ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.3.Functions_annotation.R")

##### read in object 
obj <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P16_iDT_CREpos_CREneg.rds")

##### remove low quality cells and doublets 
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

##### Broad cell type annotation 
### SingleR for broad annotation
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj, assay = "RNA"), ref = mouse.se, labels = mouse.se$label.main)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, obj))
# cluster 7: NK, NKT, T
# cluster 10, 12,14: ILCs
# cluster 9: Eosinophils
# cluster 8, 17: Neutrophils
# cluster 6: Monocytes
# cluster 3: B
# cluster 11: Endothelial
# cluster 0,1,5: DCs
# cluster 4, 15: Epithelial 
# cluster 2,13: Macrophages 

### DEGs of clusters 
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### marker gene expression
FeaturePlot(obj, features = "Xcr1") #DCs cluster 5
FeaturePlot(obj, features = "Cd209a") #DCs cluster 5
FeaturePlot(obj, features = "C1qc") #Macrophages cluster 2
FeaturePlot(obj, features = "Ly6c2") #Macrophages cluster 6
FeaturePlot(obj, features = "Ighm") #B cluster 3
FeaturePlot(obj, features = "Tff3") #Epithelial cluster 4,15
FeaturePlot(obj, features = "Cd3e") #T NK cluster 7
FeaturePlot(obj, features = "S100a8") #Neutrophils cluster 8,17
FeaturePlot(obj, features = "Syne1") #Eosinophils cluster 9
FeaturePlot(obj, features = "Il17rb") #ILCs cluster 12, 14
FeaturePlot(obj, features = "Col4a1") #Endothelial cluster 11
FeaturePlot(obj, features = "Nrgn") #EECs cluster 10

### check quality of the clusters 
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")
# low qualtiy cluster: 1, 13 (very high nFeaturs - doublets?), 0, 16 (high expression of cycling genes)

### rename clusters 
current.cluster.ids <- c(0:17)
new.cluster.ids <- c("lowQ","lowQ","Macrophages","B","Epithelial","DCs","Monocytes","T_NKT","Neutrophils","Eosinophils", "EECs","Endothelial","ILCs","lowQ","ILCs",
                     "Epithelial","lowQ","Neutrophils")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

### remove lowQ cells
Idents(obj) <- "annotation"
obj <- subset(obj, idents =  c("Macrophages","B","Epithelial","DCs","Monocytes","T_NKT","Neutrophils","Eosinophils", "ILCs","Endothelial", "EECs"))

DimPlot(obj ,group.by = "annotation", label= TRUE)

##### save objects 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P16_iDT_CREpos_CREneg_anno.rds")
