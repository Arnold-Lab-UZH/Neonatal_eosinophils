########## This code integrates and annotates CB and PB eosinophils ##########
# Figure 2 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")

##### read in object 
CB <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Cord_blood_anno.rds")
PB <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/peripheral_blood_anno.rds")

### extract eosinophils 
Idents(CB) <- "annotation"
#EPs = eosinophil progenitors 
CB <- subset(CB, idents = c("Eosinophils","EPs"))

Idents(PB) <- "annotation"
PB <- subset(PB, idents = "Eosinophils")

### add ident with type of blood 
current.cluster.ids <- c("PB7","CB1","CB2","CB3","CB4")
new.cluster.ids <- c("PB","CB","CB","CB","CB")
CB$type <- plyr::mapvalues(x = CB$condition, from = current.cluster.ids, to = new.cluster.ids)

PB$type <- "PB"

##### merge 
obj <- merge(CB, PB)
obj <- JoinLayers(obj)

##### pre-processing
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

##### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:8)
obj <- FindClusters(obj, resolution = 0.1, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:8, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",split.by  = "type",raster=FALSE, group.by = "mnn.clusters")
DimPlot(obj,reduction = "umap.mnn",split.by  = "condition",raster=FALSE, group.by = "mnn.clusters")
# 2 clusters 

### check quality of samples --> sufficiently good 
Idents(obj) <- "condition"
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
VlnPlot(obj, features = "percent.mt", pt.size = 0)

##### Annotation of the two clusters 
FeaturePlot(obj, features = c("EPX", "NFKB2","CLC","AREG",
                              "EPX","PRG2", "CLC","AREG", "ALOX15","PRSS33","PIK3R6","SMPD3"),raster=FALSE,reduction = "umap.mnn")

Idents(obj) <- "annotation"
DotPlot(obj, features = c("EPX","PRG2", "CLC","AREG", "ALOX15","PRSS33","PIK3R6","SMPD3","NFKB2" ), scale = FALSE)
# cluster 0 expresses B-like markers like the peripheral blood (PRSS33, PIK3R6, SMPD3)
# cluster 1 expressed EPX, PRG2 --> precursor markes 

obj <- JoinLayers(obj)

### rename the clusters
current.cluster.ids <- c(0,1)
new.cluster.ids <- c("B_Eos_like","Precursors")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)

p <- DimPlot(obj,reduction = "umap.mnn",group.by = "annotation",raster=FALSE, label=TRUE, cols = c("#2FB34B","#7093CC"))
ggsave(paste0("/scratch/khandl/Neonatal_eosinophils/figures/CB_PB_eos/annotated.svg"), width = 8, height = 5, plot = p)

##### eosinophil markers 
markers <- c("CLC","CCR3","SYNE1","ADGRE5","SIGLEC8","PRSS33","PIK3R6","SMPD3","EPX","PRG2")
for(i in markers) {
  p <- FeaturePlot(obj, features = i, reduction = "umap.mnn", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
  ggsave(paste0("/scratch/khandl/Neonatal_eosinophils/figures/CB_PB_eos/",i,".pdf"), width = 8, height = 5, plot = p)
} 

##### save object 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/CB_PB_eos_integrated_anno.rds")

