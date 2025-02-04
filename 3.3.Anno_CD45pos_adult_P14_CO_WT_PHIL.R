########## This code clusters and annotates cells from CD45 enriched colon neo P14 and adult PHIL and WT datasets ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.Packages_and_functions.R")

##### load  object 
obj <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL.rds")

##### remove low quality cells 
Idents(obj) <- "condition"
VlnPlot(obj, features = "percent.mt")
obj <- subset(obj, subset = percent.mt < 30)

##### pre-processing and clustering
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(object = obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
ElbowPlot(obj)
#15 PCs are significant 
obj <- FindNeighbors(object = obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 5, algorithm = 2)
obj <- RunUMAP(obj, dims = 1:15, seed.use = 5)
DimPlot(obj, label = TRUE)

##### Broad cell type annotation 
### SingleR automatic annotation 
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj), ref = mouse.se, labels = mouse.se$label.main)
plotScoreHeatmap(results)
cell.types <- unique(results$pruned.labels)
Idents(merged) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, obj))

### rename clusters 
current.cluster.ids <- c(0:24)
new.cluster.ids <- c("Mono_Mac","B","T_ILC_NK","B","B","DCs",
                     "Mono_Mac","T_ILC_NK","DCs","T_ILC_NK","Mono_Mac","DCs","DCs",
                     "Neutrophils", "Eosinophils", "Mono_Mac","Epithelial","T_ILC_NK","B","B","DCs","Fibroblasts",
                     "Endothelial","Mast","DCs")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

##### Sub-clustering for fine grade annotation 
### Mono_Mac 
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Mono_Mac",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mono_Mac_0","Mono_Mac_1","Mono_Mac_2","Mono_Mac_3"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA")

## marker genes from B Kang et al 2020
#continuum: Ly6c high MHCII negative --> Ly6c high MHCII positive, Cx3cr1 intermediate, Il10, Tnfa, Il23 low --> 
#Ly6c int/low, MHCII high, Cx3cr1 high; cluster 1 expresses less Adgre1 compared to the others 
macs_gut_marker <- c("Adgre1" ,"Ly6c1","Ly6c2","Il10","Il23a", #summary of genes that change their expression 
                     "Thbs1", "Ccr2","Tgm2","Emilin2","Ifi27l2a", #Cluster 1 in Kang et al 
                     "Cx3cr1","Dst", "H2-M2","Dnmt3a","Vcam1", #Cluster 2 in Kang et al
                     "Apoe","Ms4a7","Cd63","Hpgds", #Cluster 3 in Kang et al (also Vcam1)
                     "Pgf","Mmp13","Dnase1l3","Acp5", "C1qb","C1qc", #Cluster 4 in Kang et al (also H2-M2)
                     "Ccl7","Ccl2","Pf4","Cd163", #Cluster 6 in Kang et al
                     "Hspa1a","Hes1","Hspa1b","Ccl4","Jun", #Cluster 7 in Kang et al 
                     "Rsad1","Ifit2","Isg15","Ifit1","Cxcl9" #Cluster 11 in Kang et al
)

DotPlot(sub_celltype, features = macs_gut_marker,dot.scale = 10, scale = TRUE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("Kang et al marker") + theme(axis.text.x = element_text(angle = 90)) 

FeaturePlot(sub_celltype, features = c("Ly6c2","Thbs1","Mmp13","F13a1"))

## rename clusters
current.cluster.ids <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Mono_Mac_0",
                         "Mono_Mac_1","Mono_Mac_2","Mono_Mac_3","Mast","Neutrophils","T_ILC_NK")
new.cluster.ids <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Mac1",
                     "Mono_Mac","Mac2","Mac1","Mast","Neutrophils","T_ILC_NK")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

### T, ILCs and NK 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T_ILC_NK",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.3)
DimPlot(subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_ILC_NK_0","T_ILC_NK_1","T_ILC_NK_2","T_ILC_NK_3",
                                         "T_ILC_NK_4","T_ILC_NK_5","T_ILC_NK_6","T_ILC_NK_7"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

VlnPlot(sub_celltype, features = "nFeature_RNA")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000, margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

## check known marker genes 
#https://www.frontiersin.org/articles/10.3389/fmicb.2011.00111/full for Th subsets 
#https://www.nature.com/articles/s41385-020-0265-y
#ILC1 similar to Th1, ILC2 similar to Th2, ILC3 similar to Th17
#T-bet = Tbx21; Cd39 = Entpd1
t_marker <- c("Cd3g", #general marker of T cells 
              "Il17a","Il22", "Entpd1",  #gamma delta T cells don't express CD4 and CD8, share markers with NK cells 
              "Cd4","Foxp3", #Tregs (there could be two subsets, one expressing Il10 and one not)
              "Cd8a", "Pdcd1", #CD8 T cells (also Tbx21 and Gzmb)
              "Ccr7", #should only be high in naive CD8 T cells 
              "Il21", #Th17 (also Il17a)
              "Il10","Il13","Il5","Il4", #Th2
              "Il12a", #Th1
              "Mcpt1","Mcpt2", #Mast cells
              "Tbx21", #ILC1 
              "Gata3", #ILC2 (also Il4 and Il13 and Il5)
              "Rorc",#ILC3 (also Il17, also Il22)
              "Gzmb" ,#NK (also classified as killer ILCs) (also Cd8a)
              "Klrb1","Klrd1","Klrb1c" #NKT cells 
)

DotPlot(sub_celltype, features = t_marker,dot.scale = 10, scale = TRUE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("") + theme(axis.text.x = element_text(angle = 90)) 

FeaturePlot(sub_celltype, features = c("Gata3","Nkg7","Ccr7","Foxp3","Trdc","Rorc","Cd4","Cd8a", "Cd3g"))

## rename clusters
current.cluster.ids <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Mac1",
                         "Mac2","Mast","Mono_Mac", "Neutrophils", "T_ILC_NK_0", "T_ILC_NK_1", "T_ILC_NK_2",
                         "T_ILC_NK_3", "T_ILC_NK_4", "T_ILC_NK_5", "T_ILC_NK_6", "T_ILC_NK_7")
new.cluster.ids  <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Mac1",
                      "Mac2","Mast","Mono_Mac", "Neutrophils", "CD8T_NKT", "ILC2", "CD4_T",
                      "CD8_CCR7_T", "gamma_delta_T", "ILC3", "doublets", "doublets")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

### B cells 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "B",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "B_0","B_1","B_2","B_3","B_4"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

VlnPlot(sub_celltype, features = "nFeature_RNA")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

## check known marker genes
#PD-L2 = Pdcd1lg2, Cd73 = Nt5e
b_marker <- c("Ighd","Cd19", #naive B cells
              "Ighm", "Entpd1",  #mature B cells
              "Cd80","Nt5e","Pdcd1lg2","Ighg1","Igha" #memory B cells
)

DotPlot(sub_celltype, features = b_marker,dot.scale = 10, scale = TRUE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("") + theme(axis.text.x = element_text(angle = 90)) 

FeaturePlot(sub_celltype, features = c("Cd19","Ighm","Ighd","Igha"))

## rename clusters 
current.cluster.ids <- c("B_0","B_1","B_2","B_3","B_4", "CD4_T","CD8_CCR7_T","CD8T_NKT",
                         "DCs","doublets","Endothelial", "Eosinophils","Epithelial","Fibroblasts","gamma_delta_T","ILC2","ILC3",
                         "Mac1","Mac2", "Mast","Mono_Mac", "Neutrophils")
new.cluster.ids <- c("IgM_IgD_matureB","IgA_PC","IgM_IgD_matureB","IgA_PC","doublets", "CD4_T","CD8_CCR7_T","CD8T_NKT",
                     "DCs","doublets","Endothelial", "Eosinophils","Epithelial","Fibroblasts","gamma_delta_T","ILC2","ILC3",
                     "Mac1","Mac2", "Mast","Mono_Mac", "Neutrophils")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

### DCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "DCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "DCs_0","DCs_1","DCs_2","DCs_3","DCs_4","DCs_5","DCs_6","DCs_7"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

VlnPlot(sub_celltype, features = "nFeature_RNA")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

## check known marker genes
markers_dcs <- c("Siglech","Cox6a2", #pDC
                 "Cd209a","Mgl2","Clec10a","Cd7","Tnfsf9","Irf4", #cDC2
                 "Clec9a","Cd24a" ,"Irf8"#cDC1
)

DotPlot(sub_celltype, features = markers_dcs,dot.scale = 10, scale = FALSE) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 15)) +
  ggtitle("") + theme(axis.text.x = element_text(angle = 90)) 

FeaturePlot(sub_celltype, features = c("Siglech","Cd209a","Clec9a","Irf8", "Fscn1"))

## rename clusters
current.cluster.ids <- c("CD4_T","CD8_CCR7_T","CD8T_NKT","DCs_0","DCs_1","DCs_2","DCs_3","DCs_4","DCs_5","DCs_6","DCs_7",
                         "doublets","Endothelial","Eosinophils","Epithelial","Fibroblasts","gamma_delta_T","IgA_PC","IgM_IgD_matureB",
                         "ILC2","ILC3", "Mac1","Mac2", "Mast","Mono_Mac", "Neutrophils")
new.cluster.ids <- c("CD4_T","CD8_CCR7_T","CD8T_NKT","cDC1","cDC2","remove","cDC1","cDC2","remove","act_mig_DCs","pDCs",
                     "doublets","Endothelial","Eosinophils","Epithelial","Fibroblasts","gamma_delta_T","IgA_PC","IgM_IgD_matureB",
                     "ILC2","ILC3", "Mac1","Mac2", "Mast","Mono_Mac", "Neutrophils")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

### Eosinophils
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Eosinophils",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Eosinophils_0","Eosinophils_1"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

VlnPlot(sub_celltype, features = "nFeature_RNA")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

## rename clusters
current.cluster.ids <- c("act_mig_DCs", "CD4_T","CD8_CCR7_T","CD8T_NKT","cDC1","cDC2","doublets", "Endothelial","Eosinophils_0",
                         "Eosinophils_1","Epithelial","Fibroblasts","gamma_delta_T",
                         "IgA_PC","IgM_IgD_matureB", "ILC2","ILC3", "Mono_Mac", "Mac1","Mac2", "Mast","Neutrophils","pDCs", "remove")
new.cluster.ids <- c("act_mig_DCs", "CD4_T","CD8_CCR7_T","CD8T_NKT","cDC1","cDC2","doublets", "Endothelial","Eosinophils",
                     "remove","Epithelial","Fibroblasts","gamma_delta_T",
                     "IgA_PC","IgM_IgD_matureB", "ILC2","ILC3", "Mono_Mac", "Mac1","Mac2", "Mast","Neutrophils","pDCs", "remove")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

### remove doublets and mito high (named "remove")
Idents(subCl) <- "annotation"
subCl <- subset(subCl, idents = c("act_mig_DCs", "CD4_T","CD8_CCR7_T","CD8T_NKT","cDC1","cDC2",
                                  "Eosinophils","Epithelial","Fibroblasts","Endothelial","gamma_delta_T",
                                  "IgA_PC","IgM_IgD_matureB", "ILC2","ILC3", "Mono_Mac", "Mac1","Mac2", "Mast","Neutrophils","pDCs"))

##### save annotated object 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds")
