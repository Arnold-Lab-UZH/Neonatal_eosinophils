########## This code uses label transfer to integrate neo P14 eosinophils from all tissues (CO, SI, blood, BM, spleen) and adult eosinophils from Gurtner et al 2022 ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")

##### read in R objects 
neo_gut <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_colon_SI_anno.rds")
adult_eosSS <- readRDS("/data/khandl/common/Nature_paper_data/eosinophils_steadystate.rds")
neo_bm_spleen_blood <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_blood_bone_marrow_spleen_annotated.rds")

##### extract eosinophils from neo objects
Idents(neo_gut) <- "annotation"
neo_gut <- subset(neo_gut, idents = "Eosinophils")

Idents(neo_bm_spleen_blood) <- "annotation"
neo_bm_spleen_blood <- subset(neo_bm_spleen_blood, idents = "Eosinophils")

##### check mito cutoff in adult_eosSS and apply the same to the neo object 
VlnPlot(adult_eosSS, features = "percent.mt") #30% 
VlnPlot(neo_gut, features = "percent.mt")
VlnPlot(neo_bm_spleen_blood, features = "percent.mt")
neo_gut <- subset(neo_gut, subset = percent.mt < 30)

##### add condition to adult object 
current.cluster.ids <- c("colon", "small intestine")
new.cluster.ids <- c("adult_colon", "adult_small_int")
adult_eosSS$condition <- plyr::mapvalues(x = adult_eosSS$orig.ident, from = current.cluster.ids, to = new.cluster.ids)

##### rename adult_spleen in neo_bm_spleen_blood to adult_spleen_ctrl
current.cluster.ids <- c("adult_spleen", "NEO_P14_blood","NEO_P14_bm","NEO_P14_spleen")
new.cluster.ids <- c("adult_spleen_ctrl", "NEO_P14_blood","NEO_P14_bm","NEO_P14_spleen")
neo_bm_spleen_blood$condition <- plyr::mapvalues(x = neo_bm_spleen_blood$condition, from = current.cluster.ids, to = new.cluster.ids)

##### add prefix adult to conditions from adult eos sample
current.cluster.ids <- c("bonemarrow", "blood","spleen","stomach","adult_small_int","adult_colon")
new.cluster.ids <- c("adult_bm", "adult_blood","adult_spleen","adult_stomach","adult_small_int","adult_colon")
adult_eosSS$condition <- plyr::mapvalues(x = adult_eosSS$condition, from = current.cluster.ids, to = new.cluster.ids)

##### and add an annotation slot that contains the previous seurat_clusters 
adult_eosSS$annotation <- adult_eosSS$seurat_clusters

##### merge neo_gut and neo_others
neo <- merge(neo_gut, neo_bm_spleen_blood)
neo <- JoinLayers(neo)

##### preprocess 
neo <- NormalizeData(neo)
neo <- FindVariableFeatures(neo)
neo <- ScaleData(neo,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
neo <- RunPCA(neo)
ElbowPlot(neo)
neo <- FindNeighbors(neo, reduction = "pca", dims = 1:10)
neo <- FindClusters(neo, resolution = 0.4, algorithm = 2)
neo <- RunUMAP(neo, reduction = "pca", dims = 1:10,return.model=TRUE)
DimPlot(neo,reduction = "umap",group.by = "seurat_clusters",combine = FALSE, label.size = 2, split.by = "condition")

##### marker expression
Idents(neo) <- "seurat_clusters"
DotPlot(neo, features = c("Mki67","Epx","Ear1","Ear2","Prg2",
                               "Alox15","Aldh2","S100a9","S100a6","S100a10","Il5",
                               "Retnla","Ccl9","Il1rl1","Cd24a",
                               "Mmp9","Icosl","Il4","Tgfb1","Pirb","Rara",
                               "Cd80","Ccl3","Il1b","Cxcl2","Cd274","Tnf","Il16"), scale = TRUE,cols = c("white","darkred"), dot.scale = 8) + theme(axis.text.x = element_text(angle = 90)) 

##### apply Eos score based on Gurtner et al 
### read in object 
adult_eosSS <- readRDS("/data/khandl/common/Nature_paper_data/eosinophils_steadystate.rds")
adult_eosSS$annotation <- adult_eosSS$seurat_clusters

### DEGs
adult_eosSS <- NormalizeData(adult_eosSS, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(adult_eosSS) <- "annotation"
markers <- FindAllMarkers(object = adult_eosSS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
markers_df <- as.data.frame(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

obj <- neo
## add the score to the seurat object 
mouse_eos_genes <- list((markers_df[markers_df$cluster %in% "eosinophil progenitors",])$gene)
obj <-AddModuleScore(obj, features= mouse_eos_genes,name = "precursors")
names(x = obj[[]])

mouse_eos_genes <- list((markers_df[markers_df$cluster %in% "intestinal eosinophils",])$gene)
obj <-AddModuleScore(obj, features= mouse_eos_genes,name = "A_eos")
names(x = obj[[]])

mouse_eos_genes <- list((markers_df[markers_df$cluster %in% "basal eosinophils",])$gene)
obj <-AddModuleScore(obj, features= mouse_eos_genes,name = "B_eos")
names(x = obj[[]])

mouse_eos_genes <- list((markers_df[markers_df$cluster %in% "circulating eosinophils",])$gene)
obj <-AddModuleScore(obj, features= mouse_eos_genes,name = "Circ_eos")
names(x = obj[[]])

mouse_eos_genes <- list((markers_df[markers_df$cluster %in% "immature eosinophils",])$gene)
obj <-AddModuleScore(obj, features= mouse_eos_genes,name = "immature")
names(x = obj[[]])

### plot per condition with fixed scale 
Idents(obj) <- "seurat_clusters"
VlnPlot(obj, features = "immature1", pt.size = 0)
VlnPlot(obj, features = "precursors1",pt.size = 0)
VlnPlot(obj, features = "A_eos1",pt.size = 0)
VlnPlot(obj, features = "B_eos1",pt.size = 0)
VlnPlot(obj, features = "Circ_eos1", pt.size = 0)
VlnPlot(obj, features = "Cd274", pt.size = 0)
VlnPlot(obj, features = "Cd80",pt.size = 0)

### DEGs 
current.cluster.ids <- c(0:8)
new.cluster.ids <- c("B_eos_like", "A_eos","B_eos_like","B_eos_like","Precursor","A_eos","B_eos_like",
                     "B_eos_like","Precursor")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =50, wt = avg_log2FC))

     