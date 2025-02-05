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

##### preprocess adult_eosSS to have a normalized and scale data slot for v5 Seurat
adult_eosSS <- NormalizeData(adult_eosSS)
adult_eosSS <- FindVariableFeatures(adult_eosSS)
adult_eosSS <- ScaleData(adult_eosSS,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
adult_eosSS <- RunPCA(adult_eosSS)
adult_eosSS <- FindNeighbors(adult_eosSS, reduction = "pca", dims = 1:30)
adult_eosSS <- FindClusters(adult_eosSS, resolution = 2, algorithm = 2)
adult_eosSS <- RunUMAP(adult_eosSS, reduction = "pca", dims = 1:30,return.model=TRUE)
DimPlot(adult_eosSS,reduction = "umap",group.by = "annotation",combine = FALSE, label.size = 2, split.by = "condition")

##### Label transfer 
neo.query <- neo
neo.query <- NormalizeData(neo.query)
anchors <- FindTransferAnchors(reference = adult_eosSS, query = neo.query, dims = 1:30,reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = adult_eosSS$annotation, dims = 1:30)
neo.query <- AddMetaData(neo.query, metadata = predictions)
#MapQuery is a combined function for TransferData(), IntegrateEmbeddings() and ProjectUMAP()
neo.query <- MapQuery(anchorset = anchors, reference = adult_eosSS, query = neo.query,
                      refdata = list(celltype = "annotation"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(adult_eosSS, reduction = "umap", group.by = "annotation", label = TRUE, label.size = 3,pt.size = 0.5,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(neo.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,pt.size = 0.5,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

## prediction ID is sufficiently high, mean 0.82
VlnPlot(neo.query, features = "predicted.celltype.score", group.by = "condition")
mean(neo.query$predicted.celltype.score)

# add a meta.data column combining annotation and predicted.celltype 
neo.query$annotation <- neo.query$predicted.celltype
# merge objects with umap and ref.umap based on the adult_ss 
integrated <- merge(adult_eosSS, neo.query)
integrated[["umap"]] <- merge(adult_eosSS[["umap"]], neo.query[["ref.umap"]])
DimPlot(integrated, group.by = "condition")
DimPlot(integrated,reduction = "umap",split.by  = "condition",combine = FALSE, label.size = 2, group.by = "predicted.celltype")

## join layers to recreate original counts and data layers  
integrated <- JoinLayers(integrated)

##### highlight colon 
Idents(integrated) <- "condition"
neo <- WhichCells(integrated, idents = c("NEO_P14_colon"))
adult <- WhichCells(integrated, idents = c("adult_colon"))
DimPlot(integrated, label=T, group.by="condition", cells.highlight= list(neo, adult), cols.highlight = c( "#270CEF","#EFA40F"), cols= "#A39F9F")

##### highlight SI 
neo <- WhichCells(integrated, idents = c("NEO_P14_small_int"))
adult <- WhichCells(integrated, idents = c("adult_small_int"))
DimPlot(integrated, label=T, group.by="condition", cells.highlight= list(neo, adult), cols.highlight = c( "#270CEF","#EFA40F"), cols= "#A39F9F")

##### highlight blood 
neo <- WhichCells(integrated, idents = c("NEO_P14_blood"))
adult <- WhichCells(integrated, idents = c("adult_blood"))
DimPlot(integrated, label=T, group.by="condition", cells.highlight= list(neo, adult), cols.highlight = c( "#270CEF","#EFA40F"), cols= "#A39F9F")

##### highlight bm 
neo <- WhichCells(integrated, idents = c("NEO_P14_bm"))
adult <- WhichCells(integrated, idents = c("adult_bm"))
DimPlot(integrated, label=T, group.by="condition", cells.highlight= list(neo, adult), cols.highlight = c( "#270CEF","#EFA40F"), cols= "#A39F9F")

##### highlight spleen 
neo <- WhichCells(integrated, idents = c("NEO_P14_spleen"))
adult <- WhichCells(integrated, idents = c("adult_spleen"))
adult_ctrl <- WhichCells(integrated, idents = c("adult_spleen_ctrl"))
DimPlot(integrated, label=T, group.by="condition", cells.highlight= list(neo, adult,adult_ctrl), cols.highlight = c( "#270CEF","#F20CDC","#EFA40F"), cols= "#A39F9F")

##### save R objects 
saveRDS(integrated, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds")


