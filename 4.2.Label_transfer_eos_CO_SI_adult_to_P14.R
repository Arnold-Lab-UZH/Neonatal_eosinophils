########## This code uses label transfer to integrate neo P14 CO and SI eosinophils and adult eosinophils from Gurtner et al 2022 ##########
# Figure 3

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.Packages_and_functions.R")

##### read in R objects 
adult_eosSS <- readRDS("/data/khandl/common/Nature_paper_data/eosinophils_steadystate.rds")
neo_eos <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_colon_SI_anno.rds")

##### extract eosinophils from neo_eos
Idents(neo_eos) <- "annotation"
neo_eos <- subset(neo_eos, idents = "Eosinophils")

##### check mito cutoff in adult_eosSS and apply the same to the neo object 
VlnPlot(adult_eosSS, features = "percent.mt") #30% 
VlnPlot(neo_eos, features = "percent.mt")

neo_eos <- subset(neo_eos, subset = percent.mt < 30)
VlnPlot(neo_eos, features = "percent.mt")

##### add condition to adult object 
current.cluster.ids <- c("colon", "small intestine")
new.cluster.ids <- c("adult_colon", "adult_small_int")
adult_eosSS$condition <- plyr::mapvalues(x = adult_eosSS$orig.ident, from = current.cluster.ids, to = new.cluster.ids)

# and add an annotation slot that contains the previous seurat_clusters 
adult_eosSS$annotation <- adult_eosSS$seurat_clusters

##### preprocess adult_eosSS to have a normalized and scale data slot for v5 Seurat
adult_eosSS <- NormalizeData(adult_eosSS)
adult_eosSS <- FindVariableFeatures(adult_eosSS)
adult_eosSS <- ScaleData(adult_eosSS,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
adult_eosSS <- RunPCA(adult_eosSS)
adult_eosSS <- FindNeighbors(adult_eosSS, reduction = "pca", dims = 1:30)
adult_eosSS <- FindClusters(adult_eosSS, resolution = 2, algorithm = 4)
adult_eosSS <- RunUMAP(adult_eosSS, reduction = "pca", dims = 1:30,return.model=TRUE)
DimPlot(adult_eosSS,reduction = "umap",group.by = "annotation",combine = FALSE, label.size = 2, split.by = "condition")

##### Label transfer 
neo.query <- neo_eos
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

## chose to only take predicted ids if they have a score higher than 0.6
#plot prediction score per condition
p <- VlnPlot(neo.query, features = "predicted.celltype.score", group.by = "condition")
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/prediction_score_neo_CO_SI.svg", plot = p, width = 10, height = 6)
neo.query$celltype <- ifelse(neo.query$prediction.score.max > 0.6, neo.query$predicted.id, "non.defined")

#add a meta.data column combining annotation and predicted.celltype 
neo.query$annotation <- neo.query$celltype
#merge objects with umap and ref.umap based on the adult_ss 
integrated <- merge(adult_eosSS, neo.query)
integrated[["umap"]] <- merge(adult_eosSS[["umap"]], neo.query[["ref.umap"]])
DimPlot(integrated, group.by = "condition")
DimPlot(integrated,reduction = "umap",split.by  = "condition",combine = FALSE, label.size = 2, group.by = "predicted.celltype")

## join layers to recreate original counts and data layers  
integrated <- JoinLayers(integrated)

## investigation of undefined cells 
Idents(integrated) <- "annotation"
markers <- FindAllMarkers(object = integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

top10 <- c( "Ear1" ,    "Ear6" ,    "Ly6c2",    "S100a9",   "Ngp" ,     "Ear2",     "Prg2" ,    "Ifitm6",   "Ltf" ,     "Gm8113",  
            "Gm45039" , "Hspa1b",   "Retnla",   "Tppp3"  ,  "Gda"  ,    "Gkn2" ,    "Tff1"  ,   "Xbp1",     "Gkn1" ,  #9 
            "Dusp1"  ,  "Pgc"     ,"Gm38126"  ,"Igkc"     ,"Ier2"     ,"Rara"  ,   "Jun"     , "Jund"    , "Lair1" ,   "Plxnd1" , 
            "Mki67"   , "Birc5"   , "Rrm2"    ,"Pclaf",    "Spc24" ,   "Ccna2"   , "Nuf2" ,    "Uhrf1"  ,  "Asf1b"  ,  "Cdk1" , 
            "Tgm2"  ,   "Thbs1"  ,  "Slc11a1",  "Il1b"   ,  "Acod1"   , "Il1r2"    ,"Uaca"    , "Il1a"    , "Clec4e"   ,"Serpine1", 
            "Pi4ka"   , "Galnt2l" , "Myo1d"   , "Ldlr"    , "Arhgap15", "Nedd9"  ,  "Hbp1"   ,  "Naga"   ,  "Pde4d"  ,  "Lyst")

#heatmap per annotation cluster 
#anno 2 = Harmony, anno3 = fastMNN
heatmap_goi_coi(integrated, "annotation",top10,"RNA", c("immature_eos","circ_eos","B_eos","precursor_eos","A_eos","undefined"), 
                c(10,9,10,10,10,10),c( "#E88A1A","#E8E81A","#10A069","#26DFED","#E81818","#A39F9F"),
                c( immature_eos= "#E88A1A",circ_eos = "#E8E81A",B_eos="#10A069",precursor_eos="#26DFED",A_eos="#E81818",undefined="#A39F9F"),F,F)

p <- VlnPlot(integrated, features = "nFeature_RNA",group.by = "annotation", cols = c("#10A069","#E8E81A","#26DFED","#E88A1A","#E81818","#F408D3"))
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/nFeature.svg", plot = p, width = 5, height = 6)

p <- VlnPlot(integrated, features = "percent.mt",group.by = "annotation",cols = c("#10A069","#E8E81A","#26DFED","#E88A1A","#E81818","#F408D3"))
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/percent.mt.svg", plot = p, width = 5, height = 6)

## rename non-defined clusters as A-eos 
current.cluster.ids <- c("basal eosinophils", "circulating eosinophils","eosinophil progenitors","immature eosinophils",
                         "intestinal eosinophils","non.defined")
new.cluster.ids <- c("B_eos", "circ_eos","precursor_eos","immature_eos",
                     "A_eos","A_eos")
integrated$annotation <- plyr::mapvalues(x = integrated$annotation, from = current.cluster.ids, to = new.cluster.ids)

p <- DimPlot(integrated, group.by = "annotation",cols = c("#E81818","#10A069","#E8E81A", "#E88A1A","#26DFED"))
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/umap_annotated.svg", plot = p, width = 10, height = 6)

p <- DimPlot(integrated, group.by = "condition",cols = c("#A39F9F","#A39F9F","#A39F9F", "#A39F9F","#870C27",
                                                         "#240FF2","#A39F9F","#A39F9F"))
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/umap_neo.svg", plot = p, width = 10, height = 6)

p <- DimPlot(integrated, group.by = "annotation",split.by= "condition",cols = c("#E81818","#10A069","#E8E81A", "#E88A1A","#26DFED"))
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/umap_split_by_cond.svg", plot = p, width = 10, height = 6)

##### save R objects 
saveRDS(integrated, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_LT.rds")




