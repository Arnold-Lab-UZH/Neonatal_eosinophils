########## This code clusters and annotates cells from CD45 negative colon neo P14 PHIL and WT datasets ##########
# Figure S5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.3.Functions_annotation.R")
source("~/Projects/Neonatal_eosinophils/1.4.Functions_DEGs.R")

##### load  object 
obj <- readRDS(file = "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL.rds")

##### remove low quality cells 
VlnPlot(obj, features = "percent.mt")
obj <- subset(obj, subset = percent.mt < 30)

##### pre-processing and clustering
options(future.globals.maxSize = 8000 * 1024^2)
obj <- SCTransform(obj, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(object = obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(object = obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 5, algorithm = 1)
obj <- RunUMAP(obj, dims = 1:10, seed.use = 5)
DimPlot(obj, label = TRUE, split.by = "condition")
DimPlot(obj, label = TRUE)

##### check QC measures per cluster 
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")

##### Cell type annotation 
### SingleR for automatic annotation 
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(neo), ref = mouse.se, labels = mouse.se$label.main)
plotScoreHeatmap(results)
cell.types <- unique(results$pruned.labels)
Idents(neo) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, neo))

### DEG analysis
neo <- NormalizeData(neo, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay = "RNA"
Idents(neo) <- "seurat_clusters"
markers <- FindAllMarkers(object = neo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### plot known marker genes in a heatmap (Sirvinskas et al iScience 2022 https://pubmed.ncbi.nlm.nih.gov/35479413/; and Drokhlyansky et al Cell 2020 https://pubmed.ncbi.nlm.nih.gov/32888429/) 
# for SMCs additionally Ruan et al 2020 (https://www.ahajournals.org/doi/10.1161/ATVBAHA.120.315107)
# for Pericytes additionally: Lu et al Vasc. Health Risk manag. 2023 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10596204/) 
# and Kim et al Nat. Com 2020 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6969052/)
markers <- c("Epcam","Krt8","Krt18","Guca2a","Cldn3", #epithelial 5
             "Pclaf","Clca3b","Mki67","Lgr5", #Stem-TA 4
             "Car1","Cyp2c55","Mpp6","Selenbp1", #Colonocytes 1 4
             "Aqp4","Muc3","Ggh", #Colonocytes 2.   3
             "Ccl6","Spink4","Muc2","Reg4", #Goblet 4
             "Tuba1a","Dclk1","Chga","Chgb", #Tuft-EE 4
             "Elavl3","Elavl4","Snap25","Uchl1","Ret", #neurons 5
             "S100b","Sox10","Erbb3","Plp1","Gfap", #glial 5
             "Cck","Gfra3","Vwa5b2","Neurod1", #EEC 4
             "Lyz1","Lyz2" #Paneth cells 
)

heatmap_goi_coi(obj, "seurat_clusters",markers,"SCT", c("epithelial","StemTA","Colonocytes1","Colonocytes2","Goblet","TuftEE",
                                                 "Neurons","Glial","EEC","Paneth"), c(5,4,4,3,4,4,5,5,4,2),   
                c("#EF0A3C",  "#EF0AA3","#DE0AEF",  "#680AEF","#0A41EF",  "#0ADFEF", "#0AEF97",  "#0AEF40","#5FEF0A","#EF0AA3"),
                c(epithelial="#EF0A3C", StemTA= "#EF0AA3",Colonocytes1="#DE0AEF", Colonocytes2= "#680AEF",Goblet="#0A41EF", 
                  TuftEE= "#0ADFEF",Neurons="#0AEF97", Glial= "#0AEF40",EEC="#5FEF0A",Paneth= "#EF0AA3"),F,F)

### rename clusters 
current.cluster.ids <- c(0:14)
new.cluster.ids <- c("Epithelial_lowQ","Epithelial_lowQ","SMCs","Colonocytes","Fibroblasts","TA",
                     "TA","Goblet","Endothelial","Fibroblasts","lowQ",
                     "Myofibroblasts","Pericytes","EECs","Colonocytes")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(obj, group.by = "annotation",label = TRUE)

##### remove low quality 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("SMCs","Colonocytes",
                              "Fibroblasts","TA","Goblet","Endothelial",
                              "Myofibroblasts","Pericytes","EECs"))
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE)

##### DEG analysis of annotated clusters 
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay = "RNA"
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
write.csv(markers, "/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg_neo_P14_colon_annotation.csv" )

## plot marker genes in DotPlot
markers <- c("Actg2","Cnn1","Acta2","Des", #SMCs
             "Car1","Cyp2c55","Muc3", #Colonocytes
             "Dcn","Pdgfra","Col1a1","Lum", #Fibroblasts
             "Gpx2","Clca3b","Pclaf","Lgr5", #TA/Enterocytes
             "Fcgbp","Tff3","Muc2","Clca1", #Goblet
             "Cdh5","Plvap","Cd36","Fabp4", #Endothelial
             "Rcan2","Crispld2","Adra1a","Map3k7cl","Foxc1", #myofibroblasts
             "Rgs5","Rgs4","Pdgfrb", #Pericytes
             "Chgb","Chga","Neurod1","Gfra3" #EECs
)
Idents(obj) <- "annotation"
p <- DotPlot(obj, features = markers, scale = FALSE,cols = c("white","darkred"), dot.scale = 5) + theme(axis.text.x = element_text(angle = 90)) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/DotPlot_markers_CD45neg_neo_P14_annotated.svg", width = 8, height = 5, plot = p)

##### save object 
saveRDS(obj, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL_anno.rds")
