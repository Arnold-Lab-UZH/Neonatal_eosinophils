########## This code clusters and annotates cells from CD45 negative colon neo P14 PHIL and WT datasets ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir,"01.05.Functions_annotation.R"))
source(file.path(base_dir,"01.06.Functions_DEGs.R"))

##### load  object 
obj <- readRDS(file = file.path(seurat_objects_dir,"Neo_P14_CD45neg_colon_WT_PHIL.rds"))

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
results <- SingleR(test = as.SingleCellExperiment(obj), ref = mouse.se, labels = mouse.se$label.main)
plotScoreHeatmap(results)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, obj))

### DEG analysis
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay = "RNA"
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,slot = "data")
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

DotPlot(obj, features = markers, scale = FALSE,cols = c("white","darkred"), dot.scale = 5) + theme(axis.text.x = element_text(angle = 90)) 

##### refine annotation based on Brügger et al https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001032
markers <- c("Epcam","Lrg5","Axin2", "Mki67","Atoh1","Muc2", "Reg4", "Chga","Chgb", 
             "Krt19","Guca2a","Alpi", "S100b", "Rgs5","Pecam1","Myh11","Acta2","Vim","Col1a1","Gli1","Foxl1","Cd34","Cd90","Pdgfra"
)
Idents(obj) <- "annotation"
DotPlot(obj, features = markers, scale = TRUE,cols = c("white","darkred"), dot.scale = 5) + theme(axis.text.x = element_text(angle = 90)) 

### rename 
current.cluster.ids <- c("SMCs","Colonocytes","Fibroblasts","TA","Goblet","Endothelial","Myofibroblasts","Pericytes","EECs")
new.cluster.ids <- c("Stromal","Enterocytes","Stromal","TA","Goblet_Reg4","Endothelial","Stromal","Pericytes","EECs")
obj$annotation1 <- plyr::mapvalues(x = obj$annotation, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(obj, group.by = "annotation1",label = TRUE, label.size = 5)

#Pericytes 
FeaturePlot(obj, features = "Rgs5")

### subcluster TA 
Idents(obj) <- "annotation1"
subCl <- FindSubCluster(obj,cluster = "TA",graph.name = "SCT_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c("TA_0","TA_1","TA_2"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster", label.size = 5)
DimPlot(sub_celltype ,label = TRUE, group.by = "sub.cluster", split.by = "condition")

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

FeaturePlot(sub_celltype, features = c("Lgr5","Epcam","Axin2","Mki67")) 

DotPlot(sub_celltype, features = c("Lgr5","Epcam","Axin2","Mki67"),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# rename
current.cluster.ids <- c("EECs","Endothelial","Enterocytes","Goblet_Reg4","Pericytes", "Stromal","TA_0","TA_1","TA_2")
new.cluster.ids <- c("EECs","Endothelial","Enterocytes","Goblet_Reg4","Pericytes", "Stromal","Goblet_Reg4","TA","lowQ")
subCl$annotation1 <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation1", label = TRUE,raster=FALSE)

### subcluster Enterocytes 
Idents(subCl) <- "annotation1"
subCl <- FindSubCluster(subCl,cluster = "Enterocytes",graph.name = "SCT_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c("Enterocytes_0","Enterocytes_1","Enterocytes_2"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster", label.size = 5)
DimPlot(sub_celltype ,label = TRUE, group.by = "sub.cluster", split.by = "condition")

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

FeaturePlot(sub_celltype, features = c("Epcam","Krt19","Guca2a","Alpi")) 

DotPlot(sub_celltype, features = c("Epcam","Krt19","Guca2a","Alpi"),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# rename
current.cluster.ids <- c("EECs","Endothelial","Enterocytes_0","Enterocytes_1","Enterocytes_2","Goblet_Reg4","Pericytes", "Stromal","TA")
new.cluster.ids <- c("EECs","Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Enterocytes_Alpi_high","Goblet_Reg4","Pericytes", "Stromal","TA")
subCl$annotation1 <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation1", label = TRUE,raster=FALSE)

### subcluster Endothelial 
Idents(subCl) <- "annotation1"
subCl <- FindSubCluster(subCl,cluster = "Endothelial",graph.name = "SCT_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c("Endothelial_0","Endothelial_1","Endothelial_2"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster", label.size = 5)
DimPlot(sub_celltype ,label = TRUE, group.by = "sub.cluster", split.by = "condition")

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

FeaturePlot(sub_celltype, features = c("Rgs5","Pecam1","Vim")) 

DotPlot(sub_celltype, features = c("Rgs5","Pecam1","Vim"),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# rename
current.cluster.ids <- c("EECs","Endothelial_0","Endothelial_1","Endothelial_2","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Goblet_Reg4","Pericytes","Stromal","TA")
new.cluster.ids <- c("EECs","Endothelial","Endothelial","Endothelial", "Enterocytes_Alpi_low","Enterocytes_Alpi_high","Goblet_Reg4","Pericytes","Stromal","TA")
subCl$annotation1 <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation1", label = TRUE,raster=FALSE)

### subcluster Goblet_Reg 4 (secretory cells ) 
Idents(subCl) <- "annotation1"
subCl <- FindSubCluster(subCl,cluster = "Goblet_Reg4",graph.name = "SCT_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c("Goblet_Reg4_0","Goblet_Reg4_1","Goblet_Reg4_2","Goblet_Reg4_3"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster", label.size = 5)
DimPlot(sub_celltype ,label = TRUE, group.by = "sub.cluster", split.by = "condition")

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

FeaturePlot(sub_celltype, features = c("Atoh1","Muc2","Reg4")) 

DotPlot(sub_celltype, features = c("Atoh1","Muc2","Reg4"),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# rename
current.cluster.ids <- c("EECs","Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high", "Goblet_Reg4_0","Goblet_Reg4_1","Goblet_Reg4_2",
                         "Goblet_Reg4_3","lowQ","Pericytes","Stromal","TA")
new.cluster.ids <- c("EECs","Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Reg4_secretory","Goblet","Reg4_secretory",
                     "Reg4_secretory","lowQ","Pericytes","Stromal","TA")
subCl$annotation1 <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation1", label = TRUE,raster=FALSE)

### subcluster EECs 
Idents(subCl) <- "annotation1"
subCl <- FindSubCluster(subCl,cluster = "EECs",graph.name = "SCT_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c("EECs_0","EECs_1","EECs_2"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster", label.size = 5)
DimPlot(sub_celltype ,label = TRUE, group.by = "sub.cluster", split.by = "condition")

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

FeaturePlot(sub_celltype, features = c("S100b","Reg4","Chga","Chgb")) 

DotPlot(sub_celltype, features = c("S100b","Reg4","Chga","Chgb"),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# rename
current.cluster.ids <- c("EECs_0","EECs_1","EECs_2", "Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Goblet","lowQ","Pericytes", "Reg4_secretory","Stromal","TA")
new.cluster.ids <- c("lowQ","Glial","EECs", "Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Goblet","lowQ","Pericytes", "Reg4_secretory","Stromal","TA")
subCl$annotation1 <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation1", label = TRUE,raster=FALSE)

### subcluster Stromal 
Idents(subCl) <- "annotation1"
subCl <- FindSubCluster(subCl,cluster = "Stromal",graph.name = "SCT_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c("Stromal_0","Stromal_1","Stromal_2","Stromal_3","Stromal_4","Stromal_5"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster", label.size = 5)
DimPlot(sub_celltype ,label = TRUE, group.by = "sub.cluster", split.by = "condition")

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

FeaturePlot(sub_celltype, features = c("Myh11","Acta2","Gli1","Foxl1","Cd34","Cd90","Pdgfra","Rgs5")) 

DotPlot(sub_celltype, features = c("Col1a1","Myh11","Acta2","Gli1","Foxl1","Cd34","Cd90","Pdgfra","Actg2","Cnn1","Des","Rgs5"),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =100, wt = avg_log2FC))

# rename
current.cluster.ids <- c("EECs", "Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Glial", "Goblet","lowQ","Pericytes", "Reg4_secretory",
                         "Stromal_0","Stromal_1","Stromal_2","Stromal_3","Stromal_4","Stromal_5","TA")
new.cluster.ids <- c("EECs", "Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Glial", "Goblet","lowQ","Pericytes", "Reg4_secretory",
                     "Stromal_Pdgfra","SMCs","Stromal_Cd34","Pericytes","SMCs","lowQ","TA")
subCl$annotation1 <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation1", label = TRUE,raster=FALSE)

# remove lowQ 
Idents(subCl) <- "annotation1"
obj <- subset(subCl, idents = c("EECs", "Endothelial","Enterocytes_Alpi_low","Enterocytes_Alpi_high","Glial", "Goblet","Pericytes", "Reg4_secretory",
                                "Stromal_Pdgfra","SMCs","Stromal_Cd34","Pericytes","TA"))
DimPlot(obj, group.by = "annotation1", label = TRUE,raster=FALSE)

## plot marker genes in DotPlot
markers <- c("Epcam","Lgr5","Axin2", "Mki67","Atoh1","Muc2", "Reg4", "Chga","Chgb", 
              "Krt19","Guca2a","Alpi", "S100b", "Rgs5","Pecam1","Myh11","Acta2","Vim","Col1a1","Cd34","Pdgfra","Actg2","Cnn1","Des"
             )
Idents(obj) <- "annotation1"
p <- DotPlot(obj, features = markers, scale = TRUE,cols = c("white","darkred"), dot.scale = 5) + theme(axis.text.x = element_text(angle = 45)) 
ggsave(file.path(annotation_plots_dir,"DotPlot_markers_CD45neg_neo_P14_annotated.svg"), width = 8, height = 5, plot = p)

obj$annotation <- obj$annotation1

### DEGs per annotated cluster 
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

##### save objects
saveRDS(obj, file.path(seurat_objects_dir,"Neo_P14_CD45neg_colon_WT_PHIL_anno.rds"))

