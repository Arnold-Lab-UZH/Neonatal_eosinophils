########## This code compares BM derived precursors between adults and neo P14 ##########
# Figure S2 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.2.Functions_preprocessing.R")
source("~/Projects/Neonatal_eosinophils/1.3.Functions_annotation.R")
source("~/Projects/Neonatal_eosinophils/1.4.Functions_DEGs.R")

##### take BM data from Gurtner et al, and extract precursors from the BM 
adult_bm_all <- create_seurat_from_condition_old_WTA_version(
  path_to_st_file = file.path("/data/khandl/raw_data", "neo_other_organs", "Eos1_ST_SampleTag06_mm_bonemarrow_Expression_Data.st"), 
  project = "adult_bm_all", condition = "adult_bm_all",3,100)
adult_bm_all$percent.mt <- PercentageFeatureSet(adult_bm_all, pattern = "^mt.")
adult <-adult_bm_all

##### read in neo BM object 
neo <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_eos_blood_bone_marrow_spleen_annotated.rds")

#extract BM  
Idents(neo) <- "condition"
neo <- subset(neo, idents = "NEO_P14_bm")

##### pre-processing of adult 
adult <- NormalizeData(adult)
adult <- FindVariableFeatures(adult)
adult <- ScaleData(adult,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
adult <- RunPCA(adult)

##### merge adult and noe 
merged <- merge(neo, y = c(adult), add.cell.ids = c("neo","adult"))

merged <- JoinLayers(merged)

merged <- subset(merged, subset = percent.mt < 30)

##### cluster
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
merged <- RunPCA(merged)
print(ElbowPlot(merged)) 
merged <- FindNeighbors(object = merged, dims = 1:10)
merged <- FindClusters(merged, resolution = 0.5, random.seed = 5, algorithm = 1)
merged <- RunUMAP(merged, dims = 1:10, seed.use = 5)
DimPlot(merged, label = TRUE, group.by = "seurat_clusters", label.size = 10)
DimPlot(merged, label = TRUE, group.by = "condition")

##### fastMNN integration 
obj <- merged
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:10, reduction.name = "umap.mnn",seed.use = 5)
DimPlot(obj,reduction = "umap",raster=FALSE, label = TRUE) 
DimPlot(obj, label = TRUE, group.by = "condition")

obj <- JoinLayers(obj)

##### extract eosinophils and precursors 
### SingleR for broad annotation
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj, assay = "RNA"), ref = mouse.se, labels = mouse.se$label.main)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "mnn.clusters"
lapply(cell.types, function(x) project_annotation_to_umap_fastMNN(x, results, obj))

### extract eos and precursor clusters 
Idents(obj) <- "seurat_clusters"
sub <- subset(obj, idents = c(8,9,13,12,14) )
DimPlot(sub, group.by = "seurat_clusters")

##### reculster 
obj <- sub
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:10, reduction.name = "umap",seed.use = 5)
DimPlot(obj,reduction = "umap",raster=FALSE, label = TRUE) 
DimPlot(obj, label = TRUE, group.by = "age")
FeaturePlot(obj, features = "Siglecf",raster=FALSE,reduction = "umap")
obj <- JoinLayers(obj)

### SingleR for broad annotation
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj, assay = "RNA"), ref = mouse.se, labels = mouse.se$label.main)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, obj))

### DEGs 
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### analyse marker genes from publications 
## markers Kucinski et al Cell Stem cell 2024 (https://www.sciencedirect.com/science/article/pii/S1934590923004319)
FeaturePlot(obj, features = "Hoxb5",raster=FALSE,reduction = "umap")
FeaturePlot(obj, features = c("Procr","Klf1","Pf4","Csf1r","Siglech","Dntt","Vpreb3","Elane","H2-Aa","Prg2","Mcpt8","Cma1"),raster=FALSE,reduction = "umap")

## markers Tusi Nature 2018 (https://pmc.ncbi.nlm.nih.gov/articles/PMC5899604/)
marker <- c("Hba-a2","Hba-a1","Car2","Thd","Tfrc","Bpqm", #Erythroid 
            "Elane","Prtn3","Lcn2","S100a8","Ltf","Gstm1", #Granulocytic neutrophil
            "Epx","Prg2", #Eosinophil
            "Ccr9","Mzb1","Dntt","Flt3","Lsp1","Ccr7","Cd79a","Lef1",#Lymphatic
            "Csf1r","Ly6c2","Ccr2", #monocytic
            "Pf4","Pbx1","Vwf","Itga2b", #Megacaryocytic
            "Hlf","Gcnt2","Procr", #Multipotential progenitor
            "Prss34","Alox5","Mcpt8","Lmo4", #Basophilic or mast cell
            "H2-Aa","Cd74","H2-Eb1","H2-Ab1","Cst3"#Dendritic
)

DotPlot(obj, features = marker,dot.scale = 10, scale = TRUE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

## markers Jorssen et al Immunity 2024
FeaturePlot(obj, features = c("Epx","Prg3","Cd55","Hlf","Cd177","Flt3","Cd34","Ly6a","Gata1"),raster=FALSE,reduction = "umap")

##### extract cluster 4 and 5 
Idents(obj) <- "seurat_clusters"
sub <- subset(obj, idents = c(4,5) )
DimPlot(sub, group.by = "seurat_clusters" )

### rename clusters
current.cluster.ids <- c(4,5)
new.cluster.ids <- c("Granulocytic_progenitor","Multipotntial_progenitor")
sub$annotation <- plyr::mapvalues(x = sub$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
p <- DimPlot(sub, group.by = "annotation", label = TRUE, split.by = "condition")
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/BM_GMP_MMP/umap_split_by_cond.svg", plot = p, width = 10, height = 6)

FeaturePlot(sub, features = c("Hlf","Hoxa9", "Elane","Prtn3"),raster=FALSE,reduction = "umap")

sub <- JoinLayers(sub)

## save R object
saveRDS(sub,"/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_GMP_MMP_BM.rds")

##### DEG analysis neo P14 vs. adult 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_GMP_MMP_BM.rds")

Idents(obj) <- "annotation"
sub <- subset(obj, idents = "Granulocytic_progenitor")
Idents(sub) <- "condition"
DEG_to_csv_two_cond(sub,"NEO_P14_bm","adult_bm_all",FALSE,0.25,paste0("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_MMP_GMP/Granulocytic_progenitor_neo_vs_adult.csv") )

sub <- subset(obj, idents = "Multipotntial_progenitor")
Idents(sub) <- "condition"
DEG_to_csv_two_cond(sub,"NEO_P14_bm","adult_bm_all",FALSE,0.25,paste0("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_MMP_GMP/Multipotntial_progenitor_neo_vs_adult0.25.csv") )

##### plot granule protein genes in a dotplot 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_GMP_MMP_BM.rds")

Idents(obj) <- "annotation"
GMP <- subset(obj, idents = "Granulocytic_progenitor")
MPP <- subset(obj, idents = "Multipotntial_progenitor")

Idents(GMP) <- "condition"
p <- DotPlot(GMP, features = c("Prg2","Prg3",  "Epx", "Ear6", "Ear1", "Ear2"),scale = FALSE,cols = c("lightblue","darkred"), dot.scale = 15)
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/BM_GMP_MMP/GMP_granulogenesis_dotplot.svg", width = 8, height = 5, plot = p)

Idents(MPP) <- "condition"
p <- DotPlot(MPP, features = c("Prg2","Prg3",  "Epx", "Ear6", "Ear1", "Ear2"),scale = FALSE,cols = c("lightblue","darkred"), dot.scale = 15)
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/BM_GMP_MMP/MPP_granulogenesis_dotplot.svg", width = 8, height = 5, plot = p)

