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
#neo.query$predicted.celltype <- ifelse(neo.query$prediction.score.max > 0.6, neo.query$predicted.id, "non.defined")

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

##### check annotation: 
#do the label transferred subtypes in the neonatal dataset resemble the eosinophil subtypes from Gurtner et al also based on gene expression? 
# subset neonatal dataset 
Idents(integrated) <- "condition"
neo <- subset(integrated, idents = c("NEO_P14_blood","NEO_P14_bm","NEO_P14_colon","NEO_P14_small_int","NEO_P14_spleen"))
adult <- subset(integrated, idents = c("adult_blood","adult_bm","adult_colon","adult_small_int","adult_spleen"))

Idents(neo) <- "annotation"
DotPlot(neo, features = c("Mki67","Epx","Ear1","Ear2","Prg2",
                               "Alox15","Aldh2","S100a9","S100a6","S100a10","Il5",
                               "Retnla","Ccl9","Il1rl1","Cd24a",
                               "Mmp9","Icosl","Il4","Tgfb1","Pirb","Rara",
                               "Cd80","Ccl3","Il1b","Cxcl2","Cd274","Tnf","Il16"), scale = TRUE,cols = c("white","darkred"), dot.scale = 8) + theme(axis.text.x = element_text(angle = 90)) 

Idents(adult) <- "annotation"
DotPlot(adult, features = c("Mki67","Epx","Ear1","Ear2","Prg2",
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
Idents(obj) <- "annotation"
VlnPlot(obj, features = "immature1",cols = c( "#10A069","#E81818","#E88A1A", "#26DFED","#E8E81A"), pt.size = 0)
VlnPlot(obj, features = "precursors1",cols = c( "#10A069","#E81818","#E88A1A", "#26DFED","#E8E81A"), pt.size = 0)
VlnPlot(obj, features = "A_eos1",cols = c( "#10A069","#E81818","#E88A1A", "#26DFED","#E8E81A"), pt.size = 0)
VlnPlot(obj, features = "B_eos1",cols = c( "#10A069","#E81818","#E88A1A", "#26DFED","#E8E81A"), pt.size = 0)
VlnPlot(obj, features = "Circ_eos1",cols = c( "#10A069","#E81818","#E88A1A", "#26DFED","#E8E81A"), pt.size = 0)
VlnPlot(obj, features = "Cd274",cols = c( "#10A069","#E81818","#E88A1A", "#26DFED","#E8E81A"), pt.size = 0)
VlnPlot(obj, features = "Cd80",cols = c( "#10A069","#E81818","#E88A1A", "#26DFED","#E8E81A"), pt.size = 0)

##### save R objects 
saveRDS(integrated, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds")
integrated <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds")

##### DEGs B-eos neo vs. adult 
current.cluster.ids <- c("adult_blood", "adult_bm","adult_colon","adult_small_int","adult_spleen","adult_spleen_ctrl","adult_stomach",
                         "NEO_P14_blood","NEO_P14_bm","NEO_P14_colon","NEO_P14_small_int","NEO_P14_spleen")
new.cluster.ids <- c("adult", "adult","adult","adult","adult","adult","adult",
                     "neo","neo","neo","neo","neo")
integrated$condition <- plyr::mapvalues(x = integrated$condition, from = current.cluster.ids, to = new.cluster.ids)

Idents(integrated) <- "annotation"
b <- subset(integrated, idents = "basal eosinophils")
Idents(b) <- "condition"
DEG_to_csv_two_cond(b,"neo","adult",FALSE,0.25,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_basal_eos_neo_vs_adult.csv")

Idents(integrated) <- "annotation"
b <- subset(integrated, idents = "circulating eosinophils")
Idents(b) <- "condition"
DEG_to_csv_two_cond(b,"neo","adult",FALSE,0.25,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_circ_eos_neo_vs_adult.csv")


##### MDS plot of label transferred neo vs adult 
### per tissue 
obj <- integrated
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("adult_blood","adult_bm","adult_colon","adult_small_int","adult_spleen_ctrl","NEO_P14_blood","NEO_P14_bm",
                              "NEO_P14_colon","NEO_P14_small_int","NEO_P14_spleen"))
Idents(obj) <- "condition"
av_merged <- AverageExpression(obj, return.seurat = FALSE, assays = "RNA", slot = "counts")
av_merged_df <- as.data.frame(av_merged)

av_merged_mtx <- as.matrix(av_merged_df)

#Create DGEList object
#group is dataset in this case
group <- factor(c(colnames(av_merged_df)))
y <- DGEList(counts=av_merged_mtx,group=group)
y$samples
#library size is number of feature genes in this case 
print(barplot(y$samples$lib.size,names=colnames(y),las=2, main = "Barplot of library sizes"))
#Remove lowly expressed genes
keep <- filterByExpr(y, group=y$samples$group)
y <- y[keep,]
#Normalization
y <- calcNormFactors(y)
y$samples
#we want to compare the groups to each other, which is the eos subtype in our case 
design <- model.matrix(~group)
group <- y$samples$group

# matching colors
group_labels <- levels(group)

group_colors <- c("#1132F4","#F411E9","#11F421",  "#F4ED11","#F26B0F",
                  "#147BC1","#F1C2F2","#086621", "#EDE894","#9B7240")

pdf("/scratch/khandl/eos_datasets_comp/MDS_comp_neo_adult.pdf",  width = 10, height = 7)
plotMDS(y, pch = 16, col = group_colors, main = "MDS Plot",cex = 2)
# Add legend (customize 'group_labels' and 'group_colors' as needed)
#legend("topright", legend = group_labels, col = group_colors, pch = 16,cex = 1)
loc <- plotMDS(y, plot = FALSE)  # extract coordinates without plotting
text(loc$x, loc$y, labels = colnames(y), pos = 4, cex = 0.6)
dev.off() 

### per annotation 
obj <- integrated
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("adult_blood","adult_bm","adult_colon","adult_small_int","adult_spleen_ctrl","NEO_P14_blood","NEO_P14_bm",
                              "NEO_P14_colon","NEO_P14_small_int","NEO_P14_spleen"))

current.cluster.ids <- c("adult_blood","adult_bm","adult_colon","adult_small_int","adult_spleen_ctrl","NEO_P14_blood","NEO_P14_bm",
                         "NEO_P14_colon","NEO_P14_small_int","NEO_P14_spleen")
new.cluster.ids <- c("adult","adult","adult","adult","adult","neo","neo",
                     "neo","neo","neo")
obj$age <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)
obj$anno_cond <- paste0(obj$annotation, "_",obj$age)

Idents(obj) <- "anno_cond"
av_merged <- AverageExpression(obj, return.seurat = FALSE, assays = "RNA", slot = "counts")
av_merged_df <- as.data.frame(av_merged)

av_merged_mtx <- as.matrix(av_merged_df)

#Create DGEList object
#group is dataset in this case
group <- factor(c(colnames(av_merged_df)))
y <- DGEList(counts=av_merged_mtx,group=group)
y$samples
#library size is number of feature genes in this case 
print(barplot(y$samples$lib.size,names=colnames(y),las=2, main = "Barplot of library sizes"))
#Remove lowly expressed genes
keep <- filterByExpr(y, group=y$samples$group)
y <- y[keep,]
#Normalization
y <- calcNormFactors(y)
y$samples
#we want to compare the groups to each other, which is the eos subtype in our case 
design <- model.matrix(~group)
group <- y$samples$group

# matching colors
group_labels <- levels(group)

group_colors <- c("#10A069","#11F421", 
                  "#E8E81A","#EDE894",
                  "#26DFED", "#1132F4",
                  "#E88A1A","#9B7240",
                  "#E81818","#ED8597")

pdf("/scratch/khandl/eos_datasets_comp/MDS_comp_neo_adult_anno.pdf",  width = 10, height = 7)
plotMDS(y, pch = 16, col = group_colors, main = "MDS Plot",cex = 2)
# Add legend (customize 'group_labels' and 'group_colors' as needed)
legend("topright", legend = group_labels, col = group_colors, pch = 16,cex = 1)
loc <- plotMDS(y, plot = FALSE)  # extract coordinates without plotting
text(loc$x, loc$y, labels = colnames(y), pos = 4, cex = 0.6)
dev.off() 


