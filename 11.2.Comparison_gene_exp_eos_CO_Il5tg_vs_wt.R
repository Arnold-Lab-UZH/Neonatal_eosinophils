########## This code compares CO eosinophils between Il5-tg and Wt mice ##########
# Figure S5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.4.Functions_DEGs.R")

##### load objects 
wt_all_cells <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds")
il5tg_all_cells <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_LT.rds")

##### subset only CO eos 
Idents(wt_all_cells) <- "annotation"
wt <- subset(wt_all_cells, idents = "Eosinophils")
Idents(wt) <- "condition"
wt <- subset(wt, idents = c("adult_wt","neo_P14_wt"))

Idents(il5tg_all_cells) <- "condition"
il5tg <- subset(il5tg_all_cells, idents =c("adult_colon","NEO_P14_colon")) 
#change identity of condition
current.cluster.ids <- c("adult_colon","NEO_P14_colon")
new.cluster.ids <- c("adult_il5tg","neo_P14_il5tg")
il5tg$condition <- plyr::mapvalues(x = il5tg$condition, from = current.cluster.ids, to = new.cluster.ids)
il5tg <- JoinLayers(il5tg)

##### Label transfer to transfer labels from il5-tg to wt 
### pre-processing and clustering of il5tg
il5tg <- NormalizeData(il5tg)
il5tg <- FindVariableFeatures(il5tg)
il5tg <- ScaleData(il5tg,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
il5tg <- RunPCA(il5tg)
il5tg <- FindNeighbors(il5tg, reduction = "pca", dims = 1:30)
il5tg <- FindClusters(il5tg, resolution = 2, algorithm = 4)
il5tg <- RunUMAP(il5tg, reduction = "pca", dims = 1:30,return.model=TRUE)
DimPlot(il5tg,reduction = "umap",group.by = "annotation",combine = FALSE, label.size = 2, split.by = "condition")

##### Label transfer 
query <- wt
query <- NormalizeData(query)
anchors <- FindTransferAnchors(reference = il5tg, query = query, dims = 1:30,reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = il5tg$annotation, dims = 1:30)
query <- AddMetaData(query, metadata = predictions)
#MapQuery is a combined function for TransferData(), IntegrateEmbeddings() and ProjectUMAP()
query <- MapQuery(anchorset = anchors, reference = il5tg, query = query,
                  refdata = list(celltype = "annotation"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(il5tg, reduction = "umap", group.by = "annotation", label = TRUE, label.size = 3,pt.size = 0.5,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,pt.size = 0.5,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

# plot prediction score per condition - all cells are above 0.5, good mapping
VlnPlot(query, features = "predicted.celltype.score", group.by = "condition")

#add a meta.data column combining annotation and predicted.celltype 
query$annotation <- query$predicted.celltype
#merge objects with umap and ref.umap 
integrated <- merge(il5tg, query)
integrated[["umap"]] <- merge(il5tg[["umap"]], query[["ref.umap"]])
DimPlot(integrated, group.by = "condition")
DimPlot(integrated,reduction = "umap",split.by  = "condition",combine = FALSE, label.size = 2, group.by = "condition")

##### plotting 
Idents(integrated) <- "condition"
sub <- subset(integrated, idents = c("adult_il5tg","adult_wt"))
p <- DimPlot(sub,reduction = "umap",group.by  = "condition", cols = c( "#A39F9F","#870C27"), pt.size = 1)
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/eos_comparison/umap_adult_il5tg_wt_CO.svg", width = 10, height = 6)

Idents(integrated) <- "condition"
sub <- subset(integrated, idents = c("neo_P14_il5tg","neo_P14_wt"))
p <- DimPlot(sub,reduction = "umap",group.by  = "condition", cols = c( "#A39F9F","#240FF2"),pt.size = 1)
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/eos_comparison/umap_adult_il5tg_wt_SI.svg", width = 10, height = 6)

## join layers 
integrated <- JoinLayers(integrated)

##### Compare number of exprssed genes 
average_expression <- AverageExpression(integrated, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "condition")
average_expression_df <- as.data.frame(average_expression)
average_expression_df$gene <- rownames(average_expression_df)

### adult 
## Il5-Tg  
il5tg <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.adult.il5tg","gene")] 
#extract genes with >0 gene expression 
il5tg <- il5tg[il5tg$RNA >0,]

## wt 
wt <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.adult.wt","gene")] 
#extract genes with >0 gene expression 
wt <- wt[wt$RNA >0,]

## venn diagram 
x <- list("il5-tg" = il5tg$gene, "wt" =wt$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

### neo 
## Il5-Tg  
il5tg <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.neo.P14.il5tg","gene")] 
#extract genes with >0 gene expression 
il5tg <- il5tg[il5tg$RNA >0,]

## wt 
wt <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.neo.P14.wt","gene")] 
#extract genes with >0 gene expression 
wt <- wt[wt$RNA >0,]

## venn diagram 
x <- list("il5-tg" = il5tg$gene, "wt" =wt$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

##### comparison of DEG analysis from Il5-tg neo P14 vs. adult from 8.2 and WT neo P14 vs. adult 
### WT 
Idents(wt_all_cells) <- "annotation"
wt_eos <- subset(wt_all_cells, idents = "Eosinophils")
Idents(wt_eos) <- "condition"
DEG_to_csv_two_cond(wt_eos,"neo_P14_wt","adult_wt",FALSE,0.1,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_wt.csv")

## genes of interest from 8.2 
genes_of_interest <- c("Siglecf","Lilra6","Clec12a","Siglece", "Cd300a","Cd300lb", # Inhibitory receptors (ITIM)
                       "Thbs1","Alox5ap","Alox5","Alox15","Pla2g7", #enzymes
                       "Vegfa","Ccl9","Il1rl1","Osm", "Ccl6", "Il5","Csf2ra","Cxcl2","Ccl3",#chemokines/cytokines
                       "Ccr1","Fgfr1","Ccr3","Cxcr4", #receptors 
                       "Adam8", "Wnt11","Mmp9","Mmp25",  #MMPs 
                       "Nfkb1","Nfkb2","Nfkbie","Nfkbid","Nfkbib", #NFkb 
                       "Map2k4","Mapk3", #MapK
                       "Slc16a3","Acod1","Slc27a4","Slc7a11","Dgat1", #metabolic activity
                       "mt-Cytb","mt-Nd2","mt-Co1","mt-Nd5","mt-Nd1",#ATP synthesis 
                       "Arpc5","Arpc1b","Arpc2","Arpc4","Arpc3", #Actin dynamics
                       "Arf4","Rab2a","Rab11fip1","Arf1","Arf5", #Vesicualar traficking 
                       "C3ar1","C5ar1", "Fcgr3" #complement and fc gamma 
)

## plot genes and analyse their significance
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Siglecf",0,0,0,0,1)
gene2 <- c("Clec12a",0,0,0,0,1)
gene3 <- c("Siglece",0,0,0,0,1)
gene4 <- c("Il1rl1",0,0,0,0,1)
gene5 <- c("Osm",0,0,0,0,1)
gene6 <- c("Il5",0,0,0,0,1)
gene7 <- c("Mmp9",0,0,0,0,1)
gene8 <- c("Mapk3",0,0,0,0,1)
gene9 <- c("mt-Nd1",0,0,0,0,1)
gene10 <- c("Arpc2",0,0,0,0,1)
gene11 <- c("C3ar1",0,0,0,0,1)
gene12 <- c("C5ar1",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8,gene9,gene10,gene11,gene12)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_wt_for_plotting.csv")

# neo high and adult low 
heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_wt_for_plotting.csv",
                  genes_of_interest,
                  "condition",genes_of_interest,c("ITIM", "Enzymes", "Chemo_Cyto", "Receptors", "MMPs",
                                                  "Nfkb","MapK","metabolic","ATP","Actin","Vesicular", "Fcg"),
                  c(6,5,9,4,4,5,2,5,5,5,5,3),
                  c("#F4062E",  "#AE06F4", "#0F9EED","#06F406", "#F4D206","#F908A3", "#06F406","#F4D206","#F908A3","#06F406","#777575","#0808F9"),
                  c(ITIM="#F4062E", Enzymes= "#AE06F4", Chemo_Cyto= "#0F9EED", 
                    Receptors="#06F406",MMPs= "#F4D206",Nfkb="#F908A3", MapK="#06F406",metabolic= "#F4D206",ATP="#F908A3",Actin="#06F406",
                    Vesicular="#777575",Fcg="#0808F9"))


### il5tg 
Idents(il5tg_all_cells) <- "condition"
il5tg <- subset(il5tg_all_cells, idents =c("adult_colon","NEO_P14_colon")) 
Idents(il5tg) <- "condition"
DEG_to_csv_two_cond(il5tg,"NEO_P14_colon","adult_colon",FALSE,0.1,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_il5tg.csv")

df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_il5tg.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
rownames(df2) <- df2$X
df2$X <- NULL

write.csv(df2,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_il5tg_for_plotting.csv")

# neo high and adult low 
heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_wt_eos/CO_eos_neo_vs_adult_il5tg_for_plotting.csv",
                  genes_of_interest,
                  "condition",genes_of_interest,c("ITIM", "Enzymes", "Chemo_Cyto", "Receptors", "MMPs",
                                                  "Nfkb","MapK","metabolic","ATP","Actin","Vesicular", "Fcg"),
                  c(6,5,9,4,4,5,2,5,5,5,5,3),
                  c("#F4062E",  "#AE06F4", "#0F9EED","#06F406", "#F4D206","#F908A3", "#06F406","#F4D206","#F908A3","#06F406","#777575","#0808F9"),
                  c(ITIM="#F4062E", Enzymes= "#AE06F4", Chemo_Cyto= "#0F9EED", 
                    Receptors="#06F406",MMPs= "#F4D206",Nfkb="#F908A3", MapK="#06F406",metabolic= "#F4D206",ATP="#F908A3",Actin="#06F406",
                    Vesicular="#777575",Fcg="#0808F9"))



