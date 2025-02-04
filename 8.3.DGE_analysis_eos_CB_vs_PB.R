########## This code applies DGE analysis between CB and PB B-eos-liek eosinophils ##########
# Figure 2

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.Packages_and_functions.R")

### ##load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/CB_PB_eos_integrated_anno.rds")

### extract only B-eos-like because there are not Precursors in PB
Idents(obj) <- "annotation"
obj <- subset(obj, idents = "B_Eos_like")

##### DEG analysis 
### generate pseudobulks per sample 
pb <- AggregateExpression(obj, assays = "RNA", return.seurat = T, group.by = c("type","condition"))

### run DEG analysis with DESeq2
Idents(pb) <- "type"
DEG_two_cond_pb_DESeq2(pb, "CB", "PB", "/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CB_PB/")

##### plot GOI in heatmap 
degs <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CB_PB/DESeq2_CB_PB.csv")

goi <- c("IL1RL1","TNFAIP3","IRF1", "NFKBIA","NFKB2","RELB","NFKBID","NFKBIB","MAP3K8", #Immune signaling
                       "THBS1","AREG","PLAUR","ADAM17",#Tissue remodeling
                       "S100A8","CLEC12A","RNASE2", "MBP","CLC", #AMPs and granule proteins
                       "CD69","CD63","SIGLEC8", #activation and secretory activity 
                       "CD58","CD53","HLA-B","CD55" #immune regulation 
)

degs_goi <- degs[degs$X %in% goi,]

# heatmap 
degs_goi_for_plotting <- degs_goi[,c(1,3)]
rownames(degs_goi_for_plotting) <- degs_goi_for_plotting$X
degs_goi_for_plotting$X <- NULL

ComplexHeatmap::Heatmap(degs_goi_for_plotting, name=paste0("CB vs. PB"),
                        column_names_gp = grid::gpar(fontsize = 10),
                        row_names_gp = grid::gpar(fontsize = 10),
                        column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                        column_names_rot = 45,cluster_rows = FALSE, cluster_columns = TRUE,column_title = "genes",
                        col = circlize::colorRamp2(c(-2, 0, 6), c("blue", "white", "red")),
                        heatmap_legend_param = list(
                          color_bar = "continuous"
                        ))

## print significant genes 
sig_genes <- (degs_goi[degs_goi$p_val_adj <= 0.05,])$X
print(sig_genes)

## print non-significant gens 
non_sig_genes <- (degs_goi[degs_goi$p_val_adj > 0.05,])$X
print(non_sig_genes)

