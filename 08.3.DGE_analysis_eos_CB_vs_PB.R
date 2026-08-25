########## This code applies DGE analysis between CB and PB B-eos-liek eosinophils ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir, "01.06.Functions_DEGs.R"))

### ##load objects 
obj <- readRDS(file.path(seurat_objects_dir,"CB_PB_eos_integrated_anno.rds"))

### extract only B-eos-like because there are not Precursors in PB
Idents(obj) <- "annotation"
obj <- subset(obj, idents = "B_Eos_like")

##### DEG analysis 
### generate pseudobulks per sample 
# type = CB/PB, condition = sample id, experiment = Exp1/Exp6
pb <- AggregateExpression(obj, assays = "RNA", return.seurat = TRUE, group.by = c("type","condition"))

counts <- GetAssayData(pb, assay = "RNA", layer = "counts")
meta <- pb@meta.data
meta$type <- factor(meta$type, levels = c("CB","PB"))
meta$experiment <- factor(meta$experiment)

### inspect rank / confounding before fitting
cat("Cross-tabulation of experiment x type:\n")
print(table(meta$type))
# if 3 here, every coefficient is estimatable, if it would be lower than 3 it would mean that f.e. experiment and type are collinear and 
# therefore perfectly confounded, it would not be possible to estimate both
# it is 3 because f.e. PB7 is a PB donor inside Exp1, Experiment and type are not perfectly aligned 
cat("Design matrix rank (should equal number of columns = 3):\n")
print(Matrix::rankMatrix(model.matrix(~ type, data = meta)))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ type)

### filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
# defined CB vs. PB, so positive log2FC = high in CB 
res <- results(dds, contrast = c("type","CB","PB"))
res <- as.data.frame(res)
res$X <- rownames(res)
res$avg_log2FC <- res$log2FoldChange
res$p_val_adj <- res$padj
res$p_val <- res$pvalue
res <- res[,c("X","baseMean","avg_log2FC","p_val_adj","p_val")]
write.csv(res, file.path(DEGs_tables_dir,"DESeq2_CB_vs_PB.csv"))

##### plot GOI in heatmap 
degs <- read.csv(file.path( DEGs_tables_dir,"DESeq2_CB_vs_PB.csv"))

goi <- c("IL1RL1","TNFAIP3","IRF1", "NFKBIA","NFKB2","RELB","NFKBID","NFKBIB","MAP3K8", #Immune signaling
                       "THBS1","AREG","PLAUR","ADAM17",#Tissue remodeling
                       "S100A8","CLEC12A","RNASE2", "MBP","CLC", #AMPs and granule proteins
                       "CD69","CD63","SIGLEC8", #activation and secretory activity 
                       "CD58","CD53","HLA-B","CD55" #immune regulation 
)

degs_goi <- degs[degs$X %in% goi,]

# heatmap 
degs_goi_for_plotting <- degs_goi[,c(2,4)]
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

