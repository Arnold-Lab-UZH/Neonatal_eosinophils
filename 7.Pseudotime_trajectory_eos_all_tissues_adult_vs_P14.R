########## This code infers pseudotime trajectory onto the adult and neo P14 eosinophils and analyses DEGs along the pseudotime between ages ##########
# Figure 2 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds")

# remove spleen and stomach 
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("adult_blood","adult_bm","adult_colon","adult_small_int","NEO_P14_blood","NEO_P14_bm", "NEO_P14_colon","NEO_P14_small_int"))

# ad age ident into meta.data 
current.cluster.ids <- c("adult_blood","adult_bm","adult_colon","adult_small_int","NEO_P14_blood","NEO_P14_bm", "NEO_P14_colon","NEO_P14_small_int")
new.cluster.ids <- c("adult","adult","adult","adult","NEO_P14","NEO_P14", "NEO_P14","NEO_P14")
obj$age <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

##### integrate and cluster again with only these conditions 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(object = obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$age)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:5)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:5, reduction.name = "umap")
DimPlot(obj,reduction = "umap",group.by = "age",raster=FALSE)
p <- DimPlot(obj,reduction = "umap",group.by = "annotation",raster=FALSE, cols = c("#10A069","#10A069", "#26DFED","#E88A1A", "#E81818" ))
ggsave(paste0("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/umap_annotated.svg"), width = 8, height = 5, plot = p)

obj <- JoinLayers(obj)

## plot marker genes 
p <- FeaturePlot(obj, features = "Cd80", reduction = "umap", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/Cd80.svg", width = 8, height = 5, plot = p)
p <- FeaturePlot(obj, features = "Cd274", reduction = "umap", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/Cd27.svg", width = 8, height = 5, plot = p)
p <- FeaturePlot(obj, features = "Ear1", reduction = "umap", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/Ear1.svg", width = 8, height = 5, plot = p)
p <- FeaturePlot(obj, features = "Epx", reduction = "umap", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/Epx.svg", width = 8, height = 5, plot = p)

## degs from only the neonatal dataset 
Idents(obj) <- "age"
sub <- subset(obj, idents = "NEO_P14")

current.cluster.ids <- c("basal eosinophils","intestinal eosinophils","immature eosinophils","eosinophil progenitors","circulating eosinophils")
new.cluster.ids <-  c("basal-like eosinophils","intestinal eosinophils","immature eosinophils","eosinophil progenitors","basal-like eosinophils")
sub$annotation <- plyr::mapvalues(x = sub$annotation, from = current.cluster.ids, to = new.cluster.ids)

sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")

top10 <- unique(as.data.frame((markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)))$gene)

Idents(sub) <- "annotation"
p <- DotPlot(sub, features = top10,dot.scale = 10, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/DEGs_neo.svg", width = 20, height = 6, plot = p)

##### Follow the condiment workflow 
# following this vignette: https://hectorrdb.github.io/condimentsPaper/articles/TGFB.html
# convert Seurat object to a SCE object 
sce <- sceasy::convertFormat(obj, from="seurat", to="sce")

# plot UMAP grouped by age 
df <- bind_cols(
  as.data.frame(reducedDims(sce)$UMAP),
  sce$age
) %>%
  sample_frac(1)
p1 <- ggplot(df, aes(x = umap_1, y = umap_2, col = ...3)) +
  geom_point(size = .7) +
  scale_color_brewer(palette = "Dark2") +
  labs(col = "age")
p1

### trajectory inference 
sling <- slingshot(sce, clusterLabels = colData(sce)$annotation, reducedDim = "UMAP"
                   ,omega = TRUE, approx_points = 100)

# plot trajectory on UMAP 
df <- bind_cols(
  as.data.frame(reducedDims(sling)$UMAP),
  sling$slingPseudotime_1,
  sling$annotation,
  sling$age
) %>%
  sample_frac(1)
colnames(df) <- c("umap_1","umap_2","slingPseudotime_1","annotation","age")
curve <- slingCurves(sling)[[1]]
p <- ggplot(df, aes(x = umap_1, y = umap_2, col = slingPseudotime_1)) +
  geom_point(size = .7) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
            col = "black", size = 1.5) +  theme_classic(base_size = 25) 
ggsave(paste0("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/umap_pseudotime1.svg"), width = 8, height = 5, plot = p)

### Differential progression between neo P14 and adult 
## density plot 
p <- ggplot(df, aes(x = slingPseudotime_1)) +
  geom_density(alpha = .8, aes(fill = age), col = "transparent") +
  geom_density(aes(col = age), fill = "transparent",
               guide = FALSE, size = 1.5) +
  labs(x = "Pseudotime", fill = "age") +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("#7F7F7C","#8A181A"))
  )) +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent")
ggsave(paste0("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/density_plot.svg"), width = 8, height = 5, plot = p)

# statistical difference between ages 
# KS = Kolmogorow-Smirnow-Test
progressionTest(SlingshotDataSet(sling), conditions = sling$age, method = "KS") #p-value = 1.22e-13

### Differential expression using TradeSeq 
## fit GAM
set.seed(3)
genes = VariableFeatures(obj)
conditions = factor(sling$age)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 12

sling <- fitGAM(counts = sling, nknots = 5, 
                conditions =conditions, parallel = T, BPPARAM = BPPARAM,genes = genes)

saveRDS(sling, "/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT_pseutodime.rds")
sling <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT_pseutodime.rds")

## differential expresssion between conditions
condRes <- conditionTest(sling, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)

# extract genes with significant p values 
condRes_only_sig <- condRes[!is.na(condRes$waldStat),]
condRes_only_sig <- condRes_only_sig[condRes_only_sig$padj <= 0.05,]
condRes_only_sig$gene <- rownames(condRes_only_sig)
condRes_only_sig <- condRes_only_sig[order(condRes_only_sig$waldStat, decreasing = TRUE),]
write.csv(condRes_only_sig, "/scratch/khandl/Neonatal_eosinophils/data_files/pseudotime/DEGs_neo_P14_vs_adult_over_pseudotime_statistics.csv")

### Visiualise genes of interest along pseudotime 

## plot GOI in line graphs 
goi <- c("Prg2","Epx","Ear1","Ear2","Ear6","S100a9","S100a8","Cd63")
for (i in goi) {
  p <- plotSmoothers(sling, assays(sling)$counts,
                     gene = i,
                     alpha = 1, border = TRUE, curvesCols = c("#7F7F7C","#8A181A")) +
    scale_color_manual(values = c("#7F7F7C","#8A181A")) +
    ggtitle(i)
  ggsave(paste0("/scratch/khandl/Neonatal_eosinophils/figures/pseudotime/",i,".pdf"), width = 10, height = 6, plot = p)
}

## plot all genes in hetmap 
# based on mean smoother
yhatSmooth <- 
  predictSmooth(sling, gene = condRes_only_sig$gene, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth,1, scales::rescale))

# sort rows the same order as the data frame that is ordered by waldStat
yhatSmoothScaled_df <- as.data.frame(yhatSmoothScaled)
yhatSmoothScaled <- within(yhatSmoothScaled_df, rownames(yhatSmoothScaled_df) <- factor(rownames(yhatSmoothScaled_df), 
                                                                                        levels = condRes_only_sig$gene))
yhatSmoothScaled_df2 <- yhatSmoothScaled_df[order(match(rownames(yhatSmoothScaled_df),condRes_only_sig$gene)),]

# split for adult and neo 
yhatSmoothScaled_df2_adult <- yhatSmoothScaled_df2[,1:50]
yhatSmoothScaled_df2_neo <- yhatSmoothScaled_df2[,51:100]

heatSmooth_adult <- pheatmap(yhatSmoothScaled_df2[1:50, 1:50],
                             cluster_cols = FALSE,cluster_rows = FALSE,
                             show_rownames = TRUE, show_colnames = FALSE, main = "adult", legend = FALSE,
                             silent = TRUE
)
heatSmooth_adult

heatSmooth_neo <- 
  pheatmap(yhatSmoothScaled_df2[351:387, 51:100],
           cluster_cols = FALSE, cluster_rows = FALSE,
           show_rownames = TRUE, show_colnames = FALSE, main = "neo",
           legend = FALSE, silent = TRUE 
  )
heatSmooth_neo

## plot all genes in hetmap 
genes_of_interest <- c("Gapdh","Eno1","Enox2","Myc","Atp13a3","Hdac8","Dgkh","Dgat1","Elovl6","Prkg1", #metabolism
                       "S100a6","Slc8a1", #Calcium signaling
                       "Csf2rb","Csf3r", #Cytokine receptors
                       "Col1a1","Col3a1","Mmp8", #ECM remodeling
                       "Prg2","Epx","Ear1","Ear2","Ear6", #Granule proteins
                       "S100a8","S100a9","Lyz1","Lyz2","Acod1", #AMPs
                       "Cd63", #Secretor activit
                       "St3gal5","Pomgnt1", #Sialyltransferases
                       "mt-Cytb","mt-Nd2","mt-Co1","mt-Atp8", #Mitochondrial respoiratory chain 
                       "Vcan","Sdk1"#adhesion
)

yhatSmoothScaled_df3 <- yhatSmoothScaled_df2[rownames(yhatSmoothScaled_df2) %in% genes_of_interest,]

heatSmooth <- pheatmap(yhatSmoothScaled_df3,
                       cluster_cols = FALSE,cluster_rows = TRUE,
                       show_rownames = TRUE, show_colnames = FALSE, main = "adult/neo", legend = TRUE,
                       silent = TRUE
)
heatSmooth
