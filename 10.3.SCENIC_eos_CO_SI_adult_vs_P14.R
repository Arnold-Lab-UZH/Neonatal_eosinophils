######### This code rund SCENIC between adult and neo P14 eos from the CO and SI ##########
# Figure S3 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.7.Functions_GSEA_PROGENy_SCENIC.R")

### load seurat objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_LT.rds")

## extract only CO and SI 
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("adult_colon","adult_small_int","NEO_P14_colon","NEO_P14_small_int"))

##### download the motive databases manually and scp to science apps 
#https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/

#scp /Users/handler/Downloads/mm9-500bp-upstream-7species.mc9nr.feather khandl@cluster.s3it.uzh.ch:data/common/SCENIC
#scp /Users/handler/Downloads/mm9-tss-centered-10kb-7species.mc9nr.feather khandl@cluster.s3it.uzh.ch:data/common/SCENIC

########## all clusters combined ##########
setwd("/scratch/khandl/Neonatal_eosinophils/data_files/SCENIC/CO_SI_eos_adult_vs_neoP14/") 

### Initialize settings
exprMat <- as.matrix(obj[["RNA"]]$counts)
cellInfo <- obj@meta.data

colnames(cellInfo)[which(colnames(cellInfo)=="condition")] <- "CellType"
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

org="mgi" # or hgnc, or dmel
dbDir="/home/khandl/data/common/SCENIC" # RcisTarget databases location
myDatasetTitle="SCENIC all eos" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs

#load motive annotation and rename it 
data(list="motifAnnotations_mgi", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

### filter genes 
#Filter by the number of cells in which the gene is detected (minCountsPerGene, by default 6 UMI counts across all samples)
# and by the number of cells in which the gene is detected (minSamples, by default 1% of the cells)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

### calculating correlation
runCorrelation(exprMat_filtered, scenicOptions)

### run GENIE3: infer potential transcription factor targets based on the expression data
exprMat_filtered <- log2(exprMat_filtered+1) 
motifAnnotations_mgi_v8 <- motifAnnotations
runGenie3(exprMat_filtered, scenicOptions)

### run SCENIC to score GRN (regulons) 
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

exprMat_log <- log2(exprMat+1)
dim(exprMat)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log, skipHeatmap = TRUE, skipTsne = TRUE)

### Binarizing the network
#to make it on/off, makes it more clear to compare between conditions 
#Clustering / dim reduction on the regulon activity.
nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them

runSCENIC_4_aucell_binarize(scenicOptions, skipBoxplot = TRUE, skipHeatmaps = TRUE, skipTsne = TRUE)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

scenicOptions <- readRDS(file="int/scenicOptions.Rds") 

### plot regulons in heatmap
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                        column_names_gp = grid::gpar(fontsize = 8),
                        row_names_gp = grid::gpar(fontsize = 8),
                        column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                        column_names_rot = 45)

##### save average regulon activity and scaled average regulon activity in csv file 
write.csv(regulonActivity_byCellType, "/data/khandl/Neonatal_eosinophils/SCENIC/CO_SI_eos_adult_neoP14_regulons_actvity.csv")

##### plot the difference between adult and neo P14 CO in heatmap 
### load regulon activity file 
regulonActivity_byCellType <- read.csv("/data/khandl/Neonatal_eosinophils/SCENIC/CO_SI_eos_adult_neoP14_regulons_actvity.csv")
rownames(regulonActivity_byCellType) <- regulonActivity_byCellType$X
regulonActivity_byCellType$X <- NULL

regulonActivity_byCellType <- regulonActivity_byCellType[,colnames(regulonActivity_byCellType) %in% c("adult_colon","NEO_P14_colon")]

regulonActivity_byCellType$diff <- regulonActivity_byCellType$NEO_P14_colon - regulonActivity_byCellType$adult_colon

### select regulons of interest 
regulons_specific <- c("Cebpb (51g)","Runx3_extended (21g)","Gata2_extended (135g)","Hes1_extended (16g)","Elf2_extended (19g)",
                       "Ikzf2_extended (12g)","Ikzf1_extended (21g)", "Gata1_extended (195g)","Gatad1_extended (208g)",#differentiation and development
                       "Fosl2 (27g)","Junb (74g)","Fosb (45g)","Fos (28g)","Jdp2 (73g)","Elk4_extended (46g)","Etv5_extended (41g)",#8 differentiaito, survival and apoptosis
                       "Stat6_extended (216g)",#1 IL4/13 signaling 
                       "Relb (133g)","Nfkb2 (60g)","Nfkb1 (340g)",#3 NFkB pathway 
                       "Srf (27g)", #1 migration 
                       "Foxo1_extended (11g)","Atf4 (19g)","Bach1_extended (73g)",#3 stress response 
                       "Klf2 (77g)","Nr3c1 (22g)","Crem (41g)","Arnt (63g)","Rest (90g)","Taf1_extended (115g)","Yy1 (152g)","Rarg (44g)",
                       "Elf1_extended (474g)","Ep300_extended (42g)" #Functional regulation
)

regulonActivity_byCellType2 <- regulonActivity_byCellType[rownames( regulonActivity_byCellType) %in% regulons_specific,]
regulonActivity_byCellType2$adult_colon <- NULL
regulonActivity_byCellType2$NEO_P14_colon <- NULL

regulonActivity_byCellType2 <- regulonActivity_byCellType2 %>% arrange(factor(rownames(regulonActivity_byCellType2), 
                                                                              levels = regulons_specific))

print(ComplexHeatmap::Heatmap(regulonActivity_byCellType2, name=paste0("Diff regulon activity"),
                              column_names_gp = grid::gpar(fontsize = 10),
                              row_names_gp = grid::gpar(fontsize = 10),
                              column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                              column_names_rot = 45,cluster_rows = FALSE, cluster_columns = FALSE,column_title = "Regulons"))

##### statistical analysis between adult and neo P14 
scenicOptions <- readRDS(file="/scratch/khandl/Neonatal_eosinophils/data_files/SCENIC/CO_SI_eos_adult_vs_neoP14/int/scenicOptions.Rds") 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
df <- as.data.frame(t(getAUC(regulonAUC)))
df$Cell <- rownames(df)
rownames(df) <- NULL

#check where the colon adult and neo P14 cells are 
df_test <- df[,63:66]

df_adult_colon <- df[1:1291,]
df_neo_colon <- df[c(2172:2520,2783:7888),]
df_adult_colon$Type <- "adult_colon"
df_neo_colon$Type <- "neo_colon"

df <- rbind(df_adult_colon, df_neo_colon)

for (i in regulons_specific) {
  df_specific_regulon <- df[,colnames(df) %in% c("Cell",i,"Type")]
  p.value <- t.test(df_specific_regulon[df_specific_regulon$Type %in% "adult_colon",][[i]], 
                    df_specific_regulon[df_specific_regulon$Type %in% "neo_colon",][[i]], alternative = "two.sided")$p.value
  p.adjust <- p.adjust(p.value, method = "bonferroni", n = length(regulons_specific)) #correct for multiple comparisons 
  print(paste0(i, "; bonferroni_adj_p: ",p.adjust))
}





