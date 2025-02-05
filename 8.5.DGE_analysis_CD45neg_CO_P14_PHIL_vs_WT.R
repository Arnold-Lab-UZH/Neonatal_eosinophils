########## This code applies DGE analysis in CO CD45 negative cells neo P14 PHIL vs. WT  ##########
# Figure 5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.4.Functions_DEGs.R")

### ##load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL_anno.rds")

#normalize counts in the object  
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(obj) <- "RNA" #RNA is the default asssay for DEG analysis

##### DEG analysis 
cell_types <- c("Fibroblasts","Pericytes","SMCs")

for(i in cell_types) {
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"neo_P14_phil","neo_P14_wt",FALSE,0.25,paste0("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/", i,"_neo_phil_vs_neo_wt.csv") )
}

##### plotting of genes of interest 
### Fibroblasts
genes_of_interest <- c("Krt18","Krt8","Epcam","Cldn4","Cldn7","Dsp","Lgals4","Lgals3","S100a6",
                       "Ceacam1","Lypd8","Spint2")

#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/Fibroblasts_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Krt18",0,0,0,0,1)
gene2 <- c("S100a10",0,0,0,0,1)
gene3 <- c("Sri",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/Fibroblasts_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/Fibroblasts_neo_phil_vs_neo_wt_for_plotting.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Pericytes
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/Pericytes_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Krt18",0,0,0,0,1)
gene2 <- c("Sri",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/Fibroblasts_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/Fibroblasts_neo_phil_vs_neo_wt_for_plotting.csv",genes_of_interest,
                  "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### SMCs
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/SMCs_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Dsp",0,0,0,0,1)
gene2 <- c("S100a10",0,0,0,0,1)
gene3 <- c("Sri",0,0,0,0,1)
gene4 <- c("Ceacam1",0,0,0,0,1)
gene5 <- c("Spint2",0,0,0,0,1)
gene6<- c("Muc13",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/SMCs_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45neg/SMCs_neo_phil_vs_neo_wt_for_plotting.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

