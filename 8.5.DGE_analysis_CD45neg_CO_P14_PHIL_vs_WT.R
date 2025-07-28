########## This code applies DGE analysis in CO CD45 negative cells neo P14 PHIL vs. WT  ##########
# Figure 5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.4.Functions_DEGs.R")

### ##load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL_anno2.rds")

#normalize counts in the object  
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(obj) <- "RNA" #RNA is the default asssay for DEG analysis

##### DEG analysis 
cell_types <- (as.data.frame(table(obj$annotation)))$Var1

for(i in cell_types) {
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"neo_P14_phil","neo_P14_wt",FALSE,0.25,paste0("/scratch/khandl/Neonatal_eosinophils/CD45neg/", i,"_neo_phil_vs_neo_wt.csv") )
}

##### plotting of genes of interest 
genes_of_interest <- c("S100g","S100a16","S100a6","Cdh17","Krt7","Cldn4","Cldn7", "Angptl2","Spint2","Lgals4")

### TA
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/TA_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100a6",0,0,0,0,1)
gene2 <- c("Krt7",0,0,0,0,1)
gene3 <- c("Cldn4",0,0,0,0,1)
gene4 <- c("Cldn7",0,0,0,0,1)
gene5 <- c("Angptl2",0,0,0,0,1)
gene6 <- c("Spint2",0,0,0,0,1)
gene7 <- c("Lgals4",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/TA_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/TA_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Enterocytes Alpi hgih 
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Enterocytes_Alip_high_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100a16",0,0,0,0,1)
gene2 <- c("S100a6",0,0,0,0,1)
gene3 <- c("Cdh17",0,0,0,0,1)
gene4 <- c("Cldn4",0,0,0,0,1)
gene5 <- c("Cldn7",0,0,0,0,1)
gene6 <- c("Angptl2",0,0,0,0,1)
gene7 <- c("Spint2",0,0,0,0,1)
gene8 <- c("Lgals4",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Enterocytes_Alip_high_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Enterocytes_Alip_high_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Enterocytes Alpi low 
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Enterocytes_Alpi_low_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100a16",0,0,0,0,1)
gene2 <- c("S100a6",0,0,0,0,1)
gene3 <- c("Cdh17",0,0,0,0,1)
gene4 <- c("Cldn4",0,0,0,0,1)
gene5 <- c("Cldn7",0,0,0,0,1)
gene6 <- c("Angptl2",0,0,0,0,1)
gene7 <- c("Spint2",0,0,0,0,1)
gene8 <- c("Lgals4",0,0,0,0,1)
gene9 <- c("Krt7",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8,gene9)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Enterocytes_Alpi_low_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Enterocytes_Alpi_low_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Goblet
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Goblet_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("S100a16",0,0,0,0,1)
gene3 <- c("S100a6",0,0,0,0,1)
gene4 <- c("Krt7",0,0,0,0,1)
gene5 <- c("Cldn4",0,0,0,0,1)
gene6 <- c("Cldn7",0,0,0,0,1)
gene7 <- c("Angptl2",0,0,0,0,1)
gene8 <- c("Spint2",0,0,0,0,1)
gene9 <- c("Lgals4",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8,gene9)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Goblet_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Goblet_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Reg4_secretory
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Reg4_secretory_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100a16",0,0,0,0,1)
gene2 <- c("S100a6",0,0,0,0,1)
gene3 <- c("Krt7",0,0,0,0,1)
gene4 <- c("Cldn7",0,0,0,0,1)
gene5 <- c("Angptl2",0,0,0,0,1)
gene6 <- c("Spint2",0,0,0,0,1)
gene7 <- c("Lgals4",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Reg4_secretory_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Reg4_secretory_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### EECs
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/EECs_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("S100a16",0,0,0,0,1)
gene3 <- c("S100a6",0,0,0,0,1)
gene4 <- c("Krt7",0,0,0,0,1)
gene5 <- c("Cldn4",0,0,0,0,1)
gene6 <- c("Cldn7",0,0,0,0,1)
gene7 <- c("Angptl2",0,0,0,0,1)
gene8 <- c("Spint2",0,0,0,0,1)
gene9 <- c("Lgals4",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8,gene9)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/EECs_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/EECs_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Glial 
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Glial_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("S100a16",0,0,0,0,1)
gene3 <- c("Cdh17",0,0,0,0,1)
gene4 <- c("Krt7",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Glial_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Glial_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Endothelial
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Endothelial_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("S100a16",0,0,0,0,1)
gene3 <- c("Cdh17",0,0,0,0,1)
gene4 <- c("Krt7",0,0,0,0,1)
gene5 <- c("Cldn7",0,0,0,0,1)
gene6 <- c("Angptl2",0,0,0,0,1)
gene7 <- c("Spint2",0,0,0,0,1)
gene8 <- c("Lgals4",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Endothelial_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Endothelial_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Pericytes
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Pericytes_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("S100a16",0,0,0,0,1)
gene3 <- c("Cdh17",0,0,0,0,1)
gene4 <- c("Krt7",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Pericytes_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Pericytes_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### SMCs
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/SMCs_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("Cdh17",0,0,0,0,1)
gene3 <- c("Krt7",0,0,0,0,1)
gene4 <- c("Angptl2",0,0,0,0,1)
gene5 <- c("Spint2",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/SMCs_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/SMCs_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Stromal_Cd34
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Stromal_Cd34_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("S100a16",0,0,0,0,1)
gene3 <- c("S100a6",0,0,0,0,1)
gene4 <- c("Cdh17",0,0,0,0,1)
gene5 <- c("Krt7",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Stromal_Cd34_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Stromal_Cd34_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### Stromal_Pdgfra
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/Stromal_Pdgfra_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("S100g",0,0,0,0,1)
gene2 <- c("S100a16",0,0,0,0,1)
gene3 <- c("Cdh17",0,0,0,0,1)
gene4 <- c("Krt7",0,0,0,0,1)
gene5 <- c("Cldn7",0,0,0,0,1)
gene6 <- c("Angptl2",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5,gene6)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/CD45neg/Stromal_Pdgfra_neo_phil_vs_neo_wt.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/CD45neg/Stromal_Pdgfra_neo_phil_vs_neo_wt.csv"
                  ,genes_of_interest, "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

###### Ccl11 expression 
Idents(obj) <- "annotation1"
p <- DotPlot(obj, features = "Ccl11", scale = TRUE,cols = c("white","darkred"), dot.scale = 5) + theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/Neonatal_eosinophils/CD45neg/Ccl11.svg", width = 8, height = 5, plot = p)
