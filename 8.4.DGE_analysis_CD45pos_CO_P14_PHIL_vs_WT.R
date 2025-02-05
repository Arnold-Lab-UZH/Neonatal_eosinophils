########## This code applies DGE analysis in CO CD45 positive cells neo P14 PHIL vs. WT  ##########
# Figure 5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.4.Functions_DEGs.R")

### ##load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds")

#normalize counts in the object  
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(obj) <- "RNA" #RNA is the default asssay for DEG analysis

##### DEG analysis 
cell_types <- c("cDC1","cDC2","Mono_Mac","Mac1","Neutrophils")

for(i in cell_types) {
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"neo_P14_phil","neo_P14_wt",FALSE,0.25,paste0("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/", i,"_neo_phil_vs_neo_wt.csv") )
}

##### plot genes of interest 
### mac/mono 
genes_of_interest <- c("Itgam","Clec4a1","Clec7a","Clec4n","Fcgr2b","Cd14",#antigen presentation and phagocytosis 
                       "Vegfa","Ccl9","Tlr1","Msr1","Mgl2","Mif","Csf2rb","Csf2rb2","Il1rn","Pilra","Trem1","Il1b", #immune response and inflammation
                       "Thbs1","Mmp14","Tgfbi","Fn1","Mmp19" #extracellular matrix and tissue remodeling 
)

#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Mono_Mac_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Clec4a1",0,0,0,0,1)
gene2 <- c("Tlr1",0,0,0,0,1)
gene3 <- c("Tgfbi",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Mono_Mac_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Mono_Mac_neo_phil_vs_neo_wt_for_plotting.csv",genes_of_interest,
                  "condition",genes_of_interest,c("antigen_pres","immune_response","ECM"),
                  c(6,12,5),c("#F4062E", "#0F9EED","#36AA0B"),c(antigen_pres="#F4062E",  immune_response= "#0F9EED",  ECM= "#36AA0B"))

### Mature mac 1 
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Mac1_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Cd14",0,0,0,0,1)
gene2 <- c("Trem1",0,0,0,0,1)
gene3 <- c("Il1b",0,0,0,0,1)
gene4 <- c("Fn1",0,0,0,0,1)
gene5 <- c("Mmp19",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4,gene5)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Mac1_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Mac1_neo_phil_vs_neo_wt_for_plotting.csv",genes_of_interest,
                  "condition",genes_of_interest,c("antigen_pres","immune_response","ECM"),
                  c(6,12,5),c("#F4062E", "#0F9EED","#36AA0B"),c(antigen_pres="#F4062E",  immune_response= "#0F9EED",  ECM= "#36AA0B"))


### cDC1
genes_of_interest <- c("Itgax","Itgae","Cd24a","Clec4n","Cd44","Adgre5")

#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/cDC1_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Cd24a",0,0,0,0,1)
gene2 <- c("Clec4n",0,0,0,0,1)
gene3 <- c("Icosl",0,0,0,0,1)
gene4 <- c("Adgre5",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3,gene4)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/cDC1_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/cDC1_neo_phil_vs_neo_wt_for_plotting.csv",genes_of_interest,
                  "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### cDC2
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/cDC2_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Icosl",0,0,0,0,1)
gene2 <- c("Cd44",0,0,0,0,1)
gene3 <- c("Adgre5",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/cDC2_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/cDC2_neo_phil_vs_neo_wt_for_plotting.csv",genes_of_interest,
                  "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### neutrophils
#identify genes not present in the dataframe 
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Neutrophils_neo_phil_vs_neo_wt.csv")
df2 <- df[df$X %in% genes_of_interest,]
missing_genes <- genes_of_interest[!genes_of_interest %in% df2$X]
print(missing_genes)
#add rows with genes that are missing with 0s 
gene1 <- c("Itgae",0,0,0,0,1)
gene2 <- c("Cd24a",0,0,0,0,1)
gene3 <- c("Clec4n",0,0,0,0,1)

df_new <- rbind(df2,gene1,gene2,gene3)
rownames(df_new) <- df_new$X
df_new$X <- NULL
write.csv(df_new,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Neutrophils_neo_phil_vs_neo_wt_for_plotting.csv")

heatmap_logFC_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CD45pos/Neutrophils_neo_phil_vs_neo_wt_for_plotting.csv",genes_of_interest,
                  "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))


