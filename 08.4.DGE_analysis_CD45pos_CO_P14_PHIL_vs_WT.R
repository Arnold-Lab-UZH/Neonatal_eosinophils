########## This code applies DGE analysis in CO CD45 positive cells neo P14 PHIL vs. WT  ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir, "01.06.Functions_DEGs.R"))

### ##load objects 
obj <- readRDS(file.path(seurat_objects_dir,"Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds"))

#normalize counts in the object  
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(obj) <- "RNA" #RNA is the default asssay for DEG analysis

##### DEG analysis 
cell_types <- c("cDC1","cDC2","Mono_Mac","Mac1","Neutrophils")

for(i in cell_types) {
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"neo_P14_phil","neo_P14_wt",FALSE,0.25,paste0(file.path(DEGs_tables_CD45pos_all_cell_types_dir),"/", i,"_neo_phil_vs_neo_wt.csv") )
}

##### plot genes of interest 
### mac/mono 
genes_of_interest <- c("Itgam","Clec4a1","Clec7a","Clec4n","Fcgr2b","Cd14",#antigen presentation and phagocytosis 
                       "Vegfa","Ccl9","Tlr1","Msr1","Mgl2","Mif","Csf2rb","Csf2rb2","Il1rn","Pilra","Trem1","Il1b", #immune response and inflammation
                       "Thbs1","Mmp14","Tgfbi","Fn1","Mmp19" #extracellular matrix and tissue remodeling 
)

#identify genes not present in the dataframe 
df <- read.csv(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Mono_Mac_neo_phil_vs_neo_wt.csv"))
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
write.csv(df_new,file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Mono_Mac_neo_phil_vs_neo_wt_for_plotting.csv"))

heatmap_logFC_goi(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Mono_Mac_neo_phil_vs_neo_wt_for_plotting.csv"),genes_of_interest,
                  "condition",genes_of_interest,c("antigen_pres","immune_response","ECM"),
                  c(6,12,5),c("#F4062E", "#0F9EED","#36AA0B"),c(antigen_pres="#F4062E",  immune_response= "#0F9EED",  ECM= "#36AA0B"))

### Mature mac 1 
#identify genes not present in the dataframe 
df <- read.csv(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Mac1_neo_phil_vs_neo_wt.csv"))
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
write.csv(df_new,file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Mac1_neo_phil_vs_neo_wt_for_plotting.csv"))

heatmap_logFC_goi(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Mac1_neo_phil_vs_neo_wt_for_plotting.csv"),genes_of_interest,
                  "condition",genes_of_interest,c("antigen_pres","immune_response","ECM"),
                  c(6,12,5),c("#F4062E", "#0F9EED","#36AA0B"),c(antigen_pres="#F4062E",  immune_response= "#0F9EED",  ECM= "#36AA0B"))


### cDC1
genes_of_interest <- c("Itgax","Itgae","Cd24a","Clec4n","Cd44","Adgre5")

#identify genes not present in the dataframe 
df <- read.csv(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"cDC1_neo_phil_vs_neo_wt.csv"))
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
write.csv(df_new,file.path(DEGs_tables_CD45pos_all_cell_types_dir,"cDC1_neo_phil_vs_neo_wt_for_plotting.csv"))

heatmap_logFC_goi(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"cDC1_neo_phil_vs_neo_wt_for_plotting.csv"),genes_of_interest,
                  "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### cDC2
#identify genes not present in the dataframe 
df <- read.csv(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"cDC2_neo_phil_vs_neo_wt.csv"))
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
write.csv(df_new,file.path(DEGs_tables_CD45pos_all_cell_types_dir,"cDC2_neo_phil_vs_neo_wt_for_plotting.csv"))

heatmap_logFC_goi(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"cDC2_neo_phil_vs_neo_wt_for_plotting.csv"),genes_of_interest,
                  "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

### neutrophils
#identify genes not present in the dataframe 
df <- read.csv(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Neutrophils_neo_phil_vs_neo_wt.csv"))
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
write.csv(df_new,file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Neutrophils_neo_phil_vs_neo_wt_for_plotting.csv"))

heatmap_logFC_goi(file.path(DEGs_tables_CD45pos_all_cell_types_dir,"Neutrophils_neo_phil_vs_neo_wt_for_plotting.csv"),genes_of_interest,
                  "condition",genes_of_interest,c("genes"),
                  length(genes_of_interest),c("#F4062E"),c(genes="#F4062E"))

##### DotPlot of Il1b, Il1a, Tgfb and Tnf in CD45 pos 
obj <- readRDS(file.path(seurat_objects_dir,"Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds"))
Idents(obj) <- "condition"
sub <- subset(obj, idents = c("adult_wt","neo_P14_wt"))
sub$anno_cond <- paste0(sub$annotation, "_",sub$condition)

Idents(sub) <- "annotation"
p <- DotPlot(sub, features = c("Il1b","Il1a","Tgfb1","Tnf","Vegfa"),dot.scale = 10, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 



