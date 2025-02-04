########## This code applies DGE analysis between adult and neo P14 BM derived eosinophils ##########
# Figure S2, Figure 4

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.Packages_and_functions.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds")

##### DEG analysis 
Idents(obj) <- "condition"
DEG_to_csv_two_cond(obj,"NEO_P14_bm","adult_bm",FALSE,0.25,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_BM_eos/BM_eos_neoP14_vs_adult.csv")


##### plot GOIs in volcano plot 
goi <- c("Prg2","Prg3","Epx","Ear6","Ear1","Ear2","S100a10","S100a11","S100a6",
         "Tgfbr1","Infgr2","Stat5b","Map3k3","Map3k5","Map4","Mapk1",
         "Mapk14","Map2k4","Mapk1ip1l","S100a8","Fgfr1")
volcano_DGE_showing_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_BM_eos/BM_eos_neoP14_vs_adult.csv",
                                0.5,0.05,"neo_bm","adult_bm",goi,
                                c("#EB1D43","#0CD7F2","#7F7C7D"),c(-0.5,0.5) ,c(-8,8),
                        "/scratch/khandl/Neonatal_eosinophils/figures/DEGs_BM_eos/BM_eos_neoP14_vs_adult_volcano.pdf") 

##### plot Egfr in DotPlot 
sub <- subset(obj, idents = c("NEO_P14_colon","adult_colon", "NEO_P14_small_int","adult_small_int", "NEO_P14_blood","adult_blood"))
p <- DotPlot(sub, features = c("Egfr"), scale = TRUE,cols = c("lightblue","darkred"), dot.scale = 15) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/DEGs_BM_eos/Dotplot_egfr_colon.pdf", width = 8, height = 5, plot = p)

