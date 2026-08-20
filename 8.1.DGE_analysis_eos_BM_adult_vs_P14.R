########## This code applies DGE analysis between adult and neo P14 BM derived eosinophils ##########

##### link to libraries and functions
source("1.1.config.R")
source(file.path(base_dir,"1.3.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "1.2.Packages.R"))
source(file.path(base_dir, "1.6.Functions_DEGs.R"))

### load objects 
obj <- readRDS(file.path(seurat_objects_dir,"Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds"))

##### plot Egfr in DotPlot 
sub <- subset(obj, idents = c("NEO_P14_colon","adult_colon", "NEO_P14_small_int","adult_small_int", "NEO_P14_blood","adult_blood"))
p <- DotPlot(sub, features = c("Egfr"), scale = TRUE,cols = c("lightblue","darkred"), dot.scale = 15) 
ggsave(file.path(DEGs_plots_dir,"Dotplot_egfr_colon_si_blood.svg"), width = 8, height = 5, plot = p)

##### DEG analysis 
Idents(obj) <- "condition"
DEG_to_csv_two_cond(obj,"NEO_P14_bm","adult_bm",FALSE,0.25,file.path(DEGs_tables_dir,"BM_eos_neoP14_vs_adult.csv"))


##### plot GOIs in volcano plot 
goi <- c("Prg2","Prg3","Epx","Ear6","Ear1","Ear2","S100a10","S100a11","S100a6",
         "Tgfbr1","Ifngr2","Stat5b","Map3k3","Map3k5","Map4","Mapk1",
         "Mapk14","Map2k4","Mapk1ip1l","S100a8","Fgfr1")

volcano_DGE_showing_goi(file.path(DEGs_tables_dir,"BM_eos_neoP14_vs_adult.csv"),
                                0.5,0.05,"neo_bm","adult_bm",goi,
                                c("#EB1D43","#0CD7F2","#7F7C7D"),c(-0.5,0.5) ,c(-8,8),
                       file.path(DEGs_plots_dir,"BM_eos_neoP14_vs_adult_volcano.svg") )
