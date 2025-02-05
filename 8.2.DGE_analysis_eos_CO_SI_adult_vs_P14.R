########## This code applies DGE analysis between adult and neo P14 CO and SI derived eosinophils ##########
# Figure 3, Figure S3 and Figure 4

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.4.Functions_DEGs.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_LT.rds")

##### DEG analysis 
Idents(obj) <- "condition"

### colon CO 
DEG_to_csv_two_cond(obj,"NEO_P14_colon","adult_colon",FALSE,0.25,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CO_SI_eos/CO_eos_neoP14_vs_adult.csv")

### small intestine SI 
DEG_to_csv_two_cond(obj,"NEO_P14_small_int","adult_small_int",FALSE,0.25,"/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CO_SI_eos/SI_eos_neoP14_vs_adult.csv")


##### plot GOIs in volcano plot 
### CO
goi <- c("Cxcl2","Ccnl1","Vegfa","Adora2a","Il1r2","Ccl3","Wnt11","Cd274","Plac8","Thbs1",
                       "Cd300lf","Fth1","Mctp1","C5ar1","Cd300lb","Il5","Cfp","Aldoa","Itgb1","Lilra6",
                       "Alox5","Mmp9","C3ar1","Ccr3","Cd37","Siglece","S100a10","Ccl6","Alox5ap","Siglecf",
                       "Csf1r","Itgb2","S100a11","H2-D1","Adam8") 

volcano_DGE_showing_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CO_SI_eos/CO_eos_neoP14_vs_adult.csv",
                        0.1,0.05,"neo_colon","adult_colon",goi,
                        c("pink","lightblue","gray50"),c(-0.1,0.1) ,c(-5,5),
                        "/scratch/khandl/Neonatal_eosinophils/figures/DEGs_CO_SI_eos/CO_eos_neoP14_vs_adult_volcano.pdf") 

### SI
goi <- c("Rara","Ier2","S100a10","Cd24a","Krt80","Tgfbi","Siglecg","Cd47","Csf1r","H2-D1",
         "Il1a","Adamdec1","Il1b","Il16","Il1r2","Ccrl2","Ccl3","Cd80","Cxcl2","Csf2rb","Dusp6")

volcano_DGE_showing_goi("/scratch/khandl/Neonatal_eosinophils/data_files/DEGs_CO_SI_eos/SI_eos_neoP14_vs_adult.csv",
                        0.1,0.05,"neo_colon","adult_colon",goi,
                        c("pink","lightblue","gray50"),c(-0.1,0.1) ,c(-5,5),
                        "/scratch/khandl/Neonatal_eosinophils/figures/DEGs_CO_SI_eos/SI_eos_neoP14_vs_adult_volcano.pdf") 

##### plot GOIs in heatmap - plot average expression 
### extract only CO and SI
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("adult_colon","adult_small_int","NEO_P14_colon","NEO_P14_small_int"))

### plot
goi <- c("Siglecf","Lilra6","Clec12a","Siglece", "Cd300a","Cd300lb", # Inhibitory receptors (ITIM)
                       "Thbs1","Alox5ap","Alox5","Alox15","Pla2g7", #enzymes
                       "Vegfa","Ccl9","Il1rl1","Osm", "Ccl6", "Il5","Csf2ra","Cxcl2","Ccl3",#chemokines/cytokines
                       "Ccr1","Fgfr1","Ccr3","Cxcr4", #receptors 
                       "C3ar1","C5ar1", "Fcgr3", #complement and fc gamma 
                       "Adam8","Mmp9","Mmp25",  #MMPs 
                       "Nfkb1","Nfkb2","Nfkbie","Nfkbid","Nfkbib", #NFkb 
                       "Map2k4","Mapk3", #MapK
                       "Slc16a3","Acod1","Slc27a4","Slc7a11","Dgat1", #metabolic activity
                       "mt-Cytb","mt-Nd2","mt-Co1","mt-Nd5","mt-Nd1",#ATP synthesis 
                       "Arpc5","Arpc1b","Arpc2","Arpc4","Arpc3", #Actin dynamics
                       "Arf4","Rab2a","Rab11fip1","Arf1","Arf5" #Vesicular trafficking 
         ) 

heatmap_goi_coi(obj, "condition",goi,"RNA",c("ITIM", "Enzymes", "Chemo_Cyto", "Receptors", "MMPs",
                                                     "Nfkb","MapK","metabolic","ATP","Actin","Vesicular", "Fcg"),
                c(6,5,9,4,3,3,5,2,5,5,5,5),
                c("#F4062E",  "#AE06F4", "#0F9EED","#06F406", "#F4D206","#F908A3", "#06F406","#F4D206","#F908A3","#06F406","#777575","#0808F9"),
                c(ITIM="#F4062E", Enzymes= "#AE06F4", Chemo_Cyto= "#0F9EED", 
                  Receptors="#06F406",MMPs= "#F4D206",Nfkb="#F908A3", MapK="#06F406",metabolic= "#F4D206",ATP="#F908A3",Actin="#06F406",
                  Vesicular="#777575",Fcg="#0808F9"),F,F)

##### plot Fcgr related genes in DotPlot 
sub <- subset(obj, idents = c("NEO_P14_colon","adult_colon"))
p <- DotPlot(sub, features = c("Fcgr1","Fcgr3","Fcgr4","Fcrla","Sdhc","Fcgr2b"), scale = FALSE,cols = c("white","darkred"), dot.scale = 20) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/DEGs_CO_SI_eos/Dotplot_fcgr.svg", width = 8, height = 5, plot = p)



