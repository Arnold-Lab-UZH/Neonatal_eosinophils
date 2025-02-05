######### This code rund PROGENy between adult and neo P14 eos from the CO and SI ##########
# Figure S3

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.7.Functions_GSEA_PROGENy_SCENIC.R")

### load seurat objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_LT.rds")

##### CO
Idents(obj) <- "condition"
sub <- subset(obj, idents = c("adult_colon","NEO_P14_colon"))

PROGENy_two_cond(sub,"Mouse","condition","NEO_P14_colon","adult_colon","colon")

##### SI
sub <- subset(obj, idents = c("adult_small_int","NEO_P14_small_int"))

PROGENy_two_cond(sub,"Mouse","condition","NEO_P14_small_int","adult_small_int","small_int")
