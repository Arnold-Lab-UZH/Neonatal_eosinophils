########## This code compares gene signature scores between conditions ##########
# Figure 2, 4, 5
##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.6.Functions_signature_scores.R")

##### load objects 
eos <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds")
CD45enr <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds")


##### FcgR-mediated phagocytosis score - based on MSigDB: KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS
### BM
FcgR_phagocytosis_vln_2_cond_mm(eos, "condition",c("adult_bm", "NEO_P14_bm"),
                                "condition","adult_bm","NEO_P14_bm", 
                                c("#7F7F7C","#8A181A"),c(-0.1,0.7),
                                "/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/FcgR_BM.svg")

### blood
FcgR_phagocytosis_vln_2_cond_mm(eos, "condition",c("adult_blood", "NEO_P14_blood"),
                                "condition","adult_blood","NEO_P14_blood", 
                                c("#7F7F7C","#8A181A"),c(-0.1,0.7),
                                "/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/FcgR_blood.svg")

### CO
FcgR_phagocytosis_vln_2_cond_mm(eos, "condition",c("adult_colon", "NEO_P14_colon"),
                                "condition","adult_colon","NEO_P14_colon", 
                                c("#7F7F7C","#8A181A"),c(-0.1,0.7),
                                "/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/FcgR_CO.svg")

##### Granulogenesis score - Farifax et al 2018 and Gurtner et al 2023
### BM
granulogenesis_vln_2_cond_mm(eos, "condition",c("adult_bm", "NEO_P14_bm"),
                                "condition","adult_bm","NEO_P14_bm", 
                                c("#7F7F7C","#8A181A"),c(-0.2,6),
                                "/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/Granulogenesis_BM.svg")

### blood
granulogenesis_vln_2_cond_mm(eos, "condition",c("adult_blood", "NEO_P14_blood"),
                                "condition","adult_blood","NEO_P14_blood", 
                                c("#7F7F7C","#8A181A"),c(-0.2,4),
                                "/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/Granulogenesis_blood.svg")

##### ROS - Maas et al. 2023
### CO
ROS_vln_2_cond_mm(eos, "condition",c("adult_colon", "NEO_P14_colon"),
                  "condition","adult_colon","NEO_P14_colon", 
                  c("#7F7F7C","#8A181A"),c(-1,2),
                  "/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/CO_ROS.svg")

##### Glycolytic activity score - Mi et al 2013
### Mono/mac P14
sub <- subset(CD45enr, idents = "Mono_Mac")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#ADB7A7","#C0DBB4"),c(-1,2),
                         paste0("/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/","Mono_Mac","_glykolytic.svg"))

### mature mac 1 P14
sub <- subset(CD45enr, idents = "Mac1")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#ADE88E","#54EF0C"),c(-1,2),
                         paste0("/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/","Mac1","_glykolytic.svg"))

### cDC1 P14
sub <- subset(CD45enr, idents = "cDC1")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#EDCAE3","#F281CC"),c(-1,2),
                         paste0("/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/","cDC1","_glykolytic.svg"))

### cDC2 P14
sub <- subset(CD45enr, idents = "cDC2")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#EA88D3","#F20AB1"),c(-1,2),
                         paste0("/scratch/khandl/Neonatal_eosinophils/figures/signature_scores/","cDC2","_glykolytic.svg"))




