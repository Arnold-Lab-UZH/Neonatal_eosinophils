########## This code compares gene signature scores between conditions ##########

##### link to libraries and functions
source("1.1.config.R")
source(file.path(base_dir,"1.3.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "1.2.Packages.R"))
source(file.path(base_dir, "1.8.Functions_signature_scores.R"))

##### load objects 
eos <- readRDS(file.path(seurat_objects_dir,"Neo_P14_adult_eos_CO_SI_blood_BM_spleen_LT.rds"))
CD45enr <- readRDS(file.path(seurat_objects_dir,"Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds"))

##### FcgR-mediated phagocytosis score - based on MSigDB: KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS
### BM
FcgR_phagocytosis_vln_2_cond_mm(eos, "condition",c("adult_bm", "NEO_P14_bm"),
                                "condition","adult_bm","NEO_P14_bm", 
                                c("#7F7F7C","#8A181A"),c(-0.1,0.7),
                                file.path(Signature_scores_plots_dir,"FcgR_BM.svg"))

### blood
FcgR_phagocytosis_vln_2_cond_mm(eos, "condition",c("adult_blood", "NEO_P14_blood"),
                                "condition","adult_blood","NEO_P14_blood", 
                                c("#7F7F7C","#8A181A"),c(-0.1,0.7),
                                file.path(Signature_scores_plots_dir,"FcgR_blood.svg"))

### CO
FcgR_phagocytosis_vln_2_cond_mm(eos, "condition",c("adult_colon", "NEO_P14_colon"),
                                "condition","adult_colon","NEO_P14_colon", 
                                c("#7F7F7C","#8A181A"),c(-0.1,0.7),
                                file.path(Signature_scores_plots_dir,"FcgR_CO.svg"))

##### Granulogenesis score - Farifax et al 2018 and Gurtner et al 2023
### BM
granulogenesis_vln_2_cond_mm(eos, "condition",c("adult_bm", "NEO_P14_bm"),
                                "condition","adult_bm","NEO_P14_bm", 
                                c("#7F7F7C","#8A181A"),c(-0.2,6),
                                file.path(Signature_scores_plots_dir,"Granulogenesis_BM.svg"))

### blood
granulogenesis_vln_2_cond_mm(eos, "condition",c("adult_blood", "NEO_P14_blood"),
                                "condition","adult_blood","NEO_P14_blood", 
                                c("#7F7F7C","#8A181A"),c(-0.2,4),
                                file.path(Signature_scores_plots_dir,"Granulogenesis_blood.svg"))

##### ROS - Maas et al. 2023
### CO
ROS_vln_2_cond_mm(eos, "condition",c("adult_colon", "NEO_P14_colon"),
                  "condition","adult_colon","NEO_P14_colon", 
                  c("#7F7F7C","#8A181A"),c(-1,2),
                  file.path(Signature_scores_plots_dir,"CO_ROS.svg"))

##### Glycolytic activity score - Mi et al 2013
### Mono/mac P14
Idents(CD45enr) <- "annotation"
sub <- subset(CD45enr, idents = "Mono_Mac")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#ADB7A7","#C0DBB4"),c(-1,2),
                         file.path(Signature_scores_plots_dir,"Mono_Mac_glykolytic.svg"))

### mature mac 1 P14
sub <- subset(CD45enr, idents = "Mac1")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#ADE88E","#54EF0C"),c(-1,2),
                         file.path(Signature_scores_plots_dir,"Mac1_glykolytic.svg"))

### cDC1 P14
sub <- subset(CD45enr, idents = "cDC1")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#EDCAE3","#F281CC"),c(-1,2),
                         file.path(Signature_scores_plots_dir,"cDC1_glykolytic.svg"))

### cDC2 P14
sub <- subset(CD45enr, idents = "cDC2")
Glycolytic_vln_2_cond_mm(sub, "condition",c("neo_P14_phil","neo_P14_wt"),
                         "condition","neo_P14_phil","neo_P14_wt",
                         c("#EA88D3","#F20AB1"),c(-1,2),
                         file.path(Signature_scores_plots_dir,"cDC2_glykolytic.svg"))




