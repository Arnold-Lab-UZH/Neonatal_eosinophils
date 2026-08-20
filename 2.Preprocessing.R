########## This code generates Seurat objects from BD Rhapsody Seven Bridges outputs ##########

##### link to libraries and functions
source("1.1.config.R")
source("1.3.Output_directory_output_folder_structure_generation.R")
source(file.path(base_dir, "1.2.Packages.R"))
source(file.path(base_dir,"1.4.Functions_preprocessing.R"))

##### NEO P14 eosinophils Il5tg 
### colon and small intestine 
NEO_P14_colon_1 <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771063_NEO_P14_1_SampleTag01_mm_colon_Expression_Data.st"), 
  project = "NEO_P14_colon_1", condition = "NEO_P14_colon",3,100)

NEO_P14_colon_2 <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771065_NEO_P14_2_SampleTag05_mm_colon_Expression_Data.st"), 
  project = "NEO_P14_colon_2", condition = "NEO_P14_colon",3,100)

NEO_P14_small_int_1 <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771063_NEO_P14_1_SampleTag02_mm_small_int_Expression_Data.st"), 
  project = "NEO_P14_small_int_1", condition = "NEO_P14_small_int",3,100)

NEO_P14_small_int_2 <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771065_NEO_P14_2_SampleTag06_mm_small_int_Expression_Data.st"), 
  project = "NEO_P14_small_int_2", condition = "NEO_P14_small_int",3,100)

# merge objects
obj <- merge(NEO_P14_colon_1, y = c(NEO_P14_small_int_1,NEO_P14_colon_2,NEO_P14_small_int_2),
             add.cell.ids = c("P14_colon_1","P14_small_int_1","P14_colon_2","P14_small_int_2"))

# add mitochondrial features 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")

# assess and apply upper nFeature cutoff 
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
obj <- subset(obj, subset = nFeature_RNA < 5000)

# join layers
obj <- JoinLayers(obj)

# save object 
saveRDS(obj, file = file.path(seurat_objects_dir,"Neo_P14_eos_colon_SI.rds"))

### blood, bone marrow, spleen and spleen adult as a control 
neo_blood <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771067_P14_blood_SampleTag05_mm_Expression_Data.st"), 
  project = "neo_blood", condition = "NEO_P14_blood",3,100)

neo_bm <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771067_P14_BM_SampleTag07_mm_Expression_Data.st"), 
  project = "neo_bm", condition = "NEO_P14_bm",3,100)

neo_spleen <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771067_P14_spleen_SampleTag06_mm_Expression_Data.st"), 
  project = "neo_spleen", condition = "NEO_P14_spleen",3,100)

adult_spleen <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_il5tg_data, "GSM8771067_Adult_spleen_ctrl_SampleTag08_mm_Expression_Data.st"), 
  project = "adult_spleen", condition = "adult_spleen",3,100)

# merge objects
obj <- merge(neo_blood, y = c(neo_bm,neo_spleen,adult_spleen),
             add.cell.ids = c("neo_blood","neo_bm","neo_spleen","adult_spleen"))

# add mitochondrial features 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")

# assess and apply upper nFeature cutoff 
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
obj <- subset(obj, subset = nFeature_RNA < 6000)

# join layers
obj <- JoinLayers(obj)

# save object 
saveRDS(obj, file = file.path(seurat_objects_dir,"Neo_P14_eos_blood_bone_marrow_spleen.rds"))

##### NEO P14 and Adult CD45 enriched colon WT and PHIL - decontX 
colonCD45enr_adult_wt <- create_seurat_from_condition_DecontX(
  path_to_st_file = file.path(raw_Mm_wt_data, "GSM8771069_CD45enr_adult_colon_wt_SampleTag04_mm_Expression_Data.st"), 
  project = "colonCD45enr_adult_wt", condition = "adult_wt",3,100)

colonCD45enr_neo_wt <- create_seurat_from_condition_DecontX(
  path_to_st_file = file.path(raw_Mm_wt_data, "CGSM8771071_CD45enr_neo_P14_colon_wt_SampleTag09_mm_Expression_Data.st"), 
  project = "colonCD45enr_neo_wt", condition = "neo_P14_wt",3,100)

colonCD45enr_neo_phil <- create_seurat_from_condition_DecontX(
  path_to_st_file = file.path(raw_Mm_wt_data, "GSM8771071_CD45enr_neo_P14_colon_phil_SampleTag10_mm_Expression_Data.st"), 
  project = "colonCD45enr_neo_phil", condition = "neo_P14_phil",3,100)

# Merge objects
obj <- merge(colonCD45enr_adult_wt, y = c(colonCD45enr_neo_wt,colonCD45enr_neo_phil),
             add.cell.ids = c("colonCD45enr_adult_wt","colonCD45enr_neo_wt","colonCD45enr_neo_phil"))

# Add mitochondrial features percentage per cell 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")

# assess and apply upper nFeature cutoff 
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
obj <- subset(obj, subset = nFeature_RNA < 6000)

# join layers 
obj <- JoinLayers(obj)

# save object
saveRDS(obj, file = file.path(seurat_objects_dir,"Adult_and_Neo_P14_CD45enr_colon_WT_PHIL.rds"))

##### NEO P14 CD45 negative colon WT and PHIL
colonCD45neg_wt <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_wt_data, "GSM8771073_CD45neg_neo_P14_colon_wt_SampleTag03_mm_Expression_Data.st"), 
  project = "colonCD45neg_neo_wt", condition = "neo_P14_wt",3,100)

colonCD45neg_phil <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_wt_data, "GSM8771073_CD45neg_neo_P14_colon_phil_SampleTag04_mm_Expression_Data.st"), 
  project = "colonCD45neg_neo_phil", condition = "neo_P14_phil",3,100)

# Merge objects
obj <- merge(colonCD45neg_wt, y = colonCD45neg_phil,
             add.cell.ids = c("colonCD45neg_neo_wt","colonCD45neg_neo_phil"))

# Add mitochondrial features percentage per cell 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")

# assess and apply upper nFeature cutoff 
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
obj <- subset(obj, subset = nFeature_RNA < 6000)

# join layers 
obj <- JoinLayers(obj)

# save object
saveRDS(obj, file = file.path(seurat_objects_dir,"Neo_P14_CD45neg_colon_WT_PHIL.rds"))

##### human derived peripheral and cord blood 
### cord blood (CB) + one peripheral blood control 
CB1 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771084_CB1_hs_SampleTag01_Expression_Data.st"), 
                                                 "CB1",3,200,  "CB1","CB","Exp1","healthy")
CB2 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771085_CB2_hs_SampleTag02_Expression_Data.st"), 
                                                 "CB1",3,200,  "CB2","CB","Exp1","healthy")
CB3 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771086_CB3_hs_SampleTag03_Expression_Data.st"), 
                                                 "CB1",3,200,  "CB3","CB","Exp1","healthy")
CB4 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771087_CB4_hs_SampleTag04_Expression_Data.st"), 
                                                 "CB1",3,200,  "CB4","CB","Exp1","healthy")
PB7 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771083_PB7_hs_SampleTag05_Expression_Data.st"), 
                                                 "CB1",3,200,  "PB7","PB","Exp1","healthy")

# Merge objects
obj <- merge(CB1, y = c(CB2, CB3,CB4,PB7),add.cell.ids = c("cb1", "cb2","cb3","cb4","pb7"))

# Add mitochondrial features percentage per cell 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

# join layers 
obj <- JoinLayers(obj)

# save object 
saveRDS(obj, file = file.path(seurat_objects_dir,"Cord_blood.rds"))

### PB 
PB1 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771077_PB1_hs_SampleTag10_Expression_Data.st"), 
                                                 "PB1",3,200,  "PB1","blood_healthy","Exp6","healthy")
PB2 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771078_PB2_hs_SampleTag11_Expression_Data.st"), 
                                                 "PB2",3,200,  "PB2","blood_healthy","Exp6","healthy")
PB3 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771079_PB3_hs_SampleTag12_Expression_Data.st"), 
                                                 "PB3",3,200,  "PB3","blood_healthy","Exp6","healthy")
PB4 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771080_PB4_hs_SampleTag09_Expression_Data.st"), 
                                                 "PB4",3,200,  "PB4","blood_healthy","Exp6","healthy")
PB5 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771081_PB5_hs_SampleTag05_Expression_Data.st"), 
                                                 "PB5",3,200,  "PB5","blood_healthy","Exp6","healthy")
PB6 <- create_seurat_plus_DecontX_human_biopsies(file.path(raw_Hs_data, "GSM8771082_PB6_hs_SampleTag01_Expression_Data.st"), 
                                                 "PB6",3,200,  "PB6","blood_healthy","Exp6","healthy")

# Merge objects
obj <- merge(PB1, y = c(PB2, PB3,PB4,PB5,PB6),add.cell.ids = c("pb1", "pb2","pb3","pb4","pb5","pb6"))

# Add mitochondrial features percentage per cell 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

# join layers 
obj <- JoinLayers(obj)

# save object 
saveRDS(obj, file = file.path(seurat_objects_dir,"Peripheral_blood.rds"))

##### NEO P16 colon after eos depletion neonatally at days 5, 7, 10 (Cre+), not depleted in Cre- 
colon_iDT_Cre_pos <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_eos_Cre_data, "GSM8771075_P16_colon_iDT_Cre_pos_SampleTag09_mm_Expression_Data.st"), 
  project = "colon_iDT_Cre_pos", condition = "P16_colon_iDT_Cre_pos",3,100)

colon_iDT_Cre_neg <- create_seurat_from_condition(
  path_to_st_file = file.path(raw_Mm_eos_Cre_data, "GSM8771075_P16_colon_iDT_Cre_neg_SampleTag10_mm_Expression_Data.st"), 
  project = "colon_iDT_Cre_neg", condition = "P16_colon_iDT_Cre_neg",3,100)

# merge objects
obj <- merge(colon_iDT_Cre_pos, y = c(colon_iDT_Cre_neg),
             add.cell.ids = c("Cre_pos","Cre_neg"))

# add mitochondrial features 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")

# assess and apply upper nFeature cutoff 
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
obj <- subset(obj, subset = nFeature_RNA < 6000)

# join layers
obj <- JoinLayers(obj)

# save object 
saveRDS(obj, file =file.path(seurat_objects_dir,"Neo_P16_iDT_CREpos_CREneg.rds"))
