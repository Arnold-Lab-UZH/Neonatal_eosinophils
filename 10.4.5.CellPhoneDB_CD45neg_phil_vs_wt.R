########## This code generates a csv file for each interacting pair and compares phil and wt  ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir, "01.10.Functions_CellPhoneDB.R"))

##### read output file 
interactions_oi <- read.delim(file.path(CellPhoneDB_neo_CD45neg_phil_output_tables_dir,"statistical_analysis_significant_means_05_20_2025_153041.txt"),check.names = FALSE)
interactions_oi <- colnames(interactions_oi)
interactions_oi <- interactions_oi[15:158]

for(i in interactions_oi) {
  LR_per_2_cond_interacting_pair_oi(
    file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_significant_means_05_20_2025_100446.txt"),
    file.path(CellPhoneDB_neo_CD45neg_phil_output_tables_dir,"statistical_analysis_significant_means_05_20_2025_153041.txt"),
    i,"neo_wt","neo_phil",file.path(CellPhoneDB_CD45neg_phil_vs_wt_output_tables_dir))
}