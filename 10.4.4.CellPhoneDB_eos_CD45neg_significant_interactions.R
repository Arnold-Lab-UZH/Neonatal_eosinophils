########## This code extracts significant interactions between eosinophils and CD45 negative subsets ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir, "01.10.Functions_CellPhoneDB.R"))

##### read output file 
df <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_significant_means_05_20_2025_100446.txt"),check.names = FALSE)
df[is.na(df)] <- 0

df$sum <- rowSums(df[15:183])
df <- df[!df$sum == 0,]
write.csv(df,file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"significant_L_R_eosinophils_CD45neg.csv" ))
