########## This code counts interactions between eosinophils and CD45 negative cell types in both directions  ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir, "01.10.Functions_CellPhoneDB.R"))

##### ligand eosinophils - receptor = other CD45neg 
df <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_significant_means_05_20_2025_100446.txt"),check.names = FALSE)
df[is.na(df)] <- 0

eos_EECs <- ligand_receptor_counts(df,"Eosinophils","EECs")
eos_Endothelial <- ligand_receptor_counts(df,"Eosinophils","Endothelial")
eos_Enterocytes_Alpi_high <- ligand_receptor_counts(df,"Eosinophils","Enterocytes_Alpi_high")
eos_Enterocytes_Alpi_low <- ligand_receptor_counts(df,"Eosinophils","Enterocytes_Alpi_low")
eos_Glial <- ligand_receptor_counts(df,"Eosinophils","Glial")
eos_Goblet <- ligand_receptor_counts(df,"Eosinophils","Goblet")
eos_Pericytes <- ligand_receptor_counts(df,"Eosinophils","Pericytes")
eos_Reg4_secretory <- ligand_receptor_counts(df,"Eosinophils","Reg4_secretory")
eos_SMCs <- ligand_receptor_counts(df,"Eosinophils","SMCs")
eos_Stromal_Cd34 <- ligand_receptor_counts(df,"Eosinophils","Stromal_Cd34")
eos_Stromal_Pdgfra <- ligand_receptor_counts(df,"Eosinophils","Stromal_Pdgfra")
eos_TA <- ligand_receptor_counts(df,"Eosinophils","TA")

a <- c("Eos|EECs","Eos|Endothelial","Eos|Enterocytes_Alpi_high","Eos|Enterocytes_Alpi_low","Eos|Glial","Eos|Goblet","Eos|Pericytes",
       "Eos|Reg4_secretory","Eos|SMCs","Eos|Stromal_Cd34","Eos|Stromal_Pdgfra","Eos|TA")
b <- c(eos_EECs, eos_Endothelial, eos_Enterocytes_Alpi_high,eos_Enterocytes_Alpi_low,eos_Glial, eos_Goblet,eos_Pericytes,
       eos_Reg4_secretory,eos_SMCs,eos_Stromal_Cd34,eos_Stromal_Pdgfra,eos_TA)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("ligand = eos, receptor = other")
ggsave(file.path(Signaling_pathways_plots_dir,"LRcounts_eos_ligand_other_receptor.svg"), width = 12, height = 6, plot = p)

##### receptor eosinophils - ligand = other CD45neg 
df <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_significant_means_05_20_2025_100446.txt"),check.names = FALSE)
df[is.na(df)] <- 0

EECs_eos <- ligand_receptor_counts(df,"EECs","Eosinophils")
Endothelial_eos <- ligand_receptor_counts(df,"Endothelial","Eosinophils")
Enterocytes_Alpi_high_eos <- ligand_receptor_counts(df,"Enterocytes_Alpi_high","Eosinophils")
Enterocytes_Alpi_low_eos <- ligand_receptor_counts(df,"Enterocytes_Alpi_low","Eosinophils")
Glial_eos <- ligand_receptor_counts(df,"Glial","Eosinophils")
Goblet_eos <- ligand_receptor_counts(df,"Goblet","Eosinophils")
Pericytes_eos <- ligand_receptor_counts(df,"Pericytes","Eosinophils")
Reg4_secretory_eos <- ligand_receptor_counts(df,"Reg4_secretory","Eosinophils")
SMCs_eos <- ligand_receptor_counts(df,"SMCs","Eosinophils")
Stromal_Cd34_eos <- ligand_receptor_counts(df,"Stromal_Cd34","Eosinophils")
Stromal_Pdgfra_eos <- ligand_receptor_counts(df,"Stromal_Pdgfra","Eosinophils")
TA_eos <- ligand_receptor_counts(df,"TA","Eosinophils")

a <- c("EECs|Eos","Endothelial|Eos","Enterocytes_Alpi_high|Eos","Enterocytes_Alpi_low|Eos","Glial|Eos","Goblet|Eos","Pericytes|Eos",
       "Reg4_secretory|Eos","SMCs|Eos","Stromal_Cd34|Eos","Stromal_Pdgfra|Eos","TA|Eos")
b <- c(EECs_eos, Endothelial_eos, Enterocytes_Alpi_high_eos,Enterocytes_Alpi_low_eos,Glial_eos, Goblet_eos,Pericytes_eos,
       Reg4_secretory_eos,SMCs_eos,Stromal_Cd34_eos,Stromal_Pdgfra_eos,TA_eos)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("ligand = other, receptor = eos")
ggsave(file.path(Signaling_pathways_plots_dir,"LRcounts_eos_receptor_other_ligand.svg"), width = 12, height = 6, plot = p)

##### both
a <- c("EECs|Eos","Endothelial|Eos","Enterocytes_Alpi_high|Eos","Enterocytes_Alpi_low|Eos","Glial|Eos","Goblet|Eos","Pericytes|Eos",
       "Reg4_secretory|Eos","SMCs|Eos","Stromal_Cd34|Eos","Stromal_Pdgfra|Eos","TA|Eos")
b <- c(EECs_eos+eos_EECs, Endothelial_eos+eos_Endothelial, Enterocytes_Alpi_high_eos+eos_Enterocytes_Alpi_high,Enterocytes_Alpi_low_eos+eos_Enterocytes_Alpi_low,
       Glial_eos+eos_Glial, Goblet_eos+eos_Goblet,Pericytes_eos+eos_Pericytes,Reg4_secretory_eos+eos_Reg4_secretory,SMCs_eos+eos_SMCs,
       Stromal_Cd34_eos+eos_Stromal_Cd34,Stromal_Pdgfra_eos+eos_Stromal_Pdgfra,TA_eos+eos_TA)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("eos/CD45neg ligand-receptor both directionalities")
ggsave(file.path(Signaling_pathways_plots_dir,"LRcounts_eos_other_both.svg"), width = 12, height = 6, plot = p)

##### adhesion - adhesion 
df <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_significant_means_05_20_2025_100446.txt"),check.names = FALSE)
df[is.na(df)] <- 0

eos_EECs <- adhesion_adhesion_counts(df,"Eosinophils","EECs")
eos_Endothelial <- adhesion_adhesion_counts(df,"Eosinophils","Endothelial")
eos_Enterocytes_Alpi_high <- adhesion_adhesion_counts(df,"Eosinophils","Enterocytes_Alpi_high")
eos_Enterocytes_Alpi_low <- adhesion_adhesion_counts(df,"Eosinophils","Enterocytes_Alpi_low")
eos_Glial <- adhesion_adhesion_counts(df,"Eosinophils","Glial")
eos_Goblet <- adhesion_adhesion_counts(df,"Eosinophils","Goblet")
eos_Pericytes <- adhesion_adhesion_counts(df,"Eosinophils","Pericytes")
eos_Reg4_secretory <- adhesion_adhesion_counts(df,"Eosinophils","Reg4_secretory")
eos_SMCs <- adhesion_adhesion_counts(df,"Eosinophils","SMCs")
eos_Stromal_Cd34 <- adhesion_adhesion_counts(df,"Eosinophils","Stromal_Cd34")
eos_Stromal_Pdgfra <- adhesion_adhesion_counts(df,"Eosinophils","Stromal_Pdgfra")
eos_TA <- adhesion_adhesion_counts(df,"Eosinophils","TA")

EECs_eos <- adhesion_adhesion_counts(df,"EECs","Eosinophils")
Endothelial_eos <- adhesion_adhesion_counts(df,"Endothelial","Eosinophils")
Enterocytes_Alpi_high_eos <- adhesion_adhesion_counts(df,"Enterocytes_Alpi_high","Eosinophils")
Enterocytes_Alpi_low_eos <- adhesion_adhesion_counts(df,"Enterocytes_Alpi_low","Eosinophils")
Glial_eos <- adhesion_adhesion_counts(df,"Glial","Eosinophils")
Goblet_eos <- adhesion_adhesion_counts(df,"Goblet","Eosinophils")
Pericytes_eos <- adhesion_adhesion_counts(df,"Pericytes","Eosinophils")
Reg4_secretory_eos <- adhesion_adhesion_counts(df,"Reg4_secretory","Eosinophils")
SMCs_eos <- adhesion_adhesion_counts(df,"SMCs","Eosinophils")
Stromal_Cd34_eos <- adhesion_adhesion_counts(df,"Stromal_Cd34","Eosinophils")
Stromal_Pdgfra_eos <- adhesion_adhesion_counts(df,"Stromal_Pdgfra","Eosinophils")
TA_eos <- adhesion_adhesion_counts(df,"TA","Eosinophils")

a <- c("Eos|EECs","Eos|Endothelial","Eos|Enterocytes_Alpi_high","Eos|Enterocytes_Alpi_low","Eos|Glial","Eos|Goblet","Eos|Pericytes",
       "Eos|Reg4_secretory","Eos|SMCs","Eos|Stromal_Cd34","Eos|Stromal_Pdgfra","Eos|TA")
b <- c(EECs_eos+eos_EECs, Endothelial_eos+eos_Endothelial, Enterocytes_Alpi_high_eos+eos_Enterocytes_Alpi_high,Enterocytes_Alpi_low_eos+eos_Enterocytes_Alpi_low,
       Glial_eos+eos_Glial, Goblet_eos+eos_Goblet,Pericytes_eos+eos_Pericytes,Reg4_secretory_eos+eos_Reg4_secretory,SMCs_eos+eos_SMCs,
       Stromal_Cd34_eos+eos_Stromal_Cd34,Stromal_Pdgfra_eos+eos_Stromal_Pdgfra,TA_eos+eos_TA)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + ggtitle("adhesion - adhesion")
ggsave(file.path(Signaling_pathways_plots_dir,"LRcounts_adhesion.svg"), width = 12, height = 6, plot = p)
