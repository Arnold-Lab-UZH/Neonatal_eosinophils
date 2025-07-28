########## This code applies CellPhoneDB to predict ligand-receptor interactions between neonatal eosinophils and CD45neg cell types ##########

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.8.Functions_CellPhoneDB.R")

##### load R object 
CD45neg <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL_anno2.rds")
eos <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds")

Idents(CD45neg) <- "condition"
CD45neg <- subset(CD45neg, idents = "neo_P14_wt")

Idents(eos) <- "condition"
eos <- subset(eos, idents = "neo_P14_wt")
Idents(eos) <- "annotation"
eos <- subset(eos, idents = "Eosinophils")

DefaultAssay(eos) <- "RNA"
DefaultAssay(CD45neg) <- "RNA"

obj <- merge(CD45neg, y= c(eos))
obj <- JoinLayers(obj)

#### input file generation 
## CD45neg + eos 
Input_files_CellPhoneDB_generation_mm(obj, "annotation","neo", "/scratch/khandl/Neonatal_eosinophils/CD45neg/") 

## CD45 neg phil 
CD45neg <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL_anno2.rds")
Idents(CD45neg) <- "condition"
CD45neg <- subset(CD45neg, idents = "neo_P14_phil")
Input_files_CellPhoneDB_generation_mm(CD45neg, "annotation","CD45_neo_phil", "/scratch/khandl/Neonatal_eosinophils/CD45neg/") 

##### run in terminal within a conda environment and python 3
ssh khandl@cluster.s3it.uzh.ch
cd /data/khandl
srun --pty -n 1 -c 2 --time=300:00 --mem=300G bash -l
module load mamba
source activate cpdb
python3
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = '/data/khandl/cpdb/v5.0.0/cellphonedb.zip',
  meta_file_path = '/scratch/khandl/Neonatal_eosinophils/CD45neg/neo_meta.txt',
  counts_file_path = '/scratch/khandl/Neonatal_eosinophils/CD45neg/neo_count.txt',
  counts_data = 'ensembl',
  score_interactions = True,
  threshold = 0.1,
  output_path = '/data/khandl/cpdb/output/neo_eos_CD45neg')

cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = '/data/khandl/cpdb/v5.0.0/cellphonedb.zip',
  meta_file_path = '/scratch/khandl/Neonatal_eosinophils/CD45neg/CD45_neo_phil_meta.txt',
  counts_file_path = '/scratch/khandl/Neonatal_eosinophils/CD45neg/CD45_neo_phil_count.txt',
  counts_data = 'ensembl',
  score_interactions = True,
  threshold = 0.1,
  output_path = '/data/khandl/cpdb/output/neo_CD45neg_phil')

##### extract significant interactions between eosinophils and malignant (mouse) and epithelial (human) tumor cells 
### mouse
df <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_significant_means_05_20_2025_100446.txt",check.names = FALSE)
df[is.na(df)] <- 0

interacting_pairs <- c("EECs|Eosinophils","TA|Eosinophils","Enterocytes_Alpi_low|Eosinophils","Endothelial|Eosinophils","Stromal_Cd34|Eosinophils",
                       "SMCs|Eosinophils","Stromal_Pdgfra|Eosinophils","Reg4_secretory|Eosinophils","Glial|Eosinophils","Goblet|Eosinophils",
                       "Pericytes|Eosinophils","Enterocytes_Alip_high|Eosinophils",
                       "Eosinophils|EECs","Eosinophils|TA","Eosinophils|Enterocytes_Alpi_low","Eosinophils|Endothelial","Eosinophils|Stromal_Cd34",
                       "Eosinophils|SMCs","Eosinophils|Stromal_Pdgfra","Eosinophils|Reg4_secretory","Eosinophils|Glial","Eosinophils|Goblet",
                       "Eosinophils|Pericytes","Eosinophils|Enterocytes_Alip_high")

df2 <- df[,colnames(df) %in% c("interacting_pair","directionality", interacting_pairs)]
df2$sum <- rowSums(df2[3:26])
df2 <- df2[!df2$sum == 0,]
write.csv(df2,"/data/khandl/cpdb/output/neo_eos_CD45neg/significant_L_R_eosinophils_CD45neg.csv" )

df2 <- df
df2$sum <- rowSums(df2[15:183])
df2 <- df2[!df2$sum == 0,]
write.csv(df2,"/data/khandl/cpdb/output/neo_eos_CD45neg/significant_L_R_eosinophils_CD45neg_all_combinations.csv" )

##### compare phil and wt CD45neg interactions 
interactions_oi <- read.delim("/data/khandl/cpdb/output/neo_CD45neg_phil/statistical_analysis_significant_means_05_20_2025_153041.txt",check.names = FALSE)
interactions_oi <- colnames(interactions_oi)
interactions_oi <- interactions_oi[15:158]

for(i in interactions_oi) {
  LR_per_2_cond_interacting_pair_oi(
    "/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_significant_means_05_20_2025_100446.txt",
    "/data/khandl/cpdb/output/neo_CD45neg_phil/statistical_analysis_significant_means_05_20_2025_153041.txt",
    i,"neo_wt","neo_phil","/scratch/khandl/Neonatal_eosinophils/CellPhoneDB_CD45neg_neo_phil_vs_wt/")
}

##### counts of interactions eosinophils with CD45negative

### ligand eosinophils - receptor = other CD45neg 
df <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_significant_means_05_20_2025_100446.txt",check.names = FALSE)
df[is.na(df)] <- 0

eos_EECs <- ligand_receptor_counts(df,"Eosinophils","EECs")
eos_Endothelial <- ligand_receptor_counts(df,"Eosinophils","Endothelial")
eos_Enterocytes_Alip_high <- ligand_receptor_counts(df,"Eosinophils","Enterocytes_Alip_high")
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
b <- c(eos_EECs, eos_Endothelial, eos_Enterocytes_Alip_high,eos_Enterocytes_Alpi_low,eos_Glial, eos_Goblet,eos_Pericytes,
       eos_Reg4_secretory,eos_SMCs,eos_Stromal_Cd34,eos_Stromal_Pdgfra,eos_TA)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("ligand = eos, receptor = other")
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/LRcounts_eos_ligand_other_receptor.svg", width = 12, height = 6, plot = p)

### receptor eosinophils - ligand = other CD45neg 
df <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_significant_means_05_20_2025_100446.txt",check.names = FALSE)
df[is.na(df)] <- 0

EECs_eos <- ligand_receptor_counts(df,"EECs","Eosinophils")
Endothelial_eos <- ligand_receptor_counts(df,"Endothelial","Eosinophils")
Enterocytes_Alip_high_eos <- ligand_receptor_counts(df,"Enterocytes_Alip_high","Eosinophils")
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
b <- c(EECs_eos, Endothelial_eos, Enterocytes_Alip_high_eos,Enterocytes_Alpi_low_eos,Glial_eos, Goblet_eos,Pericytes_eos,
       Reg4_secretory_eos,SMCs_eos,Stromal_Cd34_eos,Stromal_Pdgfra_eos,TA_eos)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("ligand = other, receptor = eos")
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/LRcounts_eos_receptor_other_ligand.svg", width = 12, height = 6, plot = p)

### both
a <- c("EECs|Eos","Endothelial|Eos","Enterocytes_Alpi_high|Eos","Enterocytes_Alpi_low|Eos","Glial|Eos","Goblet|Eos","Pericytes|Eos",
       "Reg4_secretory|Eos","SMCs|Eos","Stromal_Cd34|Eos","Stromal_Pdgfra|Eos","TA|Eos")
b <- c(EECs_eos+eos_EECs, Endothelial_eos+eos_Endothelial, Enterocytes_Alip_high_eos+eos_Enterocytes_Alip_high,Enterocytes_Alpi_low_eos+eos_Enterocytes_Alpi_low,
       Glial_eos+eos_Glial, Goblet_eos+eos_Goblet,Pericytes_eos+eos_Pericytes,Reg4_secretory_eos+eos_Reg4_secretory,SMCs_eos+eos_SMCs,
       Stromal_Cd34_eos+eos_Stromal_Cd34,Stromal_Pdgfra_eos+eos_Stromal_Pdgfra,TA_eos+eos_TA)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("eos/CD45neg ligand-receptor both directionalities")
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/LRcounts_eos_other_both.svg", width = 12, height = 6, plot = p)

### adhesion - adhesion 
df <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_significant_means_05_20_2025_100446.txt",check.names = FALSE)
df[is.na(df)] <- 0

eos_EECs <- adhesion_adhesion_counts(df,"Eosinophils","EECs")
eos_Endothelial <- adhesion_adhesion_counts(df,"Eosinophils","Endothelial")
eos_Enterocytes_Alip_high <- adhesion_adhesion_counts(df,"Eosinophils","Enterocytes_Alip_high")
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
Enterocytes_Alip_high_eos <- adhesion_adhesion_counts(df,"Enterocytes_Alip_high","Eosinophils")
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
b <- c(EECs_eos+eos_EECs, Endothelial_eos+eos_Endothelial, Enterocytes_Alip_high_eos+eos_Enterocytes_Alip_high,Enterocytes_Alpi_low_eos+eos_Enterocytes_Alpi_low,
       Glial_eos+eos_Glial, Goblet_eos+eos_Goblet,Pericytes_eos+eos_Pericytes,Reg4_secretory_eos+eos_Reg4_secretory,SMCs_eos+eos_SMCs,
       Stromal_Cd34_eos+eos_Stromal_Cd34,Stromal_Pdgfra_eos+eos_Stromal_Pdgfra,TA_eos+eos_TA)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","counts")
df$counts <- as.numeric(df$counts)

p <- ggplot(df, aes(x = reorder(interacting_pair,-counts), y = counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  + ggtitle("adhesion - adhesion")
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/LRcounts_adhesion.svg", width = 12, height = 6, plot = p)

##### interactions eos with stromal cells and epithelial 
interactins_of_interest <- c("IL1A_IL1_receptor","IL1B_IL1_receptor","IL1RN_IL1_receptor","TGFB1_TGFBR3","TGFB1_TGFbeta_receptor2","TGFB1_TGFbeta_receptor1","TNF_TNFRSF1A","TNF_TNFRSF1B",
                             "VEGFA_NRP1","VEGFA_NRP2")

pval <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_pvalues_05_20_2025_100446.txt",check.names = FALSE)
means <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_interaction_scores_05_20_2025_100446.txt",check.names = FALSE)

# Eos - Stromal Cd34
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Eosinophils|Stromal_Cd34")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Eosinophils|Stromal_Cd34")]
colnames(means1) <- c("interacting_pair","means")
df1 <- merge(pval1,means1, by = "interacting_pair")
df1$celltype <- "Eos_stromalCd34"

# Eos Stromal Pdgfra
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Eosinophils|Stromal_Pdgfra")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Eosinophils|Stromal_Pdgfra")]
colnames(means1) <- c("interacting_pair","means")
df2 <- merge(pval1,means1, by = "interacting_pair")
df2$celltype <- "Eos_stromalPdgfra"

# Eos - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Eosinophils|Enterocytes_Alip_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Eosinophils|Enterocytes_Alip_high")]
colnames(means1) <- c("interacting_pair","means")
df3 <- merge(pval1,means1, by = "interacting_pair")
df3$celltype <- "Eos_Enterocytes_Alip_high"

# Eos - Enterocytes Alpi low 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Eosinophils|Enterocytes_Alpi_low")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Eosinophils|Enterocytes_Alpi_low")]
colnames(means1) <- c("interacting_pair","means")
df4 <- merge(pval1,means1, by = "interacting_pair")
df4$celltype <- "Eos_Enterocytes_Alpi_low"

# Eos - Goblet 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Eosinophils|Goblet")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Eosinophils|Goblet")]
colnames(means1) <- c("interacting_pair","means")
df5 <- merge(pval1,means1, by = "interacting_pair")
df5$celltype <- "Eos_Goblet"

# Eos - Reg4 secretory 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Eosinophils|Reg4_secretory")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Eosinophils|Reg4_secretory")]
colnames(means1) <- c("interacting_pair","means")
df6 <- merge(pval1,means1, by = "interacting_pair")
df6$celltype <- "Eos_Reg4_secretory"

# Eos - TA 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Eosinophils|TA")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Eosinophils|TA")]
colnames(means1) <- c("interacting_pair","means")
df7 <- merge(pval1,means1, by = "interacting_pair")
df7$celltype <- "Eos_TA"

# merge and plot 
df <- rbind(df1, df2)
df <- rbind(df, df3)
df <- rbind(df, df4)
df <- rbind(df, df5)
df <- rbind(df, df6)
df <- rbind(df, df7)

df <- df %>%
  mutate(pval = ifelse(pval == 0, 0.05, pval))

df$interacting_pair <- factor(df$interacting_pair, levels = interactins_of_interest) 

p <-ggplot(df, aes(x = celltype, y = interacting_pair,size =as.numeric(means), color = pval)) +  coord_flip() +
  geom_point() + scale_size(name = "score", range = c(0,10)) + theme_light() + scale_color_gradientn(colors = c("blue", "yellow", "red"), limits = c(0,0.05)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/stromal_eos_L-R.svg", width = 12, height = 6, plot = p)

##### Stromal and Eos CCL11 - CCR3
interactins_of_interest <- "CCL11_CCR3"

# Stromal Pdgfra - Eos 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Pdgfra|Eosinophils")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Pdgfra|Eosinophils")]
colnames(means1) <- c("interacting_pair","means")
df1<- merge(pval1,means1, by = "interacting_pair")
df1$celltype <- "Stromal_Pdgfra_Eos"

# Stromal Cd34 - Eos 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Cd34|Eosinophils")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Cd34|Eosinophils")]
colnames(means1) <- c("interacting_pair","means")
df2 <- merge(pval1,means1, by = "interacting_pair")
df2$celltype <- "Stromal_Cd34_Eos"

# merge and plot 
df <- rbind(df1, df2)

df <- df %>%
  mutate(pval = ifelse(pval == 0, 0.05, pval))

df$interacting_pair <- factor(df$interacting_pair, levels = interactins_of_interest) 

p <-ggplot(df, aes(x = celltype, y = interacting_pair,size =means, color = pval)) +  coord_flip() +
  geom_point() + scale_size(name = "score", range = c(7, 10)) + theme_light() + scale_color_gradientn(colors = c("blue", "yellow", "red"), limits = c(0,0.05)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/stromal_eos_CCL11_L-R.svg", width = 12, height = 6, plot = p)

##### interactions stromal CD34 and epithelial wt and phil 
interactins_of_interest <- c("BMP2_BMR1A_ACR2A", "BMP4_BMR1A_ACR2A",
                             "BMP6_BMPR1A_BMPR2","BMP6_ACVR1_BMPR2","BMP6_BMR1A_ACR2A","BMP6_ACVR_1A2A_receptor",
                             "RSPO3_LGR4","RSPO3_LGR5",
                             "WNT5A_FZD5_LRP6","WNT5A_FZD5_LRP5","WNT5A_FZD7_LRP5","WNT5A_FZD7_LRP6","WNT5A_FZD6_LRP5","WNT5A_FZD6_LRP6")

### wt 
pval <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_pvalues_05_20_2025_100446.txt",check.names = FALSE)
means <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_interaction_scores_05_20_2025_100446.txt",check.names = FALSE)

# Stromal Cd34 - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alip_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alip_high")]
colnames(means1) <- c("interacting_pair","means")
df1 <- merge(pval1,means1, by = "interacting_pair")
df1$celltype <- "Stromal_Cd34_Enterocytes_Alpi_high_wt"

# Stromal Cd34 - Enterocytes Alpi low
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_low")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_low")]
colnames(means1) <- c("interacting_pair","means")
df2 <- merge(pval1,means1, by = "interacting_pair")
df2$celltype <- "Stromal_Cd34_Enterocytes_Alpi_low_wt"

# Stromal Cd34 - Goblet
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Cd34|Goblet")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Cd34|Goblet")]
colnames(means1) <- c("interacting_pair","means")
df3 <- merge(pval1,means1, by = "interacting_pair")
df3$celltype <- "Stromal_Cd34_Goblet_wt"

# Stromal Cd34 - Reg4 secretory 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Cd34|Reg4_secretory")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Cd34|Reg4_secretory")]
colnames(means1) <- c("interacting_pair","means")
df4 <- merge(pval1,means1, by = "interacting_pair")
df4$celltype <- "Stromal_Cd34_Reg4_secretory_wt"

# Stromal Cd34 TA 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Cd34|TA")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Cd34|TA")]
colnames(means1) <- c("interacting_pair","means")
df5 <- merge(pval1,means1, by = "interacting_pair")
df5$celltype <- "Stromal_Cd34_TA_wt"

# merge and plot 
df <- rbind(df1, df2)
df <- rbind(df, df3)
df <- rbind(df, df4)
df <- rbind(df, df5)

df <- df %>%
  mutate(pval = ifelse(pval == 0, 0.05, pval))

df$interacting_pair <- factor(df$interacting_pair, levels = interactins_of_interest) 
df_wt <- df

### phil
pval <- read.delim("/data/khandl/cpdb/output/neo_CD45neg_phil/statistical_analysis_pvalues_05_20_2025_153041.txt",check.names = FALSE)
means <- read.delim("/data/khandl/cpdb/output/neo_CD45neg_phil/statistical_analysis_interaction_scores_05_20_2025_153041.txt",check.names = FALSE)

# Stromal Cd34 - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alip_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alip_high")]
colnames(means1) <- c("interacting_pair","means")
df1 <- merge(pval1,means1, by = "interacting_pair")
df1$celltype <- "Stromal_Cd34_Enterocytes_Alpi_high_phil"

# Stromal Cd34 - Enterocytes Alpi low
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_low")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_low")]
colnames(means1) <- c("interacting_pair","means")
df2 <- merge(pval1,means1, by = "interacting_pair")
df2$celltype <- "Stromal_Cd34_Enterocytes_Alpi_low_phil"

# Stromal Cd34 - Goblet
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Cd34|Goblet")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Cd34|Goblet")]
colnames(means1) <- c("interacting_pair","means")
df3 <- merge(pval1,means1, by = "interacting_pair")
df3$celltype <- "Stromal_Cd34_Goblet_phil"

# Stromal Cd34 - Reg4 secretory 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Cd34|Reg4_secretory")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Cd34|Reg4_secretory")]
colnames(means1) <- c("interacting_pair","means")
df4 <- merge(pval1,means1, by = "interacting_pair")
df4$celltype <- "Stromal_Cd34_Reg4_secretory_phil"

# Stromal Cd34 TA 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Cd34|TA")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Cd34|TA")]
colnames(means1) <- c("interacting_pair","means")
df5 <- merge(pval1,means1, by = "interacting_pair")
df5$celltype <- "Stromal_Cd34_TA_phil"

# merge and plot 
df <- rbind(df1, df2)
df <- rbind(df, df3)
df <- rbind(df, df4)
df <- rbind(df, df5)

df <- df %>%
  mutate(pval = ifelse(pval == 0, 0.05, pval))

df$interacting_pair <- factor(df$interacting_pair, levels = interactins_of_interest) 
df_phil <- df

## plot wt and phil 
df <- rbind(df_wt, df_phil)

p <-ggplot(df, aes(x = celltype, y = interacting_pair,size =as.numeric(means), color = pval)) +  coord_flip() +
  geom_point() + scale_size(name = "score", range = c(0,10)) + theme_light() + scale_color_gradientn(colors = c("blue", "yellow", "red"), limits = c(0,0.05)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/stromalCD34_epi_phil_wt_L-R.svg", width = 12, height = 6, plot = p)

##### interactions stromal aPDGFRA nd epithelial wt and phil 
interactins_of_interest <- c("BMP2_BMR1A_ACR2A","BMP2_BMPR1A_BMPR2",
                             "BMP4_BMR1A_ACR2A",
                             "BMP5_BMPR1A_BMPR2","BMP5_ACVR1_BMPR2","BMP5_BMR1A_ACR2A","BMP5_ACVR_1A2A_receptor",
                             "BMP6_BMR1A_ACR2A","BMP6_ACVR_1A2A_receptor",
                             "BMP7_ACVR1_BMPR2","BMP7_BMR1A_ACR2A",
                             "WNT4_FZD6_LRP6","WNT4_FZD6_LRP5","WNT5A_FZD6_LRP5","WNT5A_FZD6_LRP6")

### wt 
pval <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_pvalues_05_20_2025_100446.txt",check.names = FALSE)
means <- read.delim("/data/khandl/cpdb/output/neo_eos_CD45neg/statistical_analysis_interaction_scores_05_20_2025_100446.txt",check.names = FALSE)

# Stromal Pdgfra - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alip_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alip_high")]
colnames(means1) <- c("interacting_pair","means")
df6 <- merge(pval1,means1, by = "interacting_pair")
df6$celltype <- "Stromal_Pdgfra_Enterocytes_Alpi_high_wt"

# Stromal Pdgfra - Enterocytes Alpi low
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_low")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_low")]
colnames(means1) <- c("interacting_pair","means")
df7 <- merge(pval1,means1, by = "interacting_pair")
df7$celltype <- "Stromal_Pdgfra_Enterocytes_Alpi_low_wt"

# Stromal Pdgfra - Goblet
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Pdgfra|Goblet")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Pdgfra|Goblet")]
colnames(means1) <- c("interacting_pair","means")
df8 <- merge(pval1,means1, by = "interacting_pair")
df8$celltype <- "Stromal_Pdgfra_Goblet_wt"

# Stromal Pdgfra - Reg4 secretory 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Pdgfra|Reg4_secretory")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Pdgfra|Reg4_secretory")]
colnames(means1) <- c("interacting_pair","means")
df9 <- merge(pval1,means1, by = "interacting_pair")
df9$celltype <- "Stromal_Pdgfra_Reg4_secretory_wt"

# Stromal Pdgfra TA 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Pdgfra|TA")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Pdgfra|TA")]
colnames(means1) <- c("interacting_pair","means")
df10 <- merge(pval1,means1, by = "interacting_pair")
df10$celltype <- "Stromal_Pdgfra_TA_wt"

# merge and plot 
df <- rbind(df6, df7)
df <- rbind(df, df8)
df <- rbind(df, df9)
df <- rbind(df, df10)

df <- df %>%
  mutate(pval = ifelse(pval == 0, 0.05, pval))

df$interacting_pair <- factor(df$interacting_pair, levels = interactins_of_interest) 
df_wt <- df

### phil
pval <- read.delim("/data/khandl/cpdb/output/neo_CD45neg_phil/statistical_analysis_pvalues_05_20_2025_153041.txt",check.names = FALSE)
means <- read.delim("/data/khandl/cpdb/output/neo_CD45neg_phil/statistical_analysis_interaction_scores_05_20_2025_153041.txt",check.names = FALSE)

# Stromal Pdgfra - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alip_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alip_high")]
colnames(means1) <- c("interacting_pair","means")
df6 <- merge(pval1,means1, by = "interacting_pair")
df6$celltype <- "Stromal_Pdgfra_Enterocytes_Alpi_high_phil"

# Stromal Pdgfra - Enterocytes Alpi low
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_low")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_low")]
colnames(means1) <- c("interacting_pair","means")
df7 <- merge(pval1,means1, by = "interacting_pair")
df7$celltype <- "Stromal_Pdgfra_Enterocytes_Alpi_low_phil"

# Stromal Pdgfra - Goblet
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Pdgfra|Goblet")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Pdgfra|Goblet")]
colnames(means1) <- c("interacting_pair","means")
df8 <- merge(pval1,means1, by = "interacting_pair")
df8$celltype <- "Stromal_Pdgfra_Goblet_phil"

# Stromal Pdgfra - Reg4 secretory 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Pdgfra|Reg4_secretory")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Pdgfra|Reg4_secretory")]
colnames(means1) <- c("interacting_pair","means")
df9 <- merge(pval1,means1, by = "interacting_pair")
df9$celltype <- "Stromal_Pdgfra_Reg4_secretory_phil"

# Stromal Pdgfra TA 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair", "Stromal_Pdgfra|TA")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair", "Stromal_Pdgfra|TA")]
colnames(means1) <- c("interacting_pair","means")
df10 <- merge(pval1,means1, by = "interacting_pair")
df10$celltype <- "Stromal_Pdgfra_TA_phil"

# merge and plot 
df <- rbind(df6, df7)
df <- rbind(df, df8)
df <- rbind(df, df9)
df <- rbind(df, df10)

df <- df %>%
  mutate(pval = ifelse(pval == 0, 0.05, pval))

df$interacting_pair <- factor(df$interacting_pair, levels = interactins_of_interest) 
df_phil <- df

## plot wt and phil 
df <- rbind(df_wt, df_phil)

p <-ggplot(df, aes(x = celltype, y = interacting_pair,size =as.numeric(means), color = pval)) +  coord_flip() +
  geom_point() + scale_size(name = "score", range = c(0,10)) + theme_light() + scale_color_gradientn(colors = c("blue", "yellow", "red"), limits = c(0,0.05)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/scratch/khandl/eos_NEO/CellPhoneDB/stromalPDGFRA_epi_phil_wt_L-R.svg", width = 12, height = 6, plot = p)


##### DotPlot of Il1b, Il1a, Tgfb and Tnf in CD45 pos 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds")
Idents(obj) <- "condition"
sub <- subset(obj, idents = c("adult_wt","neo_P14_wt"))
sub$anno_cond <- paste0(sub$annotation, "_",sub$condition)

Idents(sub) <- "annotation"
p <- DotPlot(sub, features = c("Il1b","Il1a","Tgfb1","Tnf","Vegfa"),dot.scale = 10, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

