########## This code plots at significant interactions between eosinophils and stromal/epithelial cells in both directions ##########

##### link to libraries and functions
source("1.1.config.R")
source(file.path(base_dir,"1.3.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "1.2.Packages.R"))
source(file.path(base_dir, "1.10.Functions_CellPhoneDB.R"))

##### interactions eos with stromal cells and epithelial 
interactins_of_interest <- c("IL1A_IL1_receptor","IL1B_IL1_receptor","IL1RN_IL1_receptor","TGFB1_TGFBR3","TGFB1_TGFbeta_receptor2","TGFB1_TGFbeta_receptor1","TNF_TNFRSF1A","TNF_TNFRSF1B",
                             "VEGFA_NRP1","VEGFA_NRP2")

pval <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir),"statistical_analysis_pvalues_05_20_2025_100446.txt",check.names = FALSE)
means <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir),"statistical_analysis_interaction_scores_05_20_2025_100446.txt",check.names = FALSE)

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
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Eosinophils|Enterocytes_Alpi_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Eosinophils|Enterocytes_Alpi_high")]
colnames(means1) <- c("interacting_pair","means")
df3 <- merge(pval1,means1, by = "interacting_pair")
df3$celltype <- "Eos_Enterocytes_Alpi_high"

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
ggsave(file.path(Signaling_pathways_plots_dir,"stromal_eos_L-R.svg"), width = 12, height = 6, plot = p)


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
ggsave(file.path(Signaling_pathways_plots_dir,"stromal_eos_CCL11_L-R.svg"), width = 12, height = 6, plot = p)

