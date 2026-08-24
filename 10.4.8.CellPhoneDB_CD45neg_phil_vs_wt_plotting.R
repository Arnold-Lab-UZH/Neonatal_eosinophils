########## This code plots interactions between CD45 negative populations phil vs. wt  ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir, "01.10.Functions_CellPhoneDB.R"))

##### interactions stromal CD34 and epithelial wt and phil 
interactins_of_interest <- c("BMP2_BMR1A_ACR2A", "BMP4_BMR1A_ACR2A",
                             "BMP6_BMPR1A_BMPR2","BMP6_ACVR1_BMPR2","BMP6_BMR1A_ACR2A","BMP6_ACVR_1A2A_receptor",
                             "RSPO3_LGR4","RSPO3_LGR5",
                             "WNT5A_FZD5_LRP6","WNT5A_FZD5_LRP5","WNT5A_FZD7_LRP5","WNT5A_FZD7_LRP6","WNT5A_FZD6_LRP5","WNT5A_FZD6_LRP6")

### wt 
pval <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_pvalues_05_20_2025_100446.txt"),check.names = FALSE)
means <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_interaction_scores_05_20_2025_100446.txt"),check.names = FALSE)

# Stromal Cd34 - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_high")]
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
pval <- read.delim(file.path(CellPhoneDB_neo_CD45neg_phil_output_tables_dir,"statistical_analysis_pvalues_05_20_2025_153041.txt"),check.names = FALSE)
means <- read.delim(file.path(CellPhoneDB_neo_CD45neg_phil_output_tables_dir,"statistical_analysis_interaction_scores_05_20_2025_153041.txt"),check.names = FALSE)

# Stromal Cd34 - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Cd34|Enterocytes_Alpi_high")]
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
ggsave(file.path(Signaling_pathways_plots_dir,"stromalCD34_epi_phil_wt_L-R.svg"), width = 12, height = 6, plot = p)

##### interactions stromal aPDGFRA and epithelial wt and phil 
interactins_of_interest <- c("BMP2_BMR1A_ACR2A","BMP2_BMPR1A_BMPR2",
                             "BMP4_BMR1A_ACR2A",
                             "BMP5_BMPR1A_BMPR2","BMP5_ACVR1_BMPR2","BMP5_BMR1A_ACR2A","BMP5_ACVR_1A2A_receptor",
                             "BMP6_BMR1A_ACR2A","BMP6_ACVR_1A2A_receptor",
                             "BMP7_ACVR1_BMPR2","BMP7_BMR1A_ACR2A",
                             "WNT4_FZD6_LRP6","WNT4_FZD6_LRP5","WNT5A_FZD6_LRP5","WNT5A_FZD6_LRP6")

### wt 
pval <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_pvalues_05_20_2025_100446.txt"),check.names = FALSE)
means <- read.delim(file.path(CellPhoneDB_neo_eos_CD45neg_output_tables_dir,"statistical_analysis_interaction_scores_05_20_2025_100446.txt"),check.names = FALSE)

# Stromal Pdgfra - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_high")]
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
pval <- read.delim(file.path(CellPhoneDB_neo_CD45neg_phil_output_tables_dir,"statistical_analysis_pvalues_05_20_2025_153041.txt"),check.names = FALSE)
means <- read.delim(file.path(CellPhoneDB_neo_CD45neg_phil_output_tables_dir,"statistical_analysis_interaction_scores_05_20_2025_153041.txt"),check.names = FALSE)

# Stromal Pdgfra - Enterocytes Alpi high 
pval1 <- pval[pval$interacting_pair %in% interactins_of_interest ,colnames(pval) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_high")]
colnames(pval1) <- c("interacting_pair","pval")
means1 <- means[means$interacting_pair %in% interactins_of_interest ,colnames(means) %in% c("interacting_pair","Stromal_Pdgfra|Enterocytes_Alpi_high")]
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
ggsave(file.path(Signaling_pathways_plots_dir,"stromalPDGFRA_epi_phil_wt_L-R.svg"), width = 12, height = 6, plot = p)
