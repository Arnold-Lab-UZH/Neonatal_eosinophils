########## This code plots FACS data along a trajectory P0-Adult #####
##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.9.Functions_FACS_data_plotting.R")

##### generate plots for each measurement 
measurements_all <- c("Eos_freq","CD63_MFI","MHCII_MFI","SigF_MFI","CD11b_MFI","CD80_MFI",
                      "CD107a_MFI","CD9_MFI","CD16_MFI","CD101_MFI","PDL1_MFI","SSCA_MFI",
                      "Active_eos_freq","CD9_freq","CD63_freq","CD80_freq","CD101_freq",
                      "CD16_freq","CD107a_freq")

measurements_all <- c("PDL1_freq")

for (i in measurements_all) {
  FACS_plotting_heatmap_line_graph_and_stats_along_time_points(
    input_path = "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/1.Data/",
    measurement = i,
    output_path = "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/2.Figures/"
  )
}

for (i in measurements_all) {
  FACS_plotting_heatmap_line_graph_and_stats_along_pseudotime(
    input_path = "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/1.Data/",
    measurement = i,
    output_path = "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/2.Figures/"
  )
}

##### Heatmap for Colon 
measurements_mfi <- c("CD63_MFI","MHCII_MFI","SigF_MFI","CD11b_MFI","CD80_MFI",
                      "CD107a_MFI","CD9_MFI","CD16_MFI","CD101_MFI","PDL1_MFI","SSCA_MFI")
measurements_freq <- c("Eos_freq","Active_eos_freq","CD9_freq","CD63_freq","CD80_freq","CD101_freq",
                       "CD16_freq","CD107a_freq","PDL1_freq")

## MFI measurements 
FACS_plotting_heatmap_tissue_of_interest_non_scaled(
  input_path =  "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/1.Data/", 
  tissue_oi = "CO_" ,
  measurements = measurements_mfi,
  output_path ="/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/2.Figures/Colon/",
  id = "MFI")

FACS_plotting_heatmap_tissue_of_interest_scaled(
  input_path =  "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/1.Data/", 
  tissue_oi = "CO_" ,
  measurements = measurements_mfi,
  output_path ="/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/2.Figures/Colon/",
  id = "MFI")

## freq measurements 
FACS_plotting_heatmap_tissue_of_interest_non_scaled(
  input_path =  "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/1.Data/", 
  tissue_oi = "CO_" ,
  measurements = measurements_freq,
  output_path ="/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/2.Figures/Colon/",
  id = "Freq")

FACS_plotting_heatmap_tissue_of_interest_scaled(
  input_path =  "/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/1.Data/", 
  tissue_oi = "CO_" ,
  measurements = measurements_freq,
  output_path ="/Users/handler/Documents/10.FACS/1.FACS_data_neonatal/2.Figures/Colon/",
  id = "Freq")


