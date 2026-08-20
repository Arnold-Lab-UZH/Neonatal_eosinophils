########## This code plots FACS data along a trajectory P0-Adult #####

##### link to libraries and functions
source("1.1.config.R")
source(file.path(base_dir,"1.3.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "1.2.Packages.R"))
source(file.path(base_dir, "1.11.Functions_FACS_data_plotting.R"))

##### generate plots for each measurement 
measurements_all <- c("Eos_freq","CD63_MFI","MHCII_MFI","SigF_MFI","CD11b_MFI","CD80_MFI",
                      "CD107a_MFI","CD9_MFI","CD16_MFI","CD101_MFI","PDL1_MFI","SSCA_MFI",
                      "Active_eos_freq","CD9_freq","CD63_freq","CD80_freq","CD101_freq",
                      "CD16_freq","CD107a_freq")

for (i in measurements_all) {
  FACS_plotting_heatmap_line_graph_and_stats_along_time_points(
    input_path = file.path(FACS_data,"/"),
    measurement = i,
    output_path = file.path(FACS_plotting_plots_dir,"/")
  )
}

for (i in measurements_all) {
  FACS_plotting_heatmap_line_graph_and_stats_along_pseudotime(
    input_path = file.path(FACS_data,"/"),
    measurement = i,
    output_path = file.path(FACS_plotting_plots_dir,"/")
  )
}

##### Heatmap for Colon 
measurements_mfi <- c("CD63_MFI","MHCII_MFI","SigF_MFI","CD11b_MFI","CD80_MFI",
                      "CD107a_MFI","CD9_MFI","CD16_MFI","CD101_MFI","PDL1_MFI","SSCA_MFI")
measurements_freq <- c("Eos_freq","Active_eos_freq","CD9_freq","CD63_freq","CD80_freq","CD101_freq",
                       "CD16_freq","CD107a_freq","PDL1_freq")

## MFI measurements 
FACS_plotting_heatmap_tissue_of_interest_non_scaled(
  input_path = file.path(FACS_data,"/"),
  tissue_oi = "CO_" ,
  measurements = measurements_mfi,
  output_path = file.path(FACS_plotting_plots_dir,"/"),
  id = "MFI")

FACS_plotting_heatmap_tissue_of_interest_scaled(
  input_path = file.path(FACS_data,"/"),
  tissue_oi = "CO_" ,
  measurements = measurements_mfi,
  output_path = file.path(FACS_plotting_plots_dir,"/"),
  id = "MFI")

## freq measurements 
FACS_plotting_heatmap_tissue_of_interest_non_scaled(
  input_path = file.path(FACS_data,"/"),
  tissue_oi = "CO_" ,
  measurements = measurements_freq,
  output_path = file.path(FACS_plotting_plots_dir,"/"),
  id = "Freq")

FACS_plotting_heatmap_tissue_of_interest_scaled(
  input_path = file.path(FACS_data,"/"),
  tissue_oi = "CO_" ,
  measurements = measurements_freq,
  output_path = file.path(FACS_plotting_plots_dir,"/"),
  id = "Freq")


