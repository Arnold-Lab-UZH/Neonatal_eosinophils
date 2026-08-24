########## This code creates all directories needed to save output files based on 1.1.config.R ##########
### Function to create the folders 
create_project_dirs <- function(base_dir = getwd()) {
  # Define the base directory 
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  
  # Top-level output roots, exactly as used in the configuration 
  data_root    <- file.path(base_dir)
  results_root <- file.path(base_dir, "results")
  
  # All directories defined in 1.1.config.R
  dirs <- c(
    file.path(data_root, "seurat_objects"),
    
    file.path(results_root, "3.annotation", "plots"),
    
    file.path(results_root, "4.integration_eos", "plots"),
    file.path(results_root, "4.integration_eos", "tables"),
    
    file.path(results_root, "5.cell_type_prop", "plots"),
    file.path(results_root, "5.cell_type_prop", "tables"),
    
    file.path(results_root, "6.BM_precursors", "plots"),
    file.path(results_root, "6.BM_precursors", "tables"),
    
    file.path(results_root, "7.pseudotime_trajectory_eos", "plots"),
    file.path(results_root, "7.pseudotime_trajectory_eos", "tables"),
    
    file.path(results_root, "8.DEG_analysis", "plots"),
    file.path(results_root, "8.DEG_analysis", "tables"),
    file.path(results_root, "8.DEG_analysis","tables", "CD45pos"),
    file.path(results_root, "8.DEG_analysis","tables", "CD45neg"),
    
    file.path(results_root, "9.signature_scores", "plots"),
    file.path(results_root, "9.signature_scores", "tables"),
    
    file.path(results_root, "10.signaling_pathways", "plots"),
    file.path(results_root, "10.signaling_pathways", "tables"),
    file.path(results_root, "10.signaling_pathways", "tables","CellPhoneDB"),
    file.path(results_root, "10.signaling_pathways", "tables","CellPhoneDB","neo_eos_CD45neg_output"),
    file.path(results_root, "10.signaling_pathways", "tables","CellPhoneDB","neo_CD45neg_phil_output"),
    file.path(results_root, "10.signaling_pathways", "tables","CellPhoneDB","CD45neg_phil_vs_wt"),
    
    file.path(results_root, "11.eos_gene_expr_comp", "plots"),
    file.path(results_root, "11.eos_gene_expr_comp", "tables"),
    
    file.path(results_root, "12.bulkRNAseq", "plots"),
    file.path(results_root, "12.bulkRNAseq", "tables"),
    
    file.path(results_root, "13.FACS_plotting", "plots"),
    file.path(results_root, "13.FACS_plotting", "tables")
  )
  
  created <- vapply(dirs, function(d) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
  }, FUN.VALUE = character(1))
  
  message(
    length(created),
    " project output directories ensured under: ",
    base_dir
  )
  
  invisible(as.character(created))
}

### Create everything under a directory of your choice
create_project_dirs("/home/khandl/data/Neonatal_eosinophils")

