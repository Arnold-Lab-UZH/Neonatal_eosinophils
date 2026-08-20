########## Configuration file ##########
# Source this file at the beginning of each analysis script to define
# project-relative input and output paths.

##### Define the Project base directory (directory containing this file) to load packages and function R code 
base_dir <- tryCatch({dirname(sys.frames()[[1]]$ofile)}, error = function(e) {getwd()})

##### Configure raw data directories 
scratch_dir <- Sys.getenv("scratch_dir", "/home/khandl/scratch")

### scRNAseq data 
raw_data_dir      <- file.path(scratch_dir, "raw")
raw_Mm_il5tg_data <- file.path(raw_data_dir, "Mm_il5tg")
raw_Mm_wt_data <- file.path(raw_data_dir, "Mm_wt")
raw_Hs_data <- file.path(raw_data_dir, "Hs")
raw_Mm_eos_Cre_data <- file.path(raw_data_dir, "Mm_eos_Cre")

### bulk RNAseq data 
raw_Mm_bulk_RNAseq_data <- file.path(raw_data_dir, "Mm_bulkRNAseq")

### FACS data 
FACS_data <- file.path(raw_data_dir, "FACS")

##### Configure directories for supporting files 
SCENIC_support_files <- "/home/khandl/data/common/SCENIC"
CellPhoneDB_support_files <- "/home/khandl/data/cpdb/v5.0.0"

##### Configure data output directories (relative to repository)
data_dir          <- Sys.getenv("data_dir", "/home/khandl/data/Neonatal_eosinophils")
seurat_objects_dir <- file.path(data_dir, "seurat_objects")
seurat_objects_Gurtner_et_al <- file.path(data_dir, "seurat_objects_Gurtner_et_al")

##### Configure results output directories (plots and tables)
### Root directory 
results_dir <- file.path(data_dir, "results")

## 3. Annotation  
annotation_plots_dir <- file.path(results_dir, "3.annotation", "plots")

## 4. Integration of eosinophils from different datasets 
integration_eos_tables_dir <- file.path(results_dir, "4.integration_eos", "tables")
integration_eos_plots_dir <- file.path(results_dir, "4.integration_eos", "plots")

## 5. Analysis of cell type proportion differences between conditions
cell_type_prop_tables_dir <- file.path(results_dir, "5.cell_type_prop", "tables")
cell_type_prop_plots_dir <- file.path(results_dir, "5.cell_type_prop", "plots")

## 6. BM derived precursors
BM_precursor_tables_dir <- file.path(results_dir, "6.BM_precursors", "tables")
BM_precursor_plots_dir <- file.path(results_dir, "6.BM_precursors", "plots")

## 7. Pseudotime trajectory eos
pseudotime_trajectory_eos_tables_dir <- file.path(results_dir, "7.pseudotime_trajectory_eos", "tables")
pseudotime_trajectory_eos_plots_dir <- file.path(results_dir, "7.pseudotime_trajectory_eos", "plots")

## 8.DEG analysis 
DEGs_tables_dir <- file.path(results_dir, "8.DEG_analysis", "tables")
DEGs_plots_dir <- file.path(results_dir, "8.DEG_analysis", "plots")
DEGs_tables_CD45pos_all_cell_types_dir <- file.path(DEGs_tables_dir, "CD45pos")
DEGs_tables_CD45neg_all_cell_types_dir <- file.path(DEGs_tables_dir, "CD45neg")

## 9.Signature scores
Signature_scores_tables_dir <- file.path(results_dir, "9.signature_scores", "tables")
Signature_scores_plots_dir <- file.path(results_dir, "9.signature_scores", "plots")

## 10.Signaling pathway analyses
Signaling_pathways_tables_dir <- file.path(results_dir, "10.signaling_pathways", "tables")
Signaling_pathways_plots_dir <- file.path(results_dir, "10.signaling_pathways", "plots")
CellPhoneDB_tables_dir <- file.path(Signaling_pathways_tables_dir, "CellPhoneDB")
CellPhoneDB_neo_eos_CD45neg_output_tables_dir <- file.path(CellPhoneDB_tables_dir, "neo_eos_CD45neg_output")
CellPhoneDB_neo_CD45neg_phil_output_tables_dir <- file.path(CellPhoneDB_tables_dir, "neo_CD45neg_phil_output")
CellPhoneDB_CD45neg_phil_vs_wt_output_tables_dir <- file.path(CellPhoneDB_tables_dir, "CD45neg_phil_vs_wt")

## 11.Comparison of gene expression in eos between conditions and datasets 
Gene_expr_comp_tables_dir <- file.path(results_dir, "11.eos_gene_expr_comp", "tables")
Gene_expr_comp_plots_dir <- file.path(results_dir, "11.eos_gene_expr_comp", "plots")

## 12. Bulk RNA-seq analysis 
bulkRNAseq_tables_dir <- file.path(results_dir, "12.bulkRNAseq", "tables")
bulkRNAseq_plots_dir <- file.path(results_dir, "12.bulkRNAseq", "plots")

## 13. FACS data plotting 
FACS_plotting_tables_dir <- file.path(results_dir, "13.FACS_plotting", "tables")
FACS_plotting_plots_dir <- file.path(results_dir, "13.FACS_plotting", "plots")


