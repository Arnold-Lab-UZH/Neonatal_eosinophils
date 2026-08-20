##### Functions for pre-processing 
### Function to read in data generated from BD Seven Bridges platform (by Simona Baghai)
data_to_sparse_matrix <- function(data.st_file_path) {
  # read in file with cell index - gene name - values
  # import for one cartridge, one sample
  input <-read.table(data.st_file_path, header = T)
  # transform to matrix (data.frame actually)
  # we take as default values from the column "RSEC_Adjusted_Molecules" (= error corrected UMIs)
  mat <- input %>% pivot_wider(id_cols = Bioproduct, 
                               values_from = RSEC_Adjusted_Molecules, 
                               names_from = Cell_Index, values_fill = 0)  %>% 
    tibble::column_to_rownames("Bioproduct")
  # convert to sparse matrix (~ dgMatrix)
  sparse_mat = Matrix(as.matrix(mat),sparse=TRUE)
  return(sparse_mat)
}

### Function to generate Seurat objects for each sample 
create_seurat_from_condition <- function(
    path_to_st_file,
    project,
    condition, 
    min.cells,
    min.features
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  condition_sample$condition <- condition
  
  return(condition_sample)
}

### Function to generate Seurat objects for each sample + DecontX to account for cell free RNA 
create_seurat_from_condition_DecontX <- function(
    path_to_st_file,
    project,
    condition, 
    min.cells,
    min.features
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  sce <- as.SingleCellExperiment(condition_sample)
  ### run decontX with default settings 
  sce.delta <- decontX(sce)
  ### convert back to a Seurat object 
  seuratObject <- CreateSeuratObject(round(decontXcounts(sce.delta)))
  seuratObject$condition <- condition
  
  return(seuratObject)
}

## for human PB and CB data 
create_seurat_plus_DecontX_human_biopsies <- function(
    path_to_st_file,
    project,
    min.cells,
    min.features,
    condition_oi,
    tissue_oi,
    experiment_oi,
    phenotype_oi
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  sce <- as.SingleCellExperiment(condition_sample)
  # run decontX with default settings 
  sce.delta <- decontX(sce)
  # convert back to a Seurat object 
  seurat_object <- CreateSeuratObject(round(decontXcounts(sce.delta)))
  seurat_object$condition <- condition_oi
  seurat_object$tissue <- tissue_oi
  seurat_object$experiment <- experiment_oi
  seurat_object$phenotype <- phenotype_oi
  return(seurat_object)
}

## to pre-process data from Gurtner et al. the old WTA v1 beads were used for BD Rhapsody 
data_to_sparse_matrix_old_WTA_version <- function(data.st_file_path) {
  # read in file with cell index - gene name - values
  # import for one cartridge, one sample
  input <-read.table(data.st_file_path, header = T)
  # transform to matrix (data.frame actually)
  # we take as default values from the column "RSEC_Adjusted_Molecules" (= error corrected UMIs)
  mat <- input %>% pivot_wider(id_cols = Gene, 
                               values_from = RSEC_Adjusted_Molecules, 
                               names_from = Cell_Index, values_fill = 0)  %>% 
    tibble::column_to_rownames("Gene")
  # convert to sparse matrix (~ dgMatrix)
  sparse_mat = Matrix(as.matrix(mat),sparse=TRUE)
  return(sparse_mat)
}

create_seurat_from_condition_old_WTA_version <- function(
    path_to_st_file,
    project,
    condition, 
    min.cells,
    min.features
) {
  input_matrix <- data_to_sparse_matrix_old_WTA_version(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  condition_sample$condition <- condition
  
  return(condition_sample)
}
