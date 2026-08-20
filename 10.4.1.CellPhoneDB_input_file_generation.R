########## This code generates input files for CellPhoneDB to predict ligand-receptor interactions between neonatal eosinophils and CD45neg cell types ##########

##### link to libraries and functions
source("1.1.config.R")
source(file.path(base_dir,"1.3.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "1.2.Packages.R"))
source(file.path(base_dir, "1.10.Functions_CellPhoneDB.R"))

#### load human and mouse ensemble symbols
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

##### CD45 negative wt and eosinophils 
### load objects and merge eosinophils with CD45 negative wt dataset 
CD45neg <- readRDS(file.path(seurat_objects_dir,"Neo_P14_CD45neg_colon_WT_PHIL_anno.rds"))
eos <- readRDS(file.path(seurat_objects_dir,"Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds"))

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

### input file generation 
Input_files_CellPhoneDB_generation_mm(obj, "annotation","neo", file.path(CellPhoneDB_tables_dir,"/")) 

##### CD45 negative phil 
### load object and extract only phil  
CD45neg <- readRDS(file.path(seurat_objects_dir,"Neo_P14_CD45neg_colon_WT_PHIL_anno2.rds"))

Idents(CD45neg) <- "condition"
CD45neg <- subset(CD45neg, idents = "neo_P14_phil")

### input file generation 
Input_files_CellPhoneDB_generation_mm(CD45neg, "annotation","CD45_neo_phil",file.path(CellPhoneDB_tables_dir,"/")) 
