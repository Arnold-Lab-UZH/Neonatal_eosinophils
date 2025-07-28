##### Functions for CellPhoneDB 
### Function for generation of input files from mouse data 
#first gene symbols have to be converted from mouse to human because CellPhoneDB database only contains human L-R interactions   
#Conversion of mouse to human genes adapted partially from: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md
#and: https://www.cellphonedb.org/faq-and-troubleshooting
#outputs are two text tiles: gene counts and cell annotations (meta)
Input_files_CellPhoneDB_generation_mm <- function(
    seurat_object,
    annotation_column,
    sample_name,
    ouput_file_path
){
  ### load human and mouse ensemble symbols
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  ###generating counts file 
  #take raw data and normalize it
  count_raw_meta <- GetAssayData(object = seurat_object, layer = "counts")[,colnames(x = seurat_object)]
  count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(count_norm_meta) , mart = mouse, 
                   attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
  print(head(genesV2))
  matrixA <- count_norm_meta[match(genesV2$MGI.symbol,rownames(count_norm_meta),nomatch=T),]
  matrixB <- matrixA
  matrixB$gene <- genesV2$Gene.stable.ID
  rownames(matrixA) <- matrixB$gene
  #save count matrix as text file 
  write.table(matrixA, paste0(ouput_file_path,sample_name,"_count.txt"), sep='\t', quote=F)
  ###generating meta file based on cell type annotation of Seurat object 
  meta_data_meta <- cbind(rownames(seurat_object@meta.data), seurat_object@meta.data[,annotation_column, drop=F])  
  #save meta file as text file 
  write.table(meta_data_meta, paste0(ouput_file_path,sample_name,"_meta.txt"), sep='\t', quote=F, row.names=F)
}

### Function for generation of CellPhoneDB input files from human data 
Input_files_CellPhoneDB_generation_hs <- function(
    seurat_object,
    annotation_column,
    sample_name,
    ouput_file_path
){
  ###generating counts file 
  #take raw data and normalize it
  count_raw_meta <- GetAssayData(object = seurat_object, layer = "counts")[,colnames(x = seurat_object)]
  count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm_meta, paste0(ouput_file_path,sample_name,"_count.txt"), sep='\t', quote=F)
  ###generating meta file based on cell type annotation of Seurat object 
  meta_data_meta <- cbind(rownames(seurat_object@meta.data), seurat_object@meta.data[,annotation_column, drop=F])  
  #save meta file as text file 
  write.table(meta_data_meta, paste0(ouput_file_path,sample_name,"_meta.txt"), sep='\t', quote=F, row.names=F)
}

### Function to compare two conditions 
LR_per_2_cond_interacting_pair_oi <- function(
    file_path1,
    file_path2,
    interacting_pair_oi, 
    cond1, 
    cond2,
    output_path
){
  # read in files 
  file1 <- read.delim(file_path1, check.names = FALSE)
  file2 <- read.delim(file_path2, check.names = FALSE)
  
  # make na to 0 
  file1[is.na(file1)] <- 0
  file2[is.na(file2)] <- 0
  
  # extract interacting pair oi 
  file1_pair_oi <- file1[,colnames(file1) %in% c(interacting_pair_oi, "interacting_pair","secreted","receptor_a","receptor_b","annotation_strategy", "is_integrin","directionality","classification")]
  file2_pair_oi <- file2[,colnames(file2) %in% c(interacting_pair_oi, "interacting_pair","secreted","receptor_a","receptor_b","annotation_strategy", "is_integrin","directionality","classification")]
  
  # change colname of mean expression to condition name 
  colnames(file1_pair_oi) <- c("interacting_pair","secreted","receptor_a","receptor_b","annotation_strategy", "is_integrin","directionality","classification",cond1)
  colnames(file2_pair_oi) <- c("interacting_pair","secreted","receptor_a","receptor_b","annotation_strategy", "is_integrin","directionality","classification", cond2)
  
  # merge the two dataframes and add 0 if interaction is missing in one 
  df <- merge(file1_pair_oi, file2_pair_oi, all = TRUE)
  df[is.na(df)] <- 0
  
  #remove interactions with 0 in all 
  df$sum <- rowSums(df[9:10])
  df2 <- df[which(df$sum > 0),]
  #remove sum column
  df2$sum <- NULL
  
  # write combined output 
  write.table(df2, file = paste0(output_path,interacting_pair_oi,"all_2_cond" ,".txt"), sep = "\t", dec=".", quote = FALSE, row.names = FALSE)
}

### function to count the total of ligand and receptor interactions between two cell types of interest 
ligand_receptor_counts <- function(
    means_table, 
    interacting_partner_1,
    interacting_partner_2
){
  means <- means_table[means_table$directionality == "Ligand-Receptor",]
  eos_other <- means[,colnames(means) %in% c("interacting_pair",paste0(interacting_partner_1,"|",interacting_partner_2))]
  eos_other <- eos_other[eos_other[[paste0(interacting_partner_1,"|",interacting_partner_2)]] > 0,]
  int_1_counts <- length(rownames(eos_other))
  
  print(int_1_counts)
}

adhesion_adhesion_counts <- function(
    means_table, 
    interacting_partner_1,
    interacting_partner_2
){
  means <- means_table[means_table$directionality == "Adhesion-Adhesion",]
  eos_other <- means[,colnames(means) %in% c("interacting_pair",paste0(interacting_partner_1,"|",interacting_partner_2))]
  eos_other <- eos_other[eos_other[[paste0(interacting_partner_1,"|",interacting_partner_2)]] > 0,]
  int_1_counts <- length(rownames(eos_other))
  
  print(int_1_counts)
}

