########## This code compares gene expression between CO and SI eosinophils of adults and neo P14 ##########
# Figure S3

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")

##### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_LT.rds")

##### generate data frame with average expression of genes per condition
average_expression <- AverageExpression(obj, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "condition")
average_expression_df <- as.data.frame(average_expression)
average_expression_df$gene <- rownames(average_expression_df)

##### adult 
### CO 
colon <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.adult.colon","gene")] 
#extract genes with >0 gene expression 
colon <- colon[colon$RNA >0,]

### SI 
small_int <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.adult.small.int","gene")] 
#extract genes with >0 gene expression 
small_int <- small_int[small_int$RNA >0,]

### venn diagram 
x <- list("CO" = colon$gene, "SI" =small_int$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

##### neo P14  
### CO 
colon <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.NEO.P14.colon","gene")] 
#extract genes with >0 gene expression 
colon <- colon[colon$RNA >0,]

### SI 
small_int <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.NEO.P14.small.int","gene")] 
#extract genes with >0 gene expression 
small_int <- small_int[small_int$RNA >0,]

### venn diagram 
x <- list("CO" = colon$gene, "SI" =small_int$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 
