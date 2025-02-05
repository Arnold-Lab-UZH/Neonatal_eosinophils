########## This code compares CO eosinophils between neo P16 after eos depletion (day 5, 7, 10) vs. neo P16 control without eos depletion ##########
# Figure 6, S6

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")

##### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P16_iDT_CREpos_CREneg_anno.rds")

## extract eosinophils 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = "Eosinophils")

##### Compare number of expressed genes 
average_expression <- AverageExpression(obj, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "condition")
average_expression_df <- as.data.frame(average_expression)
average_expression_df$gene <- rownames(average_expression_df)

## Cre - (ctrl)
cre_neg <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.P16.colon.iDT.Cre.neg","gene")] 
#extract genes with >0 gene expression 
cre_neg <- cre_neg[cre_neg$RNA >0,]

## Cre - (ctrl)
cre_pos <- average_expression_df[,colnames(average_expression_df) %in% c("RNA.P16.colon.iDT.Cre.pos","gene")] 
#extract genes with >0 gene expression 
cre_pos <- cre_pos[cre_pos$RNA >0,]

## venn diagram 
x <- list("Cre+" = cre_pos$gene, "Cre-" =cre_neg$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

##### Correlation plot 
## gene expression 
df <- average_expression_df
df$gene <- NULL
#only consider genes with counts >0 in either one of the conditions 
df <- df[df$RNA.P16.colon.iDT.Cre.neg >0 | df$RNA.P16.colon.iDT.Cre.pos >0, ]

p <- ggplot(df, aes(x = RNA.P16.colon.iDT.Cre.pos , y = RNA.P16.colon.iDT.Cre.neg)) + stat_poly_eq(use_label(c("eq", "R2")))+ stat_poly_line() + 
  geom_point()+ theme_classic(base_size = 20) 
ggsave("/scratch/khandl/eos_NEO/NEO58/Gene_expr_corr.svg", width = 12, height = 6, plot = p)

