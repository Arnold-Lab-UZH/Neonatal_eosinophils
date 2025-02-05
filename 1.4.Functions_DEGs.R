##### Functions for DEG analysis 
### Function to do DGE analysis between two conditions and writing the output to a csv file 
DEG_to_csv_two_cond <- function(
    seurat_object,
    cond1,
    cond2,
    only.pos_condition,
    logfc_threshold,
    csv_file_directory
){
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "RNA")
  DefaultAssay(seurat_object) <- "RNA"
  markers <- FindMarkers(object = seurat_object, ident.1 = cond1, ident.2 = cond2, only.pos = only.pos_condition, min.pct = 0.25, 
                         logfc.threshold = logfc_threshold,slot = "data")
  write.csv(markers, file = csv_file_directory)
}

### Function to do DEG analysis between two conditions using pseudobulking and DESeq2 
DEG_two_cond_pb_DESeq2 <- function(
    pseudobulk_object,
    celltype_cond1,
    celltype_cond2,
    output_path
){
  bulk_de <- FindMarkers(object = pseudobulk_object, 
                         ident.1 = celltype_cond1,
                         ident.2 = celltype_cond2,
                         test.use = "DESeq2")
  # write output 
  write.csv(bulk_de, paste0(output_path, "DESeq2_",celltype_cond1,"_",celltype_cond2,".csv"))
}

##### Functions for plotting of genes 
### Function to plot genes of interest per condition of interest, scaled per row 
heatmap_goi_coi <- function(
    seurat_object,
    condition_oi,
    markers_oi,
    assay_type,
    groups_of_markers,
    number_of_markers_per_group,
    colors_per_group,
    groups_and_colors,
    cluster_rows_cond,
    cluster_cols_cond
){
  # average expression per cluster and condition 
  average_expression <- AverageExpression(seurat_object, return.seurat = FALSE, features = markers_oi, normalization.method = "LogNormalize",assays = assay_type, group.by = condition_oi)
  average_expression_df <- as.data.frame(average_expression)
  average_expression_df <- average_expression_df[match(markers_oi, rownames(average_expression_df)),]
  
  # prepare palette for pheatmap
  paletteLength   <- 50
  myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
  breaksList      = seq(-2, 2, by = 0.04)
  
  # prepare annotation for pheatmap
  annotation_rows             <- data.frame(markers = rep(groups_of_markers, number_of_markers_per_group))
  rownames(annotation_rows)   <- rownames(average_expression_df)
  annotation_rows$markers     <- factor(annotation_rows$markers, levels = groups_of_markers)
  
  mycolors <- colors_per_group
  names(mycolors) <- unique(annotation_rows$markers)
  mycolors <- list(category = mycolors)
  annot_colors=list(markers=groups_and_colors)
  
  p <- pheatmap(average_expression_df,scale = "row",
                color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                breaks = breaksList,
                cluster_rows = cluster_rows_cond, cluster_cols = cluster_cols_cond, 
                border_color = "black", 
                legend_breaks = -2:2, 
                cellwidth = 10, cellheight = 5,
                angle_col = "45", 
                annotation_colors = annot_colors,
                annotation_row = annotation_rows,
                fontsize = 5)
  print(p)
}

### Volcano plot after DGE analysis between two conditions 
volcano_DGE_showing_goi <- function(
    DEG_list_path, 
    avg_log2FC_cutoff,
    p_val_adj_cutoff,
    condition1,
    condition2,
    genes_of_interest,
    colors_dots,
    x_line_intercept,
    xlim_values,
    output_path
){
  DEG_list <- read.csv(DEG_list_path)
  DEG_list <- DEG_list %>% mutate(
    Expression = case_when(avg_log2FC >= avg_log2FC_cutoff & p_val_adj <= p_val_adj_cutoff ~ condition1,
                           avg_log2FC <= -avg_log2FC_cutoff & p_val_adj <= p_val_adj_cutoff ~ condition2,
                           TRUE ~ "Non sig."))
  
  rownames(DEG_list) <- DEG_list$X
  
  DEG_list$genelabels <- ifelse(rownames(DEG_list) %in% genes_of_interest,TRUE,FALSE)
  
  p <- ggplot(DEG_list, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
    geom_point(size=5, aes(color = Expression)) +
    geom_text_repel(aes(label=ifelse(DEG_list$genelabels, rownames(DEG_list),"")),size=5, max.overlaps = 1000) +
    xlab("logFC") + 
    ylab("-log10(p_val_adj)") + 
    scale_color_manual(values =c(colors_dots) ,guide = "none") + 
    theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
    theme(plot.title = element_text(size = 25, face = "bold"))  + 
    #guides(colour = guide_legend(override.aes = list(size=1.5))) + 
    theme(legend.title = element_text(size = 25), legend.text = element_text(size = 25)) + 
    geom_vline(xintercept = x_line_intercept, linetype = "dashed", color = "gray50") + #log Fold Change cutoff 
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "gray50") + #p value cutoff 
    theme_classic(base_size = 25) + xlim(xlim_values)
  print(p)
  ggsave(output_path, width = 12, height = 10, plot = p)
}

### Function to plot average log FC in heatmap and print the ranges of p-values 
heatmap_logFC_goi <- function(
    deg_file_path,
    genes_of_interest,
    condition_oi,
    markers_oi,
    groups_of_markers,
    number_of_markers_per_group,
    colors_per_group,
    groups_and_colors
){
  deg_df <- read.csv(deg_file_path)
  
  deg_df <- deg_df[deg_df$X %in% genes_of_interest,]
  
  #change order of rows based on genes of interest vector 
  deg_df <- deg_df %>% arrange(factor(X, levels = genes_of_interest))
  
  #prepare palette for pheatmap
  paletteLength   <- 50
  myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
  breaksList      = seq(-2, 2, by = 0.04)
  
  #prepare annotation for pheatmap
  annotation_rows             <- data.frame(markers = rep(groups_of_markers, number_of_markers_per_group))
  rownames(annotation_rows)   <- rownames(deg_df)
  annotation_rows$markers     <- factor(annotation_rows$markers, levels = groups_of_markers)
  
  mycolors <- colors_per_group
  names(mycolors) <- unique(annotation_rows$markers)
  mycolors <- list(category = mycolors)
  annot_colors=list(markers=groups_and_colors)
  
  deg_df2 <- deg_df[,c(1,3)]
  rownames(deg_df2) <- deg_df2$X
  deg_df2$X <- NULL
  p <- pheatmap(deg_df2,
                color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                breaks = breaksList,
                cluster_rows = F, cluster_cols = F, 
                border_color = "black", 
                legend_breaks = -2:2, 
                cellwidth = 10, cellheight = 5,
                angle_col = "45", 
                annotation_colors = annot_colors,
                annotation_row = annotation_rows,
                fontsize = 5)
  print(p)
  
  # add range for p-values adjusted
  deg_df$smaller0.05_bigger0.01 <- ifelse(deg_df$p_val_adj <= 0.05 & deg_df$p_val_adj > 0.01,TRUE,FALSE)
  deg_df$smaller0.01_bigger0.001 <- ifelse(deg_df$p_val_adj <= 0.01 & deg_df$p_val_adj > 0.001,TRUE,FALSE)
  deg_df$smaller0.001 <- ifelse(deg_df$p_val_adj <= 0.001,TRUE,FALSE)
  
  deg_df <- deg_df[deg_df$p_val_adj <= 0.05,]
  
  print(deg_df)
}


