##### Function to run GSEA of KEGG, BP, Reactome and Hallmarks using singleseqset  
singleseqset_to_csv <- function(
    seurat_object,
    condition_meta_data_column,
    p_corr_method,
    output_file
){
  #normalize data 
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize",
                                 scale.factor = 10000, margin = 1, assay = "RNA")
  DefaultAssay(seurat_object) <- "RNA"
  
  #this calculates log fold change between clusters taking the normalized counts (slot data into account)
  logfc.data <- logFC(cluster.ids=seurat_object@meta.data[[condition_meta_data_column]],expr.mat=seurat_object[["RNA"]]$data)
  names(logfc.data)
  
  ###KEGG
  #calculation of enrichment score and p values 
  gse.res_kegg <- wmw_gsea(expr.mat=seurat_object[["RNA"]]$data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],
                           gene.sets=kegg.sets)
  
  #returns two item lists: test statistics and the p-values 
  res.stats_kegg <- gse.res_kegg[["GSEA_statistics"]]
  res.pvals_kegg <- gse.res_kegg[["GSEA_p_values"]]
  
  res.pvals_kegg <- apply(res.pvals_kegg,2,p.adjust,method=p_corr_method,n = length(rownames(res.pvals_kegg))) #Correct for multiple comparisons
  #number of p-values you need to correct for is the length of res.pvals_kegg
  
  res.pvals_df_kegg <- as.data.frame(res.pvals_kegg)
  
  res.stats_df_kegg <- as.data.frame(res.stats_kegg)
  
  #combine both to one data frame to then plot in bubbleplot 
  colnames(res.stats_df_kegg) <- paste0(colnames(res.stats_df_kegg), "_zscore")
  colnames(res.pvals_df_kegg) <- paste0(colnames(res.pvals_df_kegg), "_padj") 
  
  kegg_stats_pval <- cbind(res.pvals_df_kegg,res.stats_df_kegg)
  
  ###BP
  #calculation of enrichment score and p values 
  gse.res_bp <- wmw_gsea(expr.mat=seurat_object[["RNA"]]$data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],
                         gene.sets=bp.sets)
  
  #returns two item lists: test statistics and the p-values 
  res.stats_bp <- gse.res_bp[["GSEA_statistics"]]
  res.pvals_bp <- gse.res_bp[["GSEA_p_values"]]
  
  res.pvals_bp <- apply(res.pvals_bp,2,p.adjust,method="fdr",n = length(rownames(res.pvals_bp))) #Correct for multiple comparisons
  
  res.pvals_df_bp <- as.data.frame(res.pvals_bp)
  
  res.stats_df_bp <- as.data.frame(res.stats_bp)
  
  #combine both to one data frame to then plot in bubbleplot 
  colnames(res.stats_df_bp) <- paste0(colnames(res.stats_df_bp), "_zscore")
  colnames(res.pvals_df_bp) <- paste0(colnames(res.pvals_df_bp), "_padj") 
  
  bp_stats_pval <- cbind(res.pvals_df_bp,res.stats_df_bp)
  
  ###Reactome
  #calculation of enrichment score and p values 
  gse.res_R <- wmw_gsea(expr.mat=seurat_object[["RNA"]]$data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],
                        gene.sets=R.sets)
  
  #returns two item lists: test statistics and the p-values 
  res.stats_R <- gse.res_R[["GSEA_statistics"]]
  res.pvals_R <- gse.res_R[["GSEA_p_values"]]
  
  res.pvals_R <- apply(res.pvals_R,2,p.adjust,method="fdr") #Correct for multiple comparisons
  
  res.pvals_df_R <- as.data.frame(res.pvals_R)
  
  res.stats_df_R <- as.data.frame(res.stats_R)
  
  #combine both to one data frame to then plot in bubbleplot 
  colnames(res.stats_df_R) <- paste0(colnames(res.stats_df_R), "_zscore")
  colnames(res.pvals_df_R) <- paste0(colnames(res.pvals_df_R), "_padj") 
  
  R_stats_pval <- cbind(res.pvals_df_R,res.stats_df_R)
  
  ###Hallmarks
  #calculation of enrichment score and p values 
  gse.res_H <- wmw_gsea(expr.mat=seurat_object[["RNA"]]$data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],
                        gene.sets=H.sets)
  
  #returns two item lists: test statistics and the p-values 
  res.stats_H <- gse.res_H[["GSEA_statistics"]]
  res.pvals_H <- gse.res_H[["GSEA_p_values"]]
  
  res.pvals_H <- apply(res.pvals_H,2,p.adjust,method="fdr",n = length(rownames(res.pvals_H))) #Correct for multiple comparisons
  
  res.pvals_df_H <- as.data.frame(res.pvals_H)
  
  res.stats_df_H <- as.data.frame(res.stats_H)
  
  #combine both to one data frame to then plot in bubbleplot 
  colnames(res.stats_df_H) <- paste0(colnames(res.stats_df_H), "_zscore")
  colnames(res.pvals_df_H) <- paste0(colnames(res.pvals_df_H), "_padj") 
  
  H_stats_pval <- cbind(res.pvals_df_H,res.stats_df_H)
  
  #combine
  combined <- rbind(kegg_stats_pval, bp_stats_pval)
  combined <- rbind(combined, R_stats_pval)
  combined <- rbind(combined, H_stats_pval)
  write.csv(combined,output_file)
}

### Function to run progeny analysis between conditions 
PROGENy_two_cond <- function(
    seurat_object, 
    species,
    condition_column, 
    cond1,
    cond2,
    sample_name
){
  ### compute Progeny activity score per condition
  ## ad a new assay called Progeny
  seurat_object <- progeny(seurat_object, scale=FALSE, organism=species, top=500, perm=1, return_assay = TRUE)
  
  ## scale the pathway activity scores. 
  seurat_object <- Seurat::ScaleData(seurat_object, assay = "progeny") 
  
  ## create a data frame with the specification of the cells that belong to 
  ## each cluster to match with the Progeny scores.
  Idents(seurat_object) <- condition_column
  CellsClusters <- data.frame(Cell = names(Idents(seurat_object)),Type = as.character(Idents(seurat_object)),stringsAsFactors = FALSE)
  
  ## transform Progeny scores into a data frame to better handling the results
  progeny_scores_df <- as.data.frame(t(GetAssayData(seurat_object, layer = "scale.data", assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%gather(Pathway, Activity, -Cell) 
  
  ## match Progeny scores with the cell clusters.
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
  
  ## We summarize the Progeny scores by cell population
  summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, Type) %>%
    summarise(mean = mean(Activity), std = sd(Activity))
  
  ## We prepare the data for the plot
  summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, mean) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
  
  print(ComplexHeatmap::Heatmap(summarized_progeny_scores_df, name=sample_name,
                                column_names_gp = grid::gpar(fontsize = 10),
                                row_names_gp = grid::gpar(fontsize = 10),
                                column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                                column_names_rot = 45,cluster_rows = FALSE, cluster_columns = FALSE,column_title = "GO terms"))
  
  ## caclualte the difference of scaled scores 
  summarized_progeny_scores_df_t <- as.data.frame(t(summarized_progeny_scores_df))
  print(head(summarized_progeny_scores_df_t))
  summarized_progeny_scores_df_t$diff1 <- summarized_progeny_scores_df_t[[cond1]] - summarized_progeny_scores_df_t[[cond2]]
  
  summarized_progeny_scores_df_t[[cond1]] <- NULL
  summarized_progeny_scores_df_t[[cond2]] <- NULL
  
  summarized_progeny_scores_df_t <- as.matrix(summarized_progeny_scores_df_t)
  
  print(ComplexHeatmap::Heatmap(summarized_progeny_scores_df_t, name=paste0("Diff mean scaled progeny score ",sample_name),
                                column_names_gp = grid::gpar(fontsize = 10),
                                row_names_gp = grid::gpar(fontsize = 10),
                                column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                                column_names_rot = 45,cluster_rows = FALSE, cluster_columns = FALSE,column_title = "GO terms"))
  
  ### statistics 
  ## genearate a df with means and p-values for 
  specific_signaling_df <- progeny_scores_df[progeny_scores_df$Pathway %in% "Androgen",]
  
  ## test weather the variance is different between the conditions
  var.test(Activity ~ Type, specific_signaling_df, alternative = "two.sided")
  
  ## test if data is normally distributed 
  ggqqplot(specific_signaling_df$Activity)
  
  ## apply Welch's T-test two-sided and correct p-values with Bonferroni 
  pathways <- c("Androgen","EGFR","Estrogen","Hypoxia","JAK-STAT","MAPK","NFkB","p53","PI3K","TGFb","TNFa","Trail","VEGF","WNT")
  
  for (i in pathways) {
    specific_signaling_df <- progeny_scores_df[progeny_scores_df$Pathway %in% i,]
    p.value <- t.test(specific_signaling_df[specific_signaling_df$Type %in% cond1,]$Activity, 
                      specific_signaling_df[specific_signaling_df$Type %in% cond2,]$Activity, alternative = "two.sided")$p.value
    p.adjust <- p.adjust(p.value, method = "bonferroni", n = length(rownames(summarized_progeny_scores_df_t))) #correct for multiple comparisons 
    print(paste0(i, "; bonferroni_adj_p: ",p.adjust))
  }
}
