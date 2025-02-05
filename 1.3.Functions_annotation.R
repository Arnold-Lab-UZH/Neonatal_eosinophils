##### Functions for cell type annotations 
### Function to project the cell type from SingleR result to the umap space to identify which cluster represents which cell type 
project_annotation_to_umap <- function(cell.type, singleResult, seurat_object) {
  temp <- as.data.frame(singleResult[4])
  temp$cell <- rownames(temp)
  temp <- temp%>%filter(temp$pruned.labels %in% cell.type)
  temp <- temp$cell
  print(DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 10, pt.size = 0.5, cells.highlight = temp, sizes.highlight = 2) + NoLegend() + ggtitle(cell.type))
}

## after FastMNN integration 
project_annotation_to_umap_fastMNN <- function(cell.type, singleResult, seurat_object) {
  temp <- as.data.frame(singleResult[4])
  # singleResult[5] is extracting pruned.labels from data frame
  temp$cell <- rownames(temp)
  temp <- temp%>%filter(temp$pruned.labels %in% cell.type)
  temp <- temp$cell
  print(DimPlot(seurat_object, reduction = "umap.mnn", label = TRUE, label.size = 10, pt.size = 2, cells.highlight = temp, sizes.highlight = 2,raster=FALSE) + NoLegend() + ggtitle(cell.type))
}

