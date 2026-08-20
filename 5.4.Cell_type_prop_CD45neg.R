########## This code looks at cell type proportions of CD45 negative cells from CO neo P14 PHIL and WT ##########

##### link to libraries and functions
source("1.1.config.R")
source(file.path(base_dir,"1.3.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "1.2.Packages.R"))
source(file.path(base_dir, "1.7.Functions_cell_type_prop.R"))

### load objects 
obj <- readRDS(file.path(seurat_objects_dir,"Neo_P14_CD45neg_colon_WT_PHIL_anno.rds"))

##### umap split by condition 
p <- DimPlot(obj, group.by = "annotation", pt.size = 0.5, label = FALSE, label.size = 5,split.by = "condition", 
             cols =c("#EA1351","#755307","#0CEA75","#08592E","#EA86C6","#EF9818","#8810E5","#F2D3A5","#5220EF","#119CEA","#C5E0EF","#DAEF18") )
ggsave(file.path(integration_eos_plots_dir,"umap_annotated_CD45neg.svg"), width = 20, height = 8, plot = p)

##### proportion plot 
create_table_cell_type_prop(obj, "condition","annotation",file.path(integration_eos_tables_dir),"/CD45neg")
df <- read.csv(file.path(integration_eos_tables_dir,"CD45neg_proportions_condition_annotation.csv"), header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:13)) 

#define order
df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c("TA","Enterocytes_Alpi_high","Enterocytes_Alpi_low", "Goblet","Reg4_secretory", "EECs","Glial","Endothelial",
                                                                   "Pericytes","SMCs","Stromal_Cd34" ,"Stromal_Pdgfra")))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values = c("#DAEF18","#0CEA75","#08592E","#EF9818","#F2D3A5","#EA1351","#EA86C6","#755307", "#8810E5","#5220EF","#119CEA","#C5E0EF") ) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave(file.path(integration_eos_plots_dir,"cell_type_prop_CD45neg.svg"), width = 12, height = 6, plot = p)

##### statistics 
#first control, then the condition
cell_type_prop_stats(obj,"annotation","neo_P14_wt","neo_P14_phil","condition",1.2,
                     file.path(integration_eos_plots_dir,"cell_type_prop_CD45neg_neo_wt_philstat.svg"))



