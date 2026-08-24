########## This code looks at eosinophils subtype proportions between CB and PB ##########

##### link to libraries and functions
source("01.01.config.R")
source(file.path(base_dir,"01.03.Output_directory_output_folder_structure_generation.R"))
source(file.path(base_dir, "01.02.Packages.R"))
source(file.path(base_dir, "01.07.Functions_cell_type_prop.R"))

### load objects 
obj <- readRDS(file.path(seurat_objects_dir,"CB_PB_eos_integrated_anno.rds"))

##### proportion of annotated cluster per condition
create_table_cell_type_prop(obj, "type","annotation",file.path(integration_eos_tables_dir),"/eos_CB_PB")
df <- read.csv(file.path(integration_eos_tables_dir,"eos_CB_PB_proportions_type_annotation.csv"), header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2,3)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#2AB34B","#7094CD")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave(file.path(integration_eos_plots_dir,"prop_barplot_clusters.svg"), width = 12, height = 6, plot = p)

##### statistic 
cell_type_prop_stats(obj,"annotation","PB","CB","type",1.41,
                     file.path(integration_eos_plots_dir,"stats_adult_cord.svg") )

