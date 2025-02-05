########## This code looks at eosinophils subtype proportions between CB and PB ##########
# Figure 2

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.5.Functions_cell_type_prop.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/CB_PB_eos_integrated_anno.rds")

##### proportion of annotated cluster per condition
create_table_cell_type_prop(obj, "type","annotation","/scratch/khandl/Neonatal_eosinophils/figures/CB_PB_eos/","eos_anno")
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/figures/CB_PB_eos/eos_anno_proportions_type_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2,3)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#2AB34B","#7094CD")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/CB_PB_eos/prop_barplot_clusters.svg", width = 12, height = 6, plot = p)

##### statistic 
cell_type_prop_stats(obj,"annotation","PB","CB","type",1.41,
                     "/scratch/khandl/Neonatal_eosinophils/figures/CB_PB_eos/stats_adult_cord.svg") 

