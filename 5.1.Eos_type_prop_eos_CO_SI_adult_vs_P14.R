########## This code looks at eosinophils subtype proportions between neo P14 CO/SI and adults (Gurtner et al) ##########
# Figure 3 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.5.Functions_cell_type_prop.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_adult_eos_CO_SI_LT.rds")

# extract only CO and SI
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("adult_small_int","adult_colon","NEO_P14_colon","NEO_P14_small_int"))

##### umap of annotated object split by condition
p <- DimPlot(obj, group.by = "annotation",split.by= "condition",cols = c("#E81818","#10A069","#E8E81A", "#E88A1A","#26DFED"))
ggsave(file = "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/umap_split_by_cond.svg", plot = p, width = 10, height = 6)

##### proportion plot 
create_table_cell_type_prop(obj, "condition","annotation","/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/","eos")
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/eos_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:6))

df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c("A_eos","B_eos","circ_eos", "immature_eos","precursor_eos")))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E81818","#10A069","#E8E81A", "#E88A1A", "#26DFED")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/cell_type_prop.svg", width = 12, height = 6, plot = p)

##### statistics 
#first control, then the condition
cell_type_prop_stats(obj,"annotation","adult_colon","NEO_P14_colon","condition",1.2,
                     "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/adult_CO_vs_neo_CO_stat.pdf") 

cell_type_prop_stats(obj,"annotation","adult_small_int","NEO_P14_small_int","condition",1.2,
                     "/scratch/khandl/Neonatal_eosinophils/figures/CO_SI_LT/adult_SI_vs_neo_SI_stat.pdf") 

