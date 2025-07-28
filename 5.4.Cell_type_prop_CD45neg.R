########## This code looks at cell type proportions of CD45 negative cells from CO neo P14 PHIL and WT ##########
# Figure 5, S5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.5.Functions_cell_type_prop.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL_anno2.rds")

##### umap split by condition 
p <- DimPlot(obj, group.by = "annotation", pt.size = 0.5, label = FALSE, label.size = 5,split.by = "condition", 
             cols =c("#EA1351","#755307","#0CEA75","#08592E","#EA86C6","#EF9818","#8810E5","#F2D3A5","#5220EF","#119CEA","#C5E0EF","#DAEF18") )
ggsave("/scratch/khandl/Neonatal_eosinophils/CD45neg/umap_annotated.svg", width = 20, height = 8, plot = p)

##### proportion plot 
create_table_cell_type_prop(obj, "condition","annotation","/scratch/khandl/Neonatal_eosinophils/CD45neg/","CD45neg")
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/CD45neg/CD45neg_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:13)) 

#define order
df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c("TA","Enterocytes_Alip_high","Enterocytes_Alpi_low", "Goblet","Reg4_secretory", "EECs","Glial","Endothelial",
                                                                   "Pericytes","SMCs","Stromal_Cd34" ,"Stromal_Pdgfra")))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values = c("#DAEF18","#0CEA75","#08592E","#EF9818","#F2D3A5","#EA1351","#EA86C6","#755307", "#8810E5","#5220EF","#119CEA","#C5E0EF") ) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/Neonatal_eosinophils/CD45neg/cell_type_prop.svg", width = 12, height = 6, plot = p)

##### statistics 
#first control, then the condition
cell_type_prop_stats(obj,"annotation","neo_P14_wt","neo_P14_phil","condition",1.2,
                     "/scratch/khandl/Neonatal_eosinophils/CD45neg/cell_type_prop_stat.pdf") 



