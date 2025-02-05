########## This code looks at cell type proportions of CD45 negative cells from CO neo P14 PHIL and WT ##########
# Figure 5, S5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.5.Functions_cell_type_prop.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Neo_P14_CD45neg_colon_WT_PHIL_anno.rds")

##### umap split by condition 
p <- DimPlot(obj, group.by = "annotation", pt.size = 0.5, label = TRUE, label.size = 5,split.by = "condition",
             cols = c("#1128CE", "#21932C", "#11BCCE","#ACEFB1","#81F20A",
                      "#937021","#F70B32", "#F294E7","#F2F20A"))
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/CD45neg_prop/umap_annotated.svg", width = 8, height = 5, plot = p)

#colors: 
#SMCs: #1128CE
#Fibroblasts: #11BCCE
#Myofibroblasts: #F70B32
#Pericytes: #F294E7
#Endothelial: #937021
#Colonocytes: #21932C
#TA/enterocytes: #ACEFB1
#Goblet: #81F20A
#EECs: #F2F20A

##### proportion plot 
create_table_cell_type_prop(obj, "condition","annotation","/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/","CD45neg")
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/CD45neg_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(3:8, 10:12)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values = c("#21932C", "#F2F20A", "#937021","#11BCCE","#81F20A",
                               "#F70B32","#F294E7", "#1128CE","#ACEFB1")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/cell_type_prop.svg", width = 12, height = 6, plot = p)

##### statistics 
#first control, then the condition
cell_type_prop_stats(obj,"annotation","neo_P14_wt","neo_P14_phil","condition",1.2,
                     "/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/cell_type_prop_stat.pdf") 



