########## This code looks at cell type proportions of CD45 positive cells from CO neo P14 PHIL and WT and adult WT ##########
# Figure 5, S5

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")
source("~/Projects/Neonatal_eosinophils/1.5.Functions_cell_type_prop.R")

### load objects 
obj <- readRDS("/data/khandl/Neonatal_eosinophils/seurat_objects/Adult_and_Neo_P14_CD45enr_colon_WT_PHIL_anno.rds")

## extract only immune cells 
Idents(obj) <- "annotation"
immune <- subset(obj, idents = c("CD4_T","CD8_CCR7_T","CD8T_NKT","cDC1","cDC2","pDCs","act_mig_DCs",
                                     "Eosinophils","gamma_delta_T","IgA_PC","IgM_IgD_matureB",
                                     "ILC2","ILC3", "Mono_Mac", "Mac1","Mac2", "Mast","Neutrophils"))

##### split umap 
p <- DimPlot(immune, group.by = "annotation", pt.size = 0.5, label = FALSE, label.size = 3,split.by = "condition",
             cols = c("#7F286D", "#270A7F", "#9BC6C9","#26DFED",
                      "#F281CC","#F20AB1","#E22B17","#AB95EF",
                      "#DDED0C","#EDB20C","#4166EF","#6899C1","#54EF0C", 
                      "#415B37","#6D7070","#C0DBB4","#950CED","#EDB4E7"))
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/umap_annotated.svg", width = 8, height = 5, plot = p)

#colors: 
#endothelial: #997F65
#epithelial: #AA8C50
#Fibroblasts: 422D03
#cDC1: #F281CC, cDC2 #F20AB1, activated DCs: #7F286D, pDCs:#EDB4E7
#Eos: #E22B17
#B:mature IghM: #EDB20C; Plasma celle IgA: #DDED0C
#T: CD4 T : #270A7F, T_gamma_delta: #AB95EF, ILC2: #4166EF, ILC3: 6899C1, CD8/NKT: #26DFED , CD8_CCR7: #9BC6C9
#Macs: 1: #54EF0C, mac_mono: #C0DBB4, 2: #415B37
#Mast: #6D7070
#Neutrophils:: #950CED
#eos: #E22B17

##### proportion plot 
create_table_cell_type_prop(immune, "condition","annotation","/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/","neo_immune")
df <- read.csv("/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/neo_immune_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:19)) 

df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c("gamma_delta_T", "CD4_T","CD8_CCR7_T", "CD8T_NKT","ILC2","ILC3","cDC2",
                                                                   "cDC1","act_mig_DCs","pDCs","Mac2","Mac1","Mono_Mac", "IgA_PC","IgM_IgD_matureB",
                                                                   "Mast","Neutrophils","Eosinophils")))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#AB95EF", "#270A7F", "#9BC6C9","#26DFED",
                               "#4166EF","#6899C1","#F20AB1","#F281CC",
                               "#7F286D","#EDB4E7","#415B37","#54EF0C","#C0DBB4", 
                               "#DDED0C","#EDB20C","#6D7070","#950CED","#E22B17")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/cell_type_prop.svg", width = 12, height = 6, plot = p)

##### statistics 
#neo wt vs adult wt 
cell_type_prop_stats(immune,"annotation","adult_wt","neo_P14_wt","condition",1.2,
                     "/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/neo_wt_vs_adult_wt_stat.pdf") 

#neo wt vs phil 
cell_type_prop_stats(immune,"annotation","neo_P14_wt","neo_P14_phil","condition",1.2,
                     "/scratch/khandl/Neonatal_eosinophils/figures/CD45enr_prop/neo_phil_vs_wt_stat.pdf") 
