##### Signature scores comparing two conditions 
### FcgR mediated phagocytosis - based on MSigDB: KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS
FcgR_phagocytosis_vln_2_cond_mm <- function(
    seurat_object,
    column_oi1,
    ident_to_subset,
    condition_column,
    first_identity_for_stats, 
    second_identity_for_stats,
    colors_to_plot,
    y_limits,
    output_path_file
){
  Idents(seurat_object) <- column_oi1
  sub <- subset(seurat_object, idents = ident_to_subset)
  
  FcgR_genes <- list(c("Akt1"   ,  "Akt2"  ,   "Arf6",     "Arpc1a",   "Arpc1b" ,  "Arpc2",    "Arpc3" ,   "Arpc4",    "Arpc5",    "Arpc5l",   "Asap1" ,   "Asap2",    "Asap3" ,   "Cdc42",    "Cfl1" ,   
                       "Cfl2"  ,   "Crk"    ,  "Crkl" ,    "Dnm1l"  ,  "Gab2"    , "Hck"   ,   "Inpp5d" ,  "Limk2" ,   "Lyn"   ,   "Map2k1" ,  "Mapk1"  ,  "Mapk3" ,   "Marcks" ,  "Myo10" ,   "Ncf1"  ,  
                       "Pak1"  ,   "Pik3ca"  , "Pik3cb",   "Pik3r1"  , "Pik3r2"   ,"Pik3r3" ,  "Pik3r5"  , "Pip5k1c",  "Pla2g4a",  "Plcg2"   , "Pld1"    , "Plpp2"  ,  "Prkcd"   , "Prkce"  ,  "Prkcg"  , 
                       "Rac1"    , "Rac2"   ,  "Raf1"   ,  "Rps6kb1"  ,"Rps6kb2",  "Sphk2"   , "Syk"  ,    "Vasp"   ,  "Vav2" ,    "Wasf2"  ,  "Wasl" ,    "Dnm2"   ,  "Fcgr2b",   "Gsn"     , "Pik3cg"  ,
                       "Pikfyve" , "Pip5k1a" , "Pip5k1b" , "Plpp1",    "Prkca"   , "Ptprc"  ,  "Dnm1"  ,   "Dock2"   , "Fcgr4" ,   "Marcksl1", "Pik3cd",   "Pip4k2b" , "Prkcb"  ,  "Vav1" ,    "Vav3"    ,
                       "Pla2g6" ,  "Sphk1"  ,  "Pld2"     ,"Akt3"  ,   "Limk1"    ,"Was"     , "Plpp3"  ,  "Plcg1"    ,"Lat"    ,  "Dnm3"   ,  "Fcgr1"  ,  "Pla2g4d"  ,"Wasf1"   , "Scin"  ,   "Pla2g4f" ,
                       "Amph"  ,   "Wasf3"   , "Pla2g4e" ))
  
  sub <-AddModuleScore(sub, features= FcgR_genes,name = "FcgR_score")
  names(x = sub[[]])
  
  #test
  Idents(sub) <- condition_column
  first <- subset(sub, idents = first_identity_for_stats)
  second <- subset(sub, idents = second_identity_for_stats)
  p_val <- wilcox.test(first$FcgR_score1, second$FcgR_score1, alternative = "two.sided")
  p_val <- p_val$p.value
  
  p <- VlnPlot(sub, features= "FcgR_score1", group.by = condition_column, pt.size = 0, cols = colors_to_plot) +  theme_classic() + 
    theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
    labs(title = paste0("cluster", ident_to_subset), y = "FcgR phago score", x="") + theme(legend.position="right") +  annotate("text",x=1,y=0.5, label = paste0("p_val wilcox = ",p_val)) +
    annotate("text",x=1,y=0.5, label = paste0("p_val wilcox = ",p_val)) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
    scale_y_continuous(limits = y_limits)
  print(p)
  ggsave(output_path_file, width = 8, height = 8, plot = p)
}

### Granulogenesis score - Farifax et al 2018 and Gurtner et al 2023
granulogenesis_vln_2_cond_mm <- function(
    seurat_object,
    column_oi1,
    ident_to_subset,
    condition_column,
    first_identity_for_stats, 
    second_identity_for_stats,
    colors_to_plot,
    y_limits,
    output_path_file
){
  Idents(seurat_object) <- column_oi1
  sub <- subset(seurat_object, idents = ident_to_subset)
  
  Granules_synthesis_list <- list(c("Prg2","Prg3",  "Epx", "Ear6", "Ear1", "Ear2"))
  
  sub <-AddModuleScore(sub, features= Granules_synthesis_list,name = "GranulesSynthesis")
  names(x = sub[[]])
  
  #test
  Idents(sub) <- condition_column
  first <- subset(sub, idents = first_identity_for_stats)
  second <- subset(sub, idents = second_identity_for_stats)
  p_val <- wilcox.test(first$GranulesSynthesis1, second$GranulesSynthesis1, alternative = "two.sided")
  p_val <- p_val$p.value
  
  p <- VlnPlot(sub, features= "GranulesSynthesis1", group.by = condition_column, pt.size = 0, cols = colors_to_plot) +  
    theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
    labs(title = paste0("cluster", ident_to_subset), y = " Granulogenesis score", x="") + theme(legend.position="right") +  
    annotate("text",x=1,y=0.5, label = paste0("p_val wilcox = ",p_val)) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
    scale_y_continuous(limits = y_limits)
  print(p)
  ggsave(output_path_file, width = 8, height = 8, plot = p)
}

### ROS - Maas et al. 2023
ROS_vln_2_cond_mm <- function(
    seurat_object,
    column_oi1,
    ident_to_subset,
    condition_column,
    first_identity_for_stats, 
    second_identity_for_stats,
    colors_to_plot,
    y_limits,
    output_path_file
){
  Idents(seurat_object) <- column_oi1
  sub <- subset(seurat_object, idents = ident_to_subset)
  
  ROS_list <- list(c("Atox1", "Cat", "Ccs", "Glrx", "Ipcef1", "Ncf1", "Ncf2", "Ncf4"))
  
  sub <-AddModuleScore(sub, features= ROS_list,name = "ROS")
  names(x = sub[[]])
  
  #test
  Idents(sub) <- condition_column
  first <- subset(sub, idents = first_identity_for_stats)
  second <- subset(sub, idents = second_identity_for_stats)
  p_val <- wilcox.test(first$ROS1, second$ROS1, alternative = "two.sided",)
  p_val <- p_val$p.value
  
  p <- VlnPlot(sub, features= "ROS1", group.by = condition_column, pt.size = 0, cols = colors_to_plot) +  
    theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
    labs(title = paste0("cluster", ident_to_subset), y = " ROS", x="") + theme(legend.position="right") +  
    annotate("text",x=1,y=0.5, label = paste0("p_val wilcox = ",p_val)) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
    scale_y_continuous(limits = y_limits)
  print(p)
  ggsave(output_path_file, width = 8, height = 8, plot = p)
}

### Glycolytic activity score - Mi et al 2013
Glycolytic_vln_2_cond_mm <- function(
    seurat_object,
    column_oi1,
    ident_to_subset,
    condition_column,
    first_identity_for_stats, 
    second_identity_for_stats,
    colors_to_plot,
    y_limits,
    output_path_file
){
  Idents(seurat_object) <- column_oi1
  sub <- subset(seurat_object, idents = ident_to_subset)
  
  Glycolytic_synthesis_list <- list(c("Gapdh","Pgk1",  "Pgam1", "Tpi1", "Aldoa", "Ldha","Pkm","Eno1", "Hk1","Hk2"))
  
  sub <-AddModuleScore(sub, features= Glycolytic_synthesis_list,name = "GlycolyticActivity")
  names(x = sub[[]])
  
  #test
  Idents(sub) <- condition_column
  first <- subset(sub, idents = first_identity_for_stats)
  second <- subset(sub, idents = second_identity_for_stats)
  p_val <- wilcox.test(first$GlycolyticActivity1, second$GlycolyticActivity1, alternative = "two.sided")
  p_val <- p_val$p.value
  
  p <- VlnPlot(sub, features= "GlycolyticActivity1", group.by = condition_column, pt.size = 0, cols = colors_to_plot) +  
    theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
    labs(title = paste0("cluster", ident_to_subset), y = " GlycolyticActivity score", x="") + theme(legend.position="right") +  
    annotate("text",x=1,y=0.5, label = paste0("p_val wilcox = ",p_val)) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
    scale_y_continuous(limits = y_limits)
  print(p)
  ggsave(output_path_file, width = 8, height = 8, plot = p)
}


