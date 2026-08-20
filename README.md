#### This repository contains the code that has been used to analyse bulk- and scRNA-seq of the following publication …. 

#### Data used to reproduce the code 
### Newly acquired bulk- and scRNA-seq data can be found under GEO accession number ... 
### Public data from Gurtner et al, Nature 2022, adult Il5-tg eos from various tissues, GEO accession number GSE182001 and processed Seurat objects in zenodo (https://zenodo.org/records/10225349)

### Code description: 

**1. Packages and functions:**
- 1.1.config.R: Configuration of input and output file paths 
- 1.2.Packages.R: All packages that need to be loaded
- 1.3.Output_directory_output_folder_structure_generation.R: Generation of output folder structure needed for this repository  
- 1.4.Functions_preprocessing.R: Functions used for integration of data into the Seurat workflow
- 1.5.Functions_annotation.R: Functions used for cell type annotation 
- 1.6.Functions_DEGs.R: Functions used for DGE analysis (pseudobulk and single-cell approach) and plotting of genes of interest in heatmaps and volcano plots 
- 1.7.Functions_cell_type_prop.R: Functions to compare cell type proportions between conditions in a barplot and statistical analysis 
- 1.8.Functions_signature_scores.R: Functions to apply different gene signature scores and analyse their different expression between conditions 
- 1.9.Functions_GSEA_PROGENy_SCENIC.R: Functions to predict signalling pathways betweenconditions with GSEA, PROGENy and SCENIC
- 1.10.Functions_CellPhoneDB.R: Functions to generate input files for CellPhoneDB analysis 
- 1.11.Functions_FACS_data_plotting.R: Functions to plot FACS data in heatmaps 

**2. Data integration into the Seurat workflow:**
- 2.Preprocessing.R

**3. Cell type broad annotation:**
- 3.1.Anno_eos_P14_CO_SI.R: Eosinophil enriched, Il5-Transgenic (Il5-Tg), neonatal postnatal day 14 (neo P14), colon (CO) and small intestine (SI) 
- 3.2.Anno_eos_P14_BM_blood_spleen.R: Eosinophil enriched, Il5-Tg, neo P14, bone marrow (BM), blood and spleen
- 3.3.Anno_CD45pos_adult_P14_CO_WT_PHIL.R: CD45 enriched, CO, adult wild-type (WT), neo P14 WT and PHIL (eosinophils depleted systemically) 
- 3.4.Anno_CD45neg_P14_CO_WT_PHIL.R: CD45 negative, CO, neo P14, WT and PHIL
- 3.5.Anno_CD45pos_eos_depletion_P16_CO.R: CD45 enriched, CO, eosinophils depleted at day 5,7, and 10 with iDT in Cre+, Cre- used as control, P16
- 3.6.Anno_CB_and_PB.R: Eosinophil enriched from human cord and peripheral blood (CB and PB) 

**4. Integration of eosinophils from different datasets:**
- 4.1.Label_transfer_eos_all_tissues_adult_to_P14.R: Label transfer to transfer eosinophil subtype labels from adult Il5-Tg (Gurtner et al., 2023) to neo P14 Il5-Tg BM, blood, CO and SI eosinophils
- 4.2.Label_transfer_eos_CO_SI_adult_to_P14.R: Label transfer to transfer eosinophil subtype labels from adult Il5-Tg (Gurtner et al., 2022) to neo P14 Il5-Tg CO and SI eosinophils
- 4.3.FastMNN_integration_CB_PB_eos.R: FastMNN integration of CB and PB eosinophils and annotation of subtypes
- 4.4.neo_P14_separate_clustering.R: Integration and new clustering of neo P14 and adult Il5-tg, bm, blood, spleen, colon and small intestine 

**5. Analysis of cell type proportion differences between conditions:**
- 5.1.Eos_type_prop_eos_CO_SI_adult_vs_P14.R: Eosinophil subtypes between Il5-Tg adult and neo P14 CO and SI
- 5.2.Eos_type_prop_eos_CB_PB.R: Eosinophil subtypes between CB and PB 
- 5.3.Cell_type_prop_CD45pos.R: Colon CD45 positive neo WT and PHIL 
- 5.4.Cell_type_prop_CD45neg.R: Colon CD45 negative neo WT and PHIL 

**6. BM derived precursors (GMP - granulocyte-monocyte progenitors and MPP - multipotent progenitors):**
- 6.BM_precursors_adult_vs_P14.R: Il5-Tg adult and neo P14 - integration, annotation, DGE analysis, Granule protein genes expression analysis

**Pseudotime trajectory analysis between adult and neo P14 Il5-Tg BM, blood, CO and SI:**
- 7.Pseudotime_trajectory_eos_all_tissues_adult_vs_P14.R: slingshot to infer pseudo time, tradeSeq and condiment to analyse DEGs between ages along the trajectory

**8. DEG analysis between conditions: **
- 8.1.DGE_analysis_eos_BM_adult_vs_P14.R: BM Eosinophils Il5-Tg adult vs. Neo P14, Wilcox test
- 8.2.DGE_analysis_eos_CO_SI_adult_vs_P14.R: CO and SI Eosinophils Il5-Tg adult vs. Neo P14, Wilcox test, Dotplot for Fcgr related genes
- 8.3.DGE_analysis_eos_CB_vs_PB.R: B-Eos-like CB vs. PB, DESeq2 on pseudobulks
- 8.4.DGE_analysis_CD45pos_CO_P14_PHIL_vs_WT.R: CO neo P14 PHIL vs. WT CD45 positive (Mono/Mac, Mac1, Neutrophils, cDC1, cDC2)
- 8.5.DGE_analysis_CD45neg_CO_P14_PHIL_vs_WT.R: CO neo P14 PHIL vs. WT CD45 negative (Fibroblasts, Pericytes, SMCs) and Ccl11 expression

**9. Signature scores:**
- 9.Signature_scores_eos_BM_blood_CO_adult_vs_P14.R: FcgR-mediated phagocytosis, granulogenesis, ROS and Glycolytic activity scores between neo and adult, phil and wt 

**10. Signalling pathway analyses between Il5-Tg CO and SI eos adult vs. Neo P14:**
- 10.1.Singleseqgset_eos_CO_SI_adult_vs_P14.R: GO term analysis using the singleseqgset package 
- 10.2.PROGENy_eos_CO_SI_adult_vs_P14.R: PROGENy analysis 
- 10.3.SCENIC_eos_CO_SI_adult_vs_P14.R: SCENIC analysis
- 10.4.1.CellPhoneDB_input_file_generation.R: Input file generation of CD45 neg + eos wt and CD45 neg phil 
- 10.4.2.CellPhoneDB_run_CD45neg_plus_eos.sh: Runs CellPhoneDB python script, CD45neg + eos wt
- 10.4.3.CellPhoneDB_run_CD45neg_phil.sh: Runs CellPhoneDB python script, CD45neg phil
- 10.4.4.CellPhoneDB_eos_CD45neg_significant_interactions.R: Generation of one csv file with only the significant interactions between eos and CD45 neg 
- 10.4.5.CellPhoneDB_CD45neg_phil_vs_wt.R: For each interaction pair one csv file is generated to compare phil and wt 
- 10.4.6.CellPhoneDB_CD45neg_eos_interaction_counts.R: Barplots are generated to compare counts of interaction between eos and CD45neg in both directionalities 
- 10.4.7.CellPhoneDB_CD45neg_eos_plotting.R: Interactions of interest are plotted in dot plots between eos and CD45 neg 
- 10.4.8.CellPhoneDB_CD45neg_phil_vs_wt_plotting.R: Interactions of interest are plotted in dot plots between CD45 neg cells comparing phil and wt 

**11. Comparison of gene expression in eos between conditions and datasets:**
- 11.1.Comparison_gene_exp_adult_neo_P14_CO_vs_SI.R: CO vs. SI Il5-Tg eosinophils adult and Neo P14, number of genes expressed
- 11.2.Comparison_gene_exp_eos_CO_Il5tg_vs_wt.R: CO Il5-Tg vs. Wt eosinophils, Label transfer, number of genes expressed, DEGs 
- 11.3.Comparison_gene_exp_eos_CO_P16_depleted_vs_ctrl.R: CO eosinophils neo P16 after eos depletion to control, number of genes expressed, linear correlation graph

**12.Bulk-RNAseq data analysis:**
- 12.Bulk_RNAseq_analysis.R: DGE analysis, in-vitro differentiated eosinophils neo P14 vs. Adult

**13.FACS data plotting:**
- 13.FACS_data_plotting.R: plotting in heatmap of FACS data between neo and adult across tissues 
