#### This repository contains the code that has been used to analyse bulk- and sc-RNA-seq of the following publication …. 
#### The processed files to reproduce this code can be found under GEO accession number: …

#### File contents: 

1. R libraries and functions 
2. Integrating processed count matrices into the Seurat workflow by applying basic quality cutoffs 
3. Annotation of various datasets
    1. Eosinophil enriched, Il5-Transgenic (Il5-Tg), neonatal postnatal day 14 (neo P14), colon (CO) and small intestine (SI)
    2. Eosinophil enriched, Il5-Tg, neo P14, bone marrow (BM), blood and spleen
    3. CD45 enriched, CO, adult wild-type (WT), neo P14 WT and PHIL (eosinophils depleted systemically)
    4. CD45 depleted, CO, neo P14, WT and PHIL 
    5. CD45 enriched, CO, eosinophils depleted at day 5,7, and 10 with iDT in Cre+, Cre- used as control 
    6. Eosinophil enriched from human cord and peripheral blood (CB and PB)
4. Integration of different datasets
    1. Label transfer to transfer eosinophil subtype labels from adult Il5-Tg (Gurtner et al., 2023) to neo P14 Il5-Tg BM, blood, CO and SI eosinophils 
    2. Label transfer to transfer eosinophil subtype labels from adult Il5-Tg (Gurtner et al., 2023) to neo P14 Il5-Tg CO and SI eosinophils (Figure 3)
    3. FastMNN integration of CB and PB eosinophils and annotation of subtypes 
5. Analysis of cell type proportion differences between conditions 
    1. Eosinophil subtypes between Il5-Tg adult and neo P14 CO and SI (Figure 3)
    2. Eosinophil subtypes between CB and PB (Figure 2)
    3. CD45 positive CO samples, neo P14 WT and PHIL, adult WT (Figure 5, S5)
    4. CD45 negative CO samples, neo P14 WT and PHIL (Figure 5, S5)
6. BM derived precursors (GMP - granulocyte-monocyte progenitors and MPP - multipotent progenitors) from Il5-Tg adult and neo P14 - integration, annotation, DGE analysis, Granule protein genes expression analysis (Figure S2) 
7. Pseudotime trajectory analysis between adult and neo P14 Il5-Tg BM, blood, CO and SI; slingshot to infer pseudo time, tradeSeq and condiment to analyse DEGs between ages along the trajectory (Figure 2)
8. DEG analysis between conditions 
    1. BM Eosinophils Il5-Tg adult vs. Neo P14, Wilcox test, Dotplot for Egfr gene expression (Figure S2, 4)
    2. CO and SI Eosinophils Il5-Tg adult vs. Neo P14, Wilcox test, Dotplot for Fcgr related genes (Figure 3, S3, 4)
    3. B-Eos-like CB vs. PB, DESeq2 on pseudobulbs (Figure 2)
    4. CO neo P14 PHIL vs. WT CD45 positive (Mono/Mac, Mac1, Neutrophils, cDC1, cDC2) (Figure 5) 
    5. CO neo P14 PHIL vs. WT CD45 negative (Fibroblasts, Pericytes, SMCs) (Figure 5)
9. Signature scores: FcgR-mediated phagocytosis in Il5-Tg eosinophils from BM, blood and CO adult vs. Neo P14; Granulogenesis score in Il5-Tg eosinophils from BM and blood adult vs. Neo P14; ROS in Il5-Tg eosinophils from CO adult vs. Neo P14; Glycolytic activity score in mono/Mac, mature mac1, cDC1, cDC2 neo P14 WT vs. PHIL (Figure 2,4,5) 
10. Signalling pathway analyses between Il5-Tg CO and SI eos adult vs. Neo P14 
    1.  GO term analysis using the singleseqgset package (Figure 3)
    2. PROGENy analyis (Figure S3)
    3. SCENIC analysis (Figure S3)
11. Comparison of gene expression between eosinophils
    1. CO vs. SI Il5-Tg eosinophils adult and Neo P14, number of genes expressed (Figure S3)
    2. CO Il5-Tg vs. Wt eosinophils, Label transfer, number of genes expressed, DEGs (Figure S5)
    3. CO eosinophils neo P16 after eos depletion to control, number of genes expressed, linear correlation graph (Figure 6, S6)
12. Bulk-RNAseq DGE analysis, in-vitro differentiated eosinophils neo P14 vs. Adult 
