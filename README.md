## This repository contains the code that has been used to analyse bulk- and scRNA-seq of the Arnold lab's project to investigate eosinophils in neonatal biology 

### Data used to reproduce the code 
#### Newly acquired bulk- and scRNA-seq data can be found under GEO accession number ... 
#### Public data from Gurtner et al, Nature 2022, adult Il5-tg eos from various tissues, GEO accession number GSE182001 and processed Seurat objects in zenodo (https://zenodo.org/records/10225349)

### Code description: 

**1. Packages and functions:**<br> 
Configuration of input and output file paths, loading of required packages to run this repository, generation of the results folder structure and functions needed to run pre-processing, cell type annotation, DGE analysis, cell type proportional differences, signature score analysis, GSEA, PROGENy, SCENIC, CellPhoneDB and FACS data plotting 

**2. Data integration into the Seurat workflow:**<br>
Integration of raw cell-gene count matrices generated from the BD WTA Analysis pipeline on the Seven Bridges server into the Seurat workflow 

**3. Cell type broad annotation:**<br>
Cell type annotation of all datasets using SingleR, marker gene expression and DEG analysis 

**4. Integration of eosinophils from different datasets:**<br>
Integration of eosinophils from different datasets and tissues from human and mouse using either label transfer to the steady state reference from Gurtner et al 2022, Nature, or FastMNN integration 

**5. Analysis of cell type proportion differences between conditions:**<br>
Cell type proportional comparisons using the scProportionTest R package between different conditions and tissues 

**6. BM derived precursors (GMP - granulocyte-monocyte progenitors and MPP - multipotent progenitors):**<br>
Integration, annotation, DGE analysis and granule protein genes expression analysis of bone marrow progenitors from adults vs. neonates 

**Pseudotime trajectory analysis between adult and neo P14 Il5-Tg BM, blood, CO and SI:**<br>
Application of slingshot to infer pseudotime, and tradeSeq/condiment to analyse DEGs between age along the trajectory

**8. DEG analysis between conditions:**<br>
Differential gene expression analysis between conditions and tissues across cell types applying a Wilcoxon Rank sum test for mouse data and DESeq2 on pseudobulks for human data (multiple samples per condition available)

**9. Signature scores:**<br>
Signature score calculation of FcgR-mediated phagocytosis, granulogenesis, ROS and Glycolytic activity scores between neo and adult, phil and wt, in eosinophils and other cell types 

**10. Signalling pathway analyses between Il5-Tg CO and SI eos adult vs. Neo P14:**<br>
Application of GSEA, PROGENy, SCENIC and CellPhoneDB between different conditions and cell types 

**11. Comparison of gene expression in eos between conditions and datasets:**<br>
Comparison of gene expression in terms of significant DEGs and average expression across different conditions 

**12.Bulk-RNAseq data analysis:**<br>
Bulk-RNAseq DGE analysis of in-vitro differentiated eosinophils from neonates and adults

**13.FACS data plotting:**<br>
Plotting of FACS MFI and percentages acquired from eosinophils from adults and neonates across tissues
