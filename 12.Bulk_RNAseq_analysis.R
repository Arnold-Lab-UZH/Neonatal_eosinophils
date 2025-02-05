########## This code analyses bulk RNAseq data ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Neonatal_eosinophils/1.1.Packages.R")

##### read in all samples and combine them in one data frame 
AD1 <- read.delim("/data/khandl/eos_NEO/bulkRNAseq/1.neo_project/AD1.txt")
AD2 <- read.delim("/data/khandl/eos_NEO/bulkRNAseq/1.neo_project/AD2.txt")
AD3 <- read.delim("/data/khandl/eos_NEO/bulkRNAseq/1.neo_project/AD3.txt")
ADBLD <- read.delim("/data/khandl/eos_NEO/bulkRNAseq/1.neo_project/ADBLD_Replacement_1.txt")
NEO1 <- read.delim("/data/khandl/eos_NEO/bulkRNAseq/1.neo_project/NEO1.txt")
NEO2 <- read.delim("/data/khandl/eos_NEO/bulkRNAseq/1.neo_project/NEO2.txt")
NEO3 <- read.delim("/data/khandl/eos_NEO/bulkRNAseq/1.neo_project/NEO3.txt")

# extract target_ids and TPM 
AD1 <- AD1[,colnames(AD1) %in% c("target_id","tpm")]
colnames(AD1) <- c("target_id","AD1")
AD2 <- AD2[,colnames(AD2) %in% c("target_id","tpm")]
colnames(AD2) <- c("target_id","AD2")
AD3 <- AD3[,colnames(AD3) %in% c("target_id","tpm")]
colnames(AD3) <- c("target_id","AD3")
ADBLD <- ADBLD[,colnames(ADBLD) %in% c("target_id","tpm")]
colnames(ADBLD) <- c("target_id","ADBLD")
NEO1 <- NEO1[,colnames(NEO1) %in% c("target_id","tpm")]
colnames(NEO1) <- c("target_id","NEO1")
NEO2 <- NEO2[,colnames(NEO2) %in% c("target_id","tpm")]
colnames(NEO2) <- c("target_id","NEO2")
NEO3 <- NEO3[,colnames(NEO3) %in% c("target_id","tpm")]
colnames(NEO3) <- c("target_id","NEO3")

merged <- AD1 %>% right_join(AD2, by=c("target_id"))
merged <- merged %>% right_join(AD3, by=c("target_id"))
merged <- merged %>% right_join(ADBLD, by=c("target_id"))
merged <- merged %>% right_join(NEO1, by=c("target_id"))
merged <- merged %>% right_join(NEO2, by=c("target_id"))
merged <- merged %>% right_join(NEO3, by=c("target_id"))

##### convert ensemble to gene IDs 
mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
attributes<-listAttributes(mart)
gene_ids <- getBM(attributes = c("ensembl_transcript_id","external_gene_name"), mart = mart)

colnames(merged) <- c("ensembl_transcript_id","AD1","AD2","AD3","ADBLD","NEO1","NEO2","NEO3")

df<-merged%>%
  left_join(dplyr::select(gene_ids,1:2))
#remove duplicates 
df<-df[!duplicated(df$external_gene_name),]
df[is.na(df)]<-0 #make all empty value to zero
#add external gene names as rownames (instead of ensemble IDs before)
rownames(df)<-df$external_gene_name
#remove unnessesary columns 
df<-dplyr::select(df,-ensembl_transcript_id,-external_gene_name)

##### DEG analysis between AD and NEO using edgeR 
#remove ADBLD 
df2 <- df[,colnames(df) %in% c("AD1","AD2","AD3","NEO1","NEO2","NEO3")]

#Create DGEList object
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=df2,group=group)
y$samples
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2, main = "Barplot of library sizes")

# Filter reads by counts: Most of the samples should have at least 10 reads, normalize the library and estimate dispersion
keep <- filterByExpr(y, min.count = 10)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~group)
y <- estimateDisp(y,design)
plotMDS(y, pch = 2, label=colnames(y)) 
logcpm <- cpm(y, log=TRUE)
head(logcpm)

# NEO over AD 
#so if the pair is c("A","B") then the comparison is B - A, 
#so genes with positive log-fold change are up-regulated in group B 
#compared with group A (and vice versa for genes with negative log-fold change).
neo_vs_ad <- exactTest(y, pair=c(1,2)) #the first group listed in the pair is the baseline for the comparison
topTags(neo_vs_ad)
summary(decideTests(neo_vs_ad))

neo_vs_ad = topTags(neo_vs_ad, n = Inf)
dim(neo_vs_ad)
head(neo_vs_ad$table)
neo_vs_ad$table$Gene <- rownames(neo_vs_ad$table)
View(neo_vs_ad$table)

write.table(neo_vs_ad$table, file = "/scratch/khandl/Neonatal_eosinophils/bulk_RNAseq/NEO_vs_AD.txt")

##### plot GOI of bulk data 
# read in excel file 
df <- neo_vs_ad$table

df <- df %>% mutate(
  Expression = case_when(logFC >= 0.25 & FDR <= 0.05 ~ "neo",
                         logFC <= -0.25 & FDR <= 0.05 ~ "adult",
                         TRUE ~ "Non sig."))

df <- as.data.frame(df)

rownames(df) <- df$Gene

genes_of_interest <- c("Mmp25","Ccr3","Cd101","Itgax","Mmp27","Adam19","S100a9","Camp","Csf3r","Ccl6","Lyz2",
                       "S100a8","Ahr","Csf2rb2","Cd9","Cd274","Mmp9","Alox5ap","Itgam","Mmp8","Alox15","Siglecf",
                       "Epx","Prg3","Prg2","Thbs1","Ear1","Ear2","Ear6")

df$genelabels <- ifelse(rownames(df) %in% genes_of_interest,TRUE,FALSE)

p <- ggplot(df, aes(x=logFC, y=-log10(FDR))) +
  geom_point(size=6, aes(color = Expression)) +
  geom_text_repel(aes(label=ifelse(df$genelabels, rownames(df),"")),size=5, max.overlaps = 1000) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + 
  scale_color_manual(values =c("#EA94A6","#EA94A6","#7F7C7D") ,guide = "none") + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  theme(legend.title = element_text(size = 25), legend.text = element_text(size = 25)) + 
  geom_vline(xintercept = c(-0.25,0.25) , linetype = "dashed", color = "gray50") + #log Fold Change cutoff 
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "gray50") + #p value cutoff 
  theme_classic(base_size = 25) 
print(p)
ggsave("/scratch/khandl/Neonatal_eosinophils/bulk_RNAseq/NEO_vs_AD.svg", width = 12, height = 10, plot = p)

