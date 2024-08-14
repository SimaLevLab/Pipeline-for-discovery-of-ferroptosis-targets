#this script analyzes the count matrix of the RNAseq used in this paper (for figure 3)

library(tidyverse)
library(factoextra)
library(ggrepel)
library(limma)
library(edgeR)
library(dendextend)
library(ggforce)
library(RColorBrewer)
library(viridis)

DEGCreate <- function (cts, name, plot = FALSE) {
  dge <- DGEList(counts = cts, genes = rownames(cts))
  groups <- colnames(dge)
  groups <- as.factor(str_sub(groups, 1, str_length(groups)-2))
  dge$samples$group <- groups
  library(org.Hs.eg.db)
    egSYMBOL <- toTable(org.Hs.egSYMBOL)
    m <- match(dge$genes$genes, egSYMBOL$symbol)
    dge$genes$gene_id <- egSYMBOL$gene_id[m]
  detach(package:org.Hs.eg.db)
  detach(package:AnnotationDbi)
  keep <- filterByExpr(dge, group = groups)
  dge <- dge[keep,, keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  CPM <- cpm(dge)
  LogCPM <- cpm(dge, log=TRUE)
  if(plot) {
    par(mfrow=c(1,2))
    plotMDS(LogCPM, labels=groups)
    title(main="MDS dim 1 & 2")
  }
  designMatrix <- model.matrix(~0 + groups)
  colnames(designMatrix) <- gsub("groups", "", colnames(designMatrix))
  V <- voom(dge, designMatrix, plot = plot)
}

#1st RNAseq
cts1 <- read.csv("Data/RNAseq set 1 - counts.csv", row.names = "Gene")
V1 <- DEGCreate(cts1, name = "Erastin and RSL3 1st set")

#2nd RNAseq
cts2 <- read.csv("Data/RNAseq set 2 - counts.csv", row.names = "Gene")
V2 <- DEGCreate(cts2, name = "Erastin and RSL3 2nd set")

#LIMMA analysis of set 1
design1 <- V1$design
vfit1 <- lmFit(V1, design1)
contr.matrix1 <- makeContrasts(
  BT549_E = BT549_E - BT549_C,
  BT549_R = BT549_R - BT549_C,
  HS578_E = HS578_E - HS578_C,
  HS578_R = HS578_R - HS578_C,
  MDA468_E = MDA468_E - MDA468_C,
  MDA468_R = MDA468_R - MDA468_C,
  HCC70_E = HCC70_E - HCC70_C,
  HCC70_R = HCC70_R - HCC70_C,
  levels = colnames(design1))

vfit1c <- contrasts.fit(vfit1, contr.matrix1)
efit1 <- eBayes(vfit1c, 0.01)
tT1 <- topTable(efit1, adjust.method = "fdr", sort.by="B", n = Inf)

decide1 <- decideTests(efit1, adjust.method = "none") 
decide1 <- as.data.frame(decide1)
decide1 <- as_tibble(decide1, rownames = "genes") %>% 
  rename_with(~paste(.x, "_pval", sep = ""), .cols = -genes)

DEG1 <- tT1 %>% left_join(decide1, by = "genes")

#LIMMA analysis of set 2
design2 <- V2$design
vfit2 <- lmFit(V2, design2)
contr.matrix2 <- makeContrasts(
  SUM159_E = SUM159_E - SUM159_C,
  SUM159_R = SUM159_R - SUM159_C,
  levels = colnames(design2))

vfit2c <- contrasts.fit(vfit2, contr.matrix2)
efit2 <- eBayes(vfit2c, 0.01)
tT2 <- topTable(efit2, adjust.method = "fdr", sort.by="B", n = Inf)

decide2 <- decideTests(efit2, adjust.method = "none") 
decide2 <- as.data.frame(decide2)
decide2 <- as_tibble(decide2, rownames = "genes") %>% 
  rename_with(~paste(.x, "_pval", sep = ""), .cols = -genes)

DEG2 <- tT2 %>% left_join(decide2, by = "genes")

#combine the two DEG files
DEG_combined <- DEG1 %>% full_join(DEG2, by = "genes") %>% 
  dplyr::select(genes, matches("_E$|R$"), matches("pval$")) 

#pavals version
tT1_pvals <- efit1[["p.value"]] %>% 
  as_tibble(rownames = "genes") %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_pvalue")) %>% 
  left_join(tT1, by = "genes")

tT2_pvals <- efit2[["p.value"]] %>% 
  as_tibble(rownames = "genes") %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_pvalue")) %>% 
  left_join(tT2, by = "genes")

DEG_combined_v2 <- full_join(tT1_pvals, tT2_pvals, by = "genes") %>% 
  dplyr::select(genes, matches("_E$|R$"), matches("pvalue$"))

#t-statistics version
tT1_tstat <- efit1[["t"]] %>% 
  as_tibble(rownames = "genes") %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_t")) %>% 
  left_join(tT1, by = "genes")

tT2_tstat <- efit2[["t"]] %>% 
  as_tibble(rownames = "genes") %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_t")) %>% 
  left_join(tT2, by = "genes")

DEG_combined_v3 <- full_join(tT1_tstat, tT2_tstat, by = "genes") %>% 
  dplyr::select(genes, matches("_E$|R$"), matches("t$"))
