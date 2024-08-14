#this script takes the AUCs data (calculated in script 1), and returns the correlations of AUCs to gene expression and dependency
#as shown in figure 2

library(tidyverse)
library(corrr)
library(factoextra)
library(Hmisc)
library(openxlsx)
library(cogena)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(VennDiagram)
library(eulerr)
library(GGally)
library(rentrez)
library(gprofiler2)

#additional data required, downloaded from DEPMAP (https://depmap.org/portal/)
#CCLE - download the file "CCLE_expression.csv" (to reproduce the results in the paper, download the 20Q4 version)
#Achilles - download the file "Achilles_gene_effect.csv" (to reproduce the results in the paper, download the 21Q1 version)

#------------------------------------------------------------------------------
# Correlation to CCLE expression and essentiality
#------------------------------------------------------------------------------

#before running this, create "FIN_AUCs" from the "1 - Dose response curves and AUC determination" script
AUCs <- FIN_AUCs %>% 
  select(inducer, cell, mean) %>% 
  spread(key = "inducer", value = "mean") 

Expression_corrMatrix <- AUCs %>% 
  left_join(CCLE, by = c("cell" = "ID")) %>% 
  dplyr::select(-`cell`, -`TNBC status`, -`PAM50 (Ron)`, -`PAM50 (papers)`,
                -`Groups for megacharts`, -`CCLE_Name`, -`Groups for charts`) %>% 
  cor(method = "pearson", use = "complete.obs") %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::select(gene, `Erastin`, `RSL3`, `FIN56`) %>% 
  filter(!(gene %in% c("Erastin", "RSL3", "FIN56"))) %>% 
  na.omit()

Depscores_corrMatrix <- AUCs %>% 
  left_join(Achilles, by = c("cell" = "ID")) %>% 
  select(-`cell`, -`TNBC status`, -`PAM50 (Ron)`, -`PAM50 (papers)`,
         -`Groups for megacharts`, -`Groups for charts`, -lineage) %>% 
  cor(method = "pearson", use = "complete.obs") %>% 
  as_tibble(rownames = "gene") %>% 
  dplyr::select(gene, `Erastin`, `RSL3`, `FIN56`) %>% 
  filter(!(gene %in% c("Erastin", "RSL3", "FIN56"))) %>% 
  na.omit()

Correlations <- full_join(
  Expression_corrMatrix %>% rename_with(.cols = -gene, .fn = ~str_c(.x, " AUC to expression")),
  Depscores_corrMatrix %>% rename_with(.cols = -gene, .fn = ~str_c(.x, " AUC to Depscore")),
  by = "gene")

#write the results to a file, to be used by the pipeline script  
write_csv(Correlations, "Correlation of gene expression and depScore to AUCs (2023).csv")

#------------------------------------------------------------------------------
# Correlations analysis
#------------------------------------------------------------------------------

#generate ranked gene list
RankedList <- Correlations %>% 
  mutate(across(-gene, ~min_rank(-.x))) %>% 
  rename_with(.cols = -gene, .fn = ~str_c(.x, " rank"))

#save excel sheet (this will create supplemental table S1)
thisFig <- left_join(Correlations, RankedList, by = "gene")

wb <- createWorkbook()
title <- "Correlations"
addWorksheet(wb, title, gridLines = FALSE)
writeData(wb, sheet = title, x = "Correlation of AUCs to gene expression and essentiality", startCol = 1, startRow = 1)
writeDataTable(wb, sheet = title, x = thisFig, startCol = 1, startRow = 3, 
               tableStyle = "TableStyleLight1", bandedRows = FALSE)
saveWorkbook(wb, file = "Correlation of AUCs to CCLE expression.xlsx", overwrite = TRUE)
  
#------------------------------------------------------------------------------
# Correlations plots
#------------------------------------------------------------------------------

#Venn diagram of the positivly correlated
Venn1E_list <- list("Erastin" = Expression_corrMatrix$gene[Expression_corrMatrix$Erastin > 0.65],
                    "RSL3" = Expression_corrMatrix$gene[Expression_corrMatrix$RSL3 > 0.65],
                    "FIN56" = Expression_corrMatrix$gene[Expression_corrMatrix$FIN56 > 0.65])
                 
Venn1E_partitions <- get.venn.partitions(Venn1E_list) 

venn.diagram(Venn1E_list, filename = "VENN of positive correlating genes (expression) to inducers AUCs.tiff", fill = c("steelblue3", "gray50", "lightgoldenrod2"),
             alpha = 0.6, cex = 6, cat.cex = 6, cat.dist = 0.07, resolution = 300, lwd = 4,
             fontfamily = "sans", cat.fontfamily = "sans", euler.d = FALSE, scaled = FALSE)

Venn1D_list <- list("Erastin" = Depscores_corrMatrix$gene[Depscores_corrMatrix$Erastin > 0.65],
                    "RSL3" = Depscores_corrMatrix$gene[Depscores_corrMatrix$RSL3 > 0.65],
                    "FIN56" = Depscores_corrMatrix$gene[Depscores_corrMatrix$FIN56 > 0.65])

Venn1D_partitions <- get.venn.partitions(Venn1D_list) 

venn.diagram(Venn1D_list, filename = "VENN of positive correlating genes (DepScores) to inducers AUCs.tiff", fill = c("steelblue3", "gray50", "lightgoldenrod2"),
             alpha = 0.6, cex = 6, cat.cex = 6, cat.dist = 0.07, resolution = 300, lwd = 4,
             fontfamily = "sans", cat.fontfamily = "sans", euler.d = FALSE, scaled = FALSE)


#comparing the negatively correlated
Venn2E_list <- list("Erastin" = Expression_corrMatrix$gene[Expression_corrMatrix$Erastin < -0.65],
                    "RSL3" = Expression_corrMatrix$gene[Expression_corrMatrix$RSL3 < -0.65],
                    "FIN56" = Expression_corrMatrix$gene[Expression_corrMatrix$FIN56 < -0.65])

Venn2E_partitions <- get.venn.partitions(Venn2E_list) 

venn.diagram(Venn2E_list, filename = "VENN of negative correlating genes (expression) to inducers AUCs.tiff", fill = c("steelblue3", "gray50", "lightgoldenrod2"),
             alpha = 0.6, cex = 6, cat.cex = 6, cat.dist = 0.07, resolution = 300, lwd = 4,
             fontfamily = "sans", cat.fontfamily = "sans", euler.d = FALSE, scaled = FALSE)

Venn2D_list <- list("Erastin" = Depscores_corrMatrix$gene[Depscores_corrMatrix$Erastin < -0.65],
                    "RSL3" = Depscores_corrMatrix$gene[Depscores_corrMatrix$RSL3 < -0.65],
                    "FIN56" = Depscores_corrMatrix$gene[Depscores_corrMatrix$FIN56 < -0.65])

Venn2D_partitions <- get.venn.partitions(Venn2D_list) 

venn.diagram(Venn2D_list, filename = "VENN of negative correlating genes (DepScores) to inducers AUCs.tiff", fill = c("steelblue3", "gray50", "lightgoldenrod2"),
             alpha = 0.6, cex = 6, cat.cex = 6, cat.dist = 0.07, resolution = 300, lwd = 4,
             fontfamily = "sans", cat.fontfamily = "sans", euler.d = FALSE, scaled = FALSE)

#------------------------------------------------------------------------------
# plot the genes and correlations, with pubmed terms
#------------------------------------------------------------------------------

#This section performes pubmed searches, which will obviously be different than the pubmed search done to create figure 2 in the paper (November 2023)
#To skip the pubmed searches, run the following line, which will load the relevant results of the pubmed search done in November 2023.
#Two files will be loaded: "consensus_genes_expression_pubmed" and "consensus_genes_depscore_pubmed"
#Than, skip all calls to the "get_pubmed_info" function later in the code
load("Data/Correlation consensus genes pubmed search.RData")

#pubmed function
Pubmed_terms <- c("Ferroptosis", "Lipid peroxidation", "Iron", "Glutathione", "GPX4", "SLC7A11")
set_entrez_key("") #the entrz key for accessing pubmed is private and should be added here by the user running this script

get_pubmed_info <- function(genes, terms) {
  pubmed_search <- list()
  for(i in seq_along(genes)) {
    print(genes[i])
    
    for(j in seq_along(Pubmed_terms)) {
      search_term <- paste0(genes[i], " AND ", Pubmed_terms[j])
      result <- entrez_search(db = "pubmed", term = search_term)
      pubmed_search[[i*10 + j]] <- tibble("gene" = genes[i], "term" = Pubmed_terms[j], "count" = result[["count"]])
    }
    Sys.sleep(0.1)
  }
  pubmed_search_tib <- purrr::reduce(pubmed_search, bind_rows)
}

#plot drawing function
plot_correlations <- function(data) {
  p1 <- ggplot(data, aes(x = gene, y = R, color = FIN)) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.6, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = "", y = "Correlation to\nFIN AUC") +
    theme_bw(base_size = 14) +
    theme(legend.position = c(0.25, 0.2),
          legend.direction = "horizontal",
          legend.background = element_rect(color = "black", size = 0.2),
          panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  p2 <- ggplot(data, aes(x = gene, y = term, fill = count)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "salmon", high = "red", na.value = "white") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    labs(x = "", y = "") +
    theme_bw(base_size = 16) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  ggarrange(p1 + rremove("x.text"), ggplot() + theme_void(), p2, ncol = 1, align = "v", heights = c(0.5, -0.06, 0.6))
}

#extract consensus genes
consensus_genes_expression <- Venn1E_partitions %>% 
  filter(Erastin == TRUE & RSL3 == TRUE & FIN56 == TRUE) %>% 
  select(..values..) %>% 
  unlist()

consensus_genes_depscores <- Venn1D_partitions %>% 
  filter(Erastin == TRUE & RSL3 == TRUE & FIN56 == TRUE) %>% 
  select(..values..) %>% 
  unlist()

#draw pubmed plot for expression
consensus_genes_expression_pubmed <- get_pubmed_info(consensus_genes_expression, Pubmed_terms)

thisFig <- consensus_genes_expression_pubmed %>% 
  left_join(Expression_corrMatrix, by = "gene") %>% 
  arrange(desc(`Erastin`)) %>% 
  mutate(gene = fct_inorder(gene)) %>% 
  gather(-c(gene, term, count), key = "FIN", value = "R") %>% 
  mutate(count = ifelse(count == 0, NA, count))
  
plot_correlations(thisFig)
ggsave("Correlation of gene expression to FIN AUC with pubmed.tiff", dpi = 300, height = 4.5, width = 9)

#draw pubmed plot for depScores
consensus_genes_depscore_pubmed <- get_pubmed_info(consensus_genes_depscores, Pubmed_terms)

thisFig <- consensus_genes_depscore_pubmed %>% 
  left_join(Depscores_corrMatrix, by = "gene") %>% 
  arrange(desc(`Erastin`)) %>% 
  mutate(gene = fct_inorder(gene)) %>% 
  gather(-c(gene, term, count), key = "FIN", value = "R") %>% 
  mutate(count = ifelse(count == 0, NA, count))

plot_correlations(thisFig)
ggsave("Correlation of gene depScores to FIN AUC with pubmed.tiff", dpi = 300, height = 4.5, width = 10)
