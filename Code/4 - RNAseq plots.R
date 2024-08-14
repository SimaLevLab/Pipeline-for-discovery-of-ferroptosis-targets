#this script will create some of the RNAseq plots shown in figure 3

library(tidyverse)
library(ggpubr)
library(cogena)
library(VennDiagram)
library(ggforce)
library(factoextra)
library(ComplexUpset)
library(rstatix)
library(factoextra)
library(GSVA)
library(cmapR)
library(fgsea)
library(ggtext)
library(glue)
library(gprofiler2)
library(openxlsx)

FerrDB_driver <- read_csv("Data/FerrDB_info_driver.csv") %>% pull(Symbol)
FerrDB_suppresor <- read_csv("Data/FerrDB_info_suppressor.csv") %>% pull(Symbol)
FerrDB_marker <- read_csv("Data/FerrDB_info_marker.csv") %>% pull(Symbol)
validated_biomarkers <-  c("DDIT4", "JDP2", "SLC38A2", "HMOX1", "VLDLR", "SLC6A9", "CHAC1", "MT2A", "OSGIN1", "PCK2",
                           "RPS6KA2", "SLC1A5", "ARHGEF2", "CALM2", "GARS", "PHGDH", "SPIN4", "TMEM154", "SACS", "NDC80",
                           "ERMAP", "ASNS", "PRR11", "MYC", "SLC7A11", "TXNRD1")

#run the RNAseq analysis script first
source("3 - RNAseq analysis.R")

#-------------------------------------------------------------------------------------
# generate excel list (this will create supplemental table S2)
#-------------------------------------------------------------------------------------

#output to excel file
tempDF <- DEG_combined_v2 %>% 
  mutate(mean = rowMeans(across(-c(genes, contains("pvalue"))), na.rm = TRUE)) %>% 
  arrange(desc(mean)) %>% 
  select(-mean)

wb <- createWorkbook()
options("openxlsx.numFmt" = "0.00")
addWorksheet(wb, "DEGs", gridLines = FALSE)
writeData(wb, sheet = "DEGs", x = "Log2FC of Erastin and RSL3 vs. DMSO", startCol = 1, startRow = 1)
writeDataTable(wb, sheet = "DEGs", x = tempDF, startCol = 1, startRow = 3, 
               tableStyle = "TableStyleLight1", bandedRows = FALSE)
negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
NAstyle <- createStyle(bgFill = "gray")
conditionalFormatting(wb, sheet = "DEGs", cols = 2:11, rows = 4:14518, type = "expression", 
                      rule = "AND(B4 > 0, L4 < 0.05)", style = posStyle)
conditionalFormatting(wb, sheet = "DEGs", cols = 2:11, rows = 4:14518, type = "expression", 
                      rule = "AND(B4 < 0, L4 < 0.05)", style = negStyle)

saveWorkbook(wb, file = "DEGs of RNAseq.xlsx", overwrite = TRUE)

#-------------------------------------------------------------------------------------
# PCA plots
#-------------------------------------------------------------------------------------

#combine the two RNAseq - just for the PCA plots 
cts_combined <- bind_cols(cts1, cts2) %>% 
  select(-matches("MDA231|HCC38"))
V_combined <- DEGCreate(cts_combined, name = "Erastin and RSL3 combined set")
normCPM <- t(V_combined$E)

#PCA plot
thisFig <- as_tibble(normCPM, rownames = "sample") %>% 
  mutate(cell = str_extract(sample, "^[:alnum:]*")) %>% 
  mutate(sample = str_replace(sample, "^[:alnum:]*_", "")) %>% 
  dplyr::select(sample, cell, everything())

PCA <- prcomp(thisFig[,c(-1,-2)])

fviz_pca_ind(PCA, axes = c(1, 2), geom = "point", axes.linetype = NA) +
  geom_point(aes(color = thisFig$cell), size = 3) +
  geom_text_repel(aes(label = str_replace(thisFig$sample, "_[:graph:]*$", ""))) +
  geom_mark_rect(aes(fill = thisFig$cell, label = thisFig$cell), size = 0.1, label.fontsize = 18, label.fontface = "plain") +
  labs(title = "", color = "Cell line") +
  theme_classic(18) +
  theme(legend.position = "none",
        legend.background = element_rect(color = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 18, color = "black"))

ggsave("PCA - Erastin and RSL3 in 5 cell lines.tiff", dpi = 200, height = 6, width = 6)

#quantify DEGs
decide_pval <- left_join(decide1, decide2, by = "genes")
Upreg <- decide_pval %>% summarise(across(-genes, ~sum(.x == 1, na.rm = TRUE))) %>% mutate(group = "upreg")
Downreg <- decide_pval %>% summarise(across(-genes, ~sum(.x == -1, na.rm = TRUE))) %>% mutate(group = "downreg")
Both <- decide_pval %>% summarise(across(-genes, ~sum(.x != 0, na.rm = TRUE))) %>% mutate(group = "both")
thisFig <- bind_rows(Upreg, Downreg) %>% 
  gather(-group, key = "treat", value = "n") %>% 
  mutate(treat = str_remove(treat, "_pval")) %>% 
  arrange(desc(n)) %>% 
  mutate(treat = fct_inorder(treat))

ggplot() +
  geom_col(data = thisFig %>% filter(group == "upreg"), aes(x = treat, y = n), fill = "salmon") +
  geom_col(data = thisFig %>% filter(group == "downreg"), aes(x = treat, y = -n), fill = "lightblue") +
  geom_text(data = thisFig %>% filter(group == "upreg"), aes(x = treat, y = n, label = n), color = "red", nudge_y = 75, size = 5) +
  geom_text(data = thisFig %>% filter(group == "downreg"), aes(x = treat, y = -n, label = n), color = "blue", nudge_y = -75, size = 5) +
  scale_y_continuous(limits = c(-1000, 1000)) +
  annotate(geom = "text", x = 5, y = 1000, label = str_c("Total # genes: ", dim(decide_pval)[1]), size = 6) +
  labs(x = "", y = "Downreg DEGs          Upreg DEGs") +
  theme_classic(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("PCA - Quantify DEGs.tiff", dpi = 200, height = 6, width = 6)

#UpSet plot
#the following are ferroptosis biomarkers validated in our recent publication (https://doi.org/10.1002/advs.202307263)
validated_biomarkers <-  c("DDIT4", "JDP2", "SLC38A2", "HMOX1", "VLDLR", "SLC6A9", "CHAC1", "MT2A", "OSGIN1", "PCK2",
                           "RPS6KA2", "SLC1A5", "ARHGEF2", "CALM2", "GARS", "PHGDH", "SPIN4", "TMEM154", "SACS", "NDC80",
                           "ERMAP", "ASNS", "PRR11", "MYC", "SLC7A11", "TXNRD1")

upset_plot <- DEG_combined %>% 
  as_tibble() %>% 
  na.omit() %>% 
  select(genes, contains("pval")) %>% 
  rename_with(.cols = -genes, .fn = ~str_remove(.x, "_pval")) %>% 
  mutate(across(.cols = -genes, .fns = ~ifelse(.x == 1, TRUE, FALSE))) %>% 
  select(genes, contains("MDA468"), contains("HCC70"), contains("BT549"), contains("HS578"), contains("SUM159")) %>% 
  mutate(genes = ifelse(genes %in% validated_biomarkers, genes, NA))

ComplexUpset::upset(upset_plot, intersect = colnames(upset_plot)[2:11], min_degree = 2, 
                    sort_intersections_by = c("degree", "cardinality"),
                    set_sizes = FALSE, sort_sets = FALSE, height_ratio = 0.8, name = "Groups",
                    annotations = list(
                      'text_panel' = ggplot(mapping = aes(y = 0, label = genes)) +
                        geom_text_repel(angle = 90, color = "red", hjust = 0, direction = "y", 
                                        min.segment.length = 10, size = 3, box.padding = 0.1) +
                        scale_y_continuous(limits = c(0, Inf)) +
                        theme_void()
                    ),
                    base_annotations = list(
                      'Intersection size'= (intersection_size() 
                                           + theme(panel.grid = element_blank(),
                                                   axis.title = element_text(color = "black", size = 16),
                                                   axis.text = element_text(color = "black", size = 16))
                                           + ylab('# Genes in intersection'))),
                    matrix = intersection_matrix() +
                      theme_bw(base_size = 18),
                    stripes = c("gray", "white"),
                    queries = list(upset_query(set = "SUM159_E", color = "red"),
                                   upset_query(set = "SUM159_R", color = "red"),
                                   upset_query(set = "HS578_E", color = "red"),
                                   upset_query(set = "HS578_R", color = "red"),
                                   upset_query(set = "BT549_E", color = "red"),
                                   upset_query(set = "BT549_R", color = "red"),
                                   upset_query(set = "MDA468_E", color = "blue"),
                                   upset_query(set = "MDA468_R", color = "blue"),
                                   upset_query(set = "HCC70_E", color = "blue"),
                                   upset_query(set = "HCC70_R", color = "blue")))

ggsave("UPSET of DEGs in all FINs.tiff", dpi = 300, height = 7, width = 14)

#-------------------------------------------------------------------------------------
# FIN heatmap
#-------------------------------------------------------------------------------------

#load data
DEG_combined_full <- DEG_combined %>% 
  na.omit() %>% 
  rename_with(.cols = -c(genes, contains("pval")), .fn = ~str_c(.x, "_exp")) %>% 
  left_join(DEG_combined_v3 %>% select(genes, contains("_t")), by = "genes") %>% 
  pivot_longer(-genes, names_to = c("treat", ".value"), names_pattern = "(.*)_(.*)$")

#heatmap ploting function - this function will plot the heatmaps shown in figure 3, 
#and will return the genes shown in the heatmap, used in further analysis
plot_heatmap <- function(Erastin = TRUE, RSL3 = TRUE, n = 75, width = 10, height = 4) {
  title <- paste0("Erastin-", Erastin, ", RSL3-", RSL3)
  thisFig_1 <- DEG_combined_full
  if(Erastin == FALSE) thisFig_1 <- thisFig_1 %>% filter(!str_detect(treat, "_E"))
  if(RSL3 == FALSE) thisFig_1 <- thisFig_1 %>% filter(!str_detect(treat, "_R"))
  trean_num <- thisFig_1 %>% distinct(treat) %>% dim()
  
  thisFig_2 <- thisFig_1 %>% 
    group_by(genes) %>% 
    mutate(n_up = sum(pval),
           median = mean(exp, na.rm = TRUE),
           median_t = mean(t, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(in_FerrDB = case_when(genes %in% FerrDB_marker ~ "Marker", 
                                 genes %in% FerrDB_driver ~ "Driver", 
                                 genes %in% FerrDB_suppresor ~ "Suppresor"))
  
  thisFig_up <- thisFig_2 %>% 
    filter(n_up > 2) %>% 
    slice_max(median, n = n * trean_num[1]) %>% 
    mutate(value = ifelse(abs(exp) > 2, 2 * (exp / abs(exp)), exp)) %>% 
    arrange(desc(n_up), desc(median)) %>% 
    mutate(genes = fct_inorder(genes)) %>% 
    mutate(color = ifelse(genes %in% validated_biomarkers, "red", "black")) %>% 
    mutate(x_label = glue("<span style='color:{color}'>{genes}</span>")) %>% 
    mutate(x_label = fct_inorder(x_label)) %>% 
    mutate(treat = fct_relevel(treat, "MDA468_JB2", "HCC70_JB2"))
  
  thisFig_down <- thisFig_2 %>% 
    filter(n_up < -1) %>% 
    slice_min(median, n = n * trean_num[1]) %>% 
    mutate(value = ifelse(abs(exp) > 2, 2 * (exp / abs(exp)), exp)) %>% 
    arrange(n_up, median) %>% 
    mutate(genes = fct_inorder(genes)) 
  
  ggplot(thisFig_up, aes(x = x_label, y = treat, fill = value)) +
    geom_tile(color = "black") +
    geom_text(aes(x = x_label, label = ifelse(pval == 1, "*", NA))) +
    geom_point(data = thisFig_up %>% filter(!is.na(in_FerrDB)) %>% distinct(genes, .keep_all = TRUE), 
               aes(x = x_label, color = in_FerrDB), y = trean_num[1] + 1) +
    scale_fill_gradient2(low = "lightblue", high = "green4", mid = "white", midpoint = 0, na.value = "grey50") +
    coord_cartesian(clip = "off") +
    labs(x = "", y = "", fill = "Log2 FC", color = "FerrDB") +
    theme_minimal(base_size = 14) +
    theme(legend.margin = unit(0, "cm"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(paste0("heatmap of upreg genes (", title, ").tiff"), dpi = 300, height = height, width = width)
  
  ggplot(thisFig_down, aes(x = genes, y = treat, fill = value)) +
    geom_tile(color = "black") +
    geom_text(aes(label = ifelse(pval == -1, "*", NA))) +
    geom_point(data = thisFig_down %>% filter(!is.na(in_FerrDB)) %>% distinct(genes, .keep_all = TRUE), 
               aes(x = genes, color = in_FerrDB), y = 10.8) +
    scale_fill_gradient2(low = "lightblue", high = "green4", mid = "white", midpoint = 0, na.value = "grey50") +
    coord_cartesian(clip = "off") +
    labs(x = "", y = "", fill = "Log2 FC", color = "FerrDB") +
    theme_minimal(base_size = 14) +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(paste0("heatmap of downreg genes (", title, ").tiff"), dpi = 300, height = height, width = width)
  
  genesets <- list("up_mean" = thisFig_up %>% distinct(genes) %>% pull(genes) %>% as.character(),
                   "down_mean" = thisFig_down %>% distinct(genes) %>% pull(genes) %>% as.character())
  return(genesets)
}

consensus_genes <- plot_heatmap(Erastin = TRUE, RSL3 = TRUE, n = 75, width = 14, height = 3)
