#this script will generate the transcriptomic similarity scores, used as predictors for the pipeline (figure 4)

library(tidyverse)
library(cmapR)
library(cogena)
library(ggpubr)
library(fgsea)
library(openxlsx)
library(archive)

#run the RNAseq analysis script first
source("3 - RNAseq analysis.R")

#---------------------------------------------------------------------------------------------------------------------------------
# CMAP analysis - prepare genesets from the CMAP data
# NOTE - this section takes ~2.5 hours to run. To skip this part, load the following file containing the results of this section:
CMAP_analysis <- read.csv(archive_read(archive = "Data/CMAP analysis.7z", file = "CMAP analysis.csv"))
# and skip to the next section
#---------------------------------------------------------------------------------------------------------------------------------

#prepare the gene sets from GSE92742. Parse parts of this huge file, keep only the "trt_sh.cgs"
#the source files needed for this section are not available in the Data folder due to their sizes (~ 20 GB). 
#They can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
ds_path <- "Data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx" 
sig_meta <- read.delim("Data/GSE92742_Broad_LINCS_sig_info.txt")
gene_data <- read.csv("Data/GSE92742_Broad_LINCS_gene_info.csv")
relevant_IDs <- sig_meta %>% 
  filter(pert_type == "trt_sh.cgs") %>% 
  select(sig_id) %>% 
  as_vector()
col_meta <- read_gctx_meta(ds_path, dim="col")
idx <- which(col_meta$id %in% relevant_IDs)
my_ds <- parse_gctx(ds_path, cid = idx)

x1 <- my_ds@mat
x1 <- as_tibble(x1, rownames = "ID")
write_csv(x1, "GSE92742 parced.csv")

TS_geneset_list_up_20 <- list()
TS_geneset_list_up_50 <- list()
TS_geneset_list_up_100 <- list()
TS_geneset_list_up_200 <- list()
TS_geneset_list_down_20 <- list()
TS_geneset_list_down_50 <- list()
TS_geneset_list_down_100 <- list()
TS_geneset_list_down_200 <- list()

for(n in 1:dim(x1)[2]) {
  currCol <- colnames(x1)[n+1]
  
  X1DF <- x1 %>% 
    select(ID, any_of(currCol)) %>% 
    arrange(desc(across(any_of(currCol)))) %>% 
    mutate(ID = as.numeric(ID)) %>% 
    left_join(gene_data, by = "ID") %>% 
    dplyr::select(name)
  
  TS_geneset_list_up_20[n] <- slice_head(X1DF, n = 20)
  TS_geneset_list_up_50[n] <- slice_head(X1DF, n = 50)
  TS_geneset_list_up_100[n] <- slice_head(X1DF, n = 100)
  TS_geneset_list_up_200[n] <- slice_head(X1DF, n = 200)
  TS_geneset_list_down_20[n] <- slice_tail(X1DF, n = 20)
  TS_geneset_list_down_50[n] <- slice_tail(X1DF, n = 50)
  TS_geneset_list_down_100[n] <- slice_tail(X1DF, n = 100)
  TS_geneset_list_down_200[n] <- slice_tail(X1DF, n = 200)
  
  names(TS_geneset_list_up_20)[n] <- currCol
  names(TS_geneset_list_up_50)[n] <- currCol
  names(TS_geneset_list_up_100)[n] <- currCol
  names(TS_geneset_list_up_200)[n] <- currCol
  names(TS_geneset_list_down_20)[n] <- currCol
  names(TS_geneset_list_down_50)[n] <- currCol
  names(TS_geneset_list_down_100)[n] <- currCol
  names(TS_geneset_list_down_200)[n] <- currCol

  print(paste0("creating: ", n, " / ", dim(x1)[2], " (", round(n/dim(x1)[2]*100,3), "%)"))
}

gmtlist2file(TS_geneset_list_up_20, "CMAP TS complete gene set up (20 genes).gmt")
gmtlist2file(TS_geneset_list_up_50, "CMAP TS complete gene set up (50 genes).gmt")
gmtlist2file(TS_geneset_list_up_100, "CMAP TS complete gene set up (100 genes).gmt")
gmtlist2file(TS_geneset_list_up_200, "CMAP TS complete gene set up (200 genes).gmt")
gmtlist2file(TS_geneset_list_down_20, "CMAP TS complete gene set down (20 genes).gmt")
gmtlist2file(TS_geneset_list_down_50, "CMAP TS complete gene set down (50 genes).gmt")
gmtlist2file(TS_geneset_list_down_100, "CMAP TS complete gene set down (100 genes).gmt")
gmtlist2file(TS_geneset_list_down_200, "CMAP TS complete gene set down (200 genes).gmt")

#load the genesets (created in the previous section, or loaded here from the saved gmt files)
GMT_files <- c("CMAP TS complete gene set up (20 genes).gmt", "CMAP TS complete gene set up (50 genes).gmt", 
              "CMAP TS complete gene set up (100 genes).gmt", "CMAP TS complete gene set up (200 genes).gmt", 
              "CMAP TS complete gene set down (20 genes).gmt", "CMAP TS complete gene set down (50 genes).gmt", 
              "CMAP TS complete gene set down (100 genes).gmt", "CMAP TS complete gene set down (200 genes).gmt")

#Function for calculating enrichment using CAMERA
CameraFunction <- function(contr.matrix1, V1, design1, idx1,
                           contr.matrix2, V2, design2, idx2) {
  cam.list <- list()
  for(n in 1:dim(contr.matrix1)[2]) {
    name <- names(contr.matrix1[,n])[contr.matrix1[,n] == 1]
    
    cam.list[[n]] <- camera(V1, idx1, design1, contrast=contr.matrix1[,n]) %>% 
      as_tibble(rownames = "pathway")   %>% 
      mutate(result = ifelse(Direction == "Up", PValue, -PValue)) %>% 
      dplyr::select(pathway, result) 
    
    names(cam.list[[n]])[names(cam.list[[n]]) == "result"] <- name
    print(n)
    
  }
  cam.object1 <- Reduce(full_join, cam.list)
  
  cam.list <- list()
  for(n in 1:dim(contr.matrix2)[2]) {
    name <- names(contr.matrix2[,n])[contr.matrix2[,n] == 1]
    
    cam.list[[n]] <- camera(V2, idx2, design2, contrast=contr.matrix2[,n]) %>% 
      as_tibble(rownames = "pathway") %>% 
      mutate(result = ifelse(Direction == "Up", PValue, -PValue)) %>% 
      dplyr::select(pathway, result) 
    
    names(cam.list[[n]])[names(cam.list[[n]]) == "result"] <- name
    print(n+8)
  }
  cam.object2 <- Reduce(full_join, cam.list)
  
  cam.object <- full_join(cam.object1, cam.object2, by = "pathway")
  cam.object.log <- cam.object %>% 
    mutate(across(where(is.numeric), ~ifelse(.x>0, -log10(.x), log10(abs(.x)))))
  
  return(cam.object.log)
}

#Calculate enrichement scores
Camera_result <- tibble("pathway" = names(gmt2list(paste0("Data/", GMT_files[1])))) 

for (n in 1:8) {
  title <-  str_extract(GMT_files[n], "gene set [:alpha:]* \\([:digit:]* genes\\)")
  title_short <- paste(str_extract(title, "up|down"), str_extract(title, "[:digit:]*(?= genes)"), sep = "_")
  CMAP_geneset <- gmt2list(paste0("Data/", GMT_files[n]))
  idx1 <- ids2indices(CMAP_geneset, identifiers = V1$genes$genes)
  idx2 <- ids2indices(CMAP_geneset, identifiers = V2$genes$genes)
  x <- CameraFunction(contr.matrix1, V1, design1, idx1, 
                                        contr.matrix2, V2, design2, idx2, excel = TRUE, title = title)
  x1 <- x %>% 
    rename_with(~paste0(.x, "_", title_short), .cols = -pathway)
  
  Camera_result <- Camera_result %>% 
    left_join(x1)
  
  print(paste0("--------------------finished: n = ", n))
}

CMAP_analysis <- Camera_result %>% 
  mutate(ID = str_replace(pathway, "^CGS[:digit:]*_", ""),
         ID = str_replace(ID, "_[:alnum:]*:", "_"),
         ID = str_replace(ID, ":[:graph:]*$", "")) %>% 
  filter(str_detect(pathway, "96H")) %>% 
  dplyr::select(ID, everything(), -pathway)

#---------------------------------------------------------------------------------------------------------------------------------
# Transcriptomic similarity scores analysis
#---------------------------------------------------------------------------------------------------------------------------------

#load the results of the previous section (if skipped)
CMAP_analysis <- read_csv("Data/CMAP analysis.csv")

#long version of results tibble
CMAP_analysis_long <- CMAP_analysis %>% 
  mutate(cell = sapply(str_split(ID, "_"), function(x) x[1])) %>% 
  mutate(gene = sapply(str_split(ID, "_"), function(x) x[2])) %>% 
  select(-ID) %>% 
  pivot_longer(-c(cell, gene), names_to = c("pert", ".value", "size"), names_pattern = "(.*_.)_(.*)_(.*)")
  
#extracting the maximum value for each gene
CMAP_analysis_max <- CMAP_analysis %>%  
  mutate(cell = sapply(str_split(ID, "_"), function(x) x[1])) %>% 
  mutate(gene = sapply(str_split(ID, "_"), function(x) x[2])) %>% 
  dplyr::select(ID, cell, gene, everything())  %>% 
  group_by(gene) %>% 
  summarise(across(contains("_"), max)) 

CMAP_analysis_max_ranked <- CMAP_analysis_max %>% 
  mutate(across(contains("_"), ~min_rank(.x)))

#combine the up and down gene sets. 
#The metric diff is UP - DOWN, so +4000 is the best, -4000 is the worse
CMAP_analysis_max_combined <- CMAP_analysis_max_ranked %>% 
  gather(-gene, key = "pert", value = "rank") %>% 
  mutate(up_down = str_extract(pert, "up|down"), .after = pert) %>% 
  mutate(pert = str_remove(pert, "_up|_down")) %>% 
  spread(key = "up_down", value = "rank") %>% 
  mutate(total = up - down) %>% 
  dplyr::select(-up, -down) %>% 
  spread(key = "pert", value = "total")

#write to file and proceed in the pipeline buildup script
write_csv(CMAP_analysis_max_combined, "CMAP analysis max combined.csv")

#---------------------------------------------------------------------------------------------------------------------------------
# plots shown in figure 4 and S3
#---------------------------------------------------------------------------------------------------------------------------------

cells <- CMAP_analysis_long %>% 
  group_by(cell) %>% 
  dplyr::count() %>% 
  filter(n > 100000) %>% 
  pull(cell)

#CMAP specific genes
thisFig <- CMAP_analysis_long %>% 
  filter(size == 20, cell %in% cells) %>% 
  filter(gene %in% c("GCH1", "GPX4")) %>% 
  group_by(gene, pert) %>% 
  mutate(is_max = ifelse(up == max(up), "max", "none")) %>% 
  ungroup()

ggplot(thisFig, aes(x = cell, y = up)) +
  geom_point(data = thisFig %>% filter(is_max != "max"), color = "black") +
  geom_point(data = thisFig %>% filter(is_max == "max"), color = "red") +
  facet_grid(gene ~ pert, scales = "free") +
  labs(x = "Cell lines in CMAP", y = "Transcriptomic similarity scores") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("CMAP scores per cell line - examples.tiff", dpi = 300, width = 14, height = 4)

#demonstrate the ranking
thisFig <- left_join(CMAP_analysis_max %>% select(gene, BT549_E_up_20, BT549_E_down_20),
                     CMAP_analysis_max_ranked %>% select(gene, BT549_E_up_20, BT549_E_down_20), by = "gene")

p1 <- ggplot(thisFig, aes(x = BT549_E_up_20.x, y = BT549_E_up_20.y)) +
  geom_point() +
  geom_point(data = thisFig %>% filter(gene %in% c("GPX4", "GCH1")), color = "red", size = 3) +
  geom_text_repel(data = thisFig %>% filter(gene %in% c("GPX4", "GCH1")), 
                  aes(label = gene), color = "red", size = 5, nudge_x = 2, nudge_y = 200) +
  labs(x = "Maximum similarity score", y = "Ranking", title = "Upregulated gene set scores") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))

p2 <- ggplot(thisFig, aes(x = BT549_E_down_20.x, y = BT549_E_down_20.y)) +
  geom_point() +
  geom_point(data = thisFig %>% filter(gene %in% c("GPX4", "GCH1")), color = "red", size = 3) +
  geom_text_repel(data = thisFig %>% filter(gene %in% c("GPX4", "GCH1")), 
                  aes(label = gene), color = "red", size = 5, nudge_x = 2, nudge_y = -200) +
  labs(x = "Maximum similarity score", y = "", title = "Downregulated gene set scores") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))

ggarrange(p1, p2 + rremove("y.text"), align = "h")
ggsave("CMAP scores ranking - examples.tiff", dpi = 300, width = 8, height = 4)

#plot the up/down scores
make_gradient <- function(deg = 45, n = 100, cols = blues9) {
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(data = rep(seq(0, 1, length.out = n) * cos(rad), n), byrow = TRUE, ncol = n) +
    matrix(data = rep(seq(0, 1, length.out = n) * sin(rad), n), byrow = FALSE, ncol = n)
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(image = mat,
                   width = unit(1, "npc"),
                   height = unit(1, "npc"), 
                   interpolate = TRUE
  )
}
cols <- colorRampPalette(c("skyblue3", "skyblue", "#BCDAA7", "#70AD47"))(100)

g1 <- make_gradient(deg = 45, n = 500, cols = cols)
g2 <- make_gradient(deg = 0, n = 500, cols = cols)

thisFig <- CMAP_analysis_max_ranked %>% 
  pivot_longer(-gene, names_to = c("pert", ".value", "size"), names_pattern = "(.*_.)_(.*)_(.*)") %>% 
  filter(size == 20) %>% 
  mutate(label = case_when(gene %in% c("GPX4", "GCH1") ~ "anti",
                           gene %in% c("ACSL4") ~ "pro",
                           TRUE ~ "none")) %>% 
  mutate(total = up - down) %>% 
  filter(pert == "HCC70_R") %>% 
  rowwise() %>% 
  mutate(y.value = ifelse(label == "none", runif(1), 0.5)) %>% 
  ungroup()

p1 <- ggplot(thisFig, aes(x = up, y = down)) +
  annotation_custom(grob = g1, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  geom_point(data = thisFig %>% filter(label == "none"), color = "gray60", alpha = 0.5, shape = 1) +
  geom_point(data = thisFig %>% filter(label == "anti"), color = "darkgreen", size = 3) +
  geom_point(data = thisFig %>% filter(label == "pro"), color = "black", size = 3) +
  geom_text_repel(data = thisFig %>% filter(label == "anti"), aes(label = gene), color = "darkgreen", size = 5, fontface = "bold") +
  geom_text_repel(data = thisFig %>% filter(label == "pro"), aes(label = gene), color = "black", size = 5, fontface = "bold") +
  geom_label(data = thisFig %>% distinct(pert), aes(label = pert), x = 50, y = 200, size = 5, hjust = 0) +
  labs(x = "Upregulated signatures ranks", y = "Downregulated signatures ranks") +
  theme_bw(base_size = 18) +
  theme(axis.text = element_text(color = "black", size = 14),
        strip.text = element_blank())

#total scores plot
p2 <- ggplot(thisFig, aes(x = total, y = y.value)) +
  annotation_custom(grob = g2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  geom_point(data = thisFig %>% filter(label == "none"), color = "gray60", alpha = 0.5, shape = 1) +
  geom_point(data = thisFig %>% filter(label == "anti"), color = "darkgreen", size = 3) +
  geom_point(data = thisFig %>% filter(label == "pro"), color = "black", size = 3) +
  geom_text_repel(data = thisFig %>% filter(label %in% c("anti", "pro")), aes(label = gene, color = label), 
                  size = 5, fontface = "bold") +
  scale_color_manual(values = c("darkgreen", "black")) +
  labs(x = "Transcriptomic similarity score", y = "") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(color = "black", size = 14))

ggarrange(p1, p2, ncol = 1, align = "hv", heights = c(0.8, 0.2))
ggsave("CMAP scores rankings per perturbation (both plots, just HCC_70).tiff", dpi = 300, width = 6, height = 6.5)

#plot FerrDB distribution in each pertutbation
FerrDB_driver <- read_csv("Data/FerrDB_info_driver.csv") %>% pull(Symbol)
FerrDB_suppresor <- read_csv("Data/FerrDB_info_suppressor.csv") %>% pull(Symbol)
FerrDB_marker <- read_csv("Data/FerrDB_info_marker.csv") %>% pull(Symbol)

FerrDB_list <- list("FerrDB_drivers" = FerrDB_driver, "FerrDB_suppresors" = FerrDB_suppresor, "FerrDB_markers" = FerrDB_marker)
thisFig <- CMAP_analysis_max_ranked %>% 
  pivot_longer(-gene, names_to = c("pert", ".value", "size"), names_pattern = "(.*_.)_(.*)_(.*)") %>% 
  filter(size == 20) %>% 
  mutate(total = up - down)
pert <- thisFig %>% distinct(pert) %>% pull(pert) %>% as.character()

FGSEA_list <- list()
for(i in seq_along(pert)) {
  tempDF <- thisFig %>% 
    filter(pert == pert[i])
  ranks <- setNames(tempDF$total, tempDF$gene)
  FGSEA <- fgsea(pathways = FerrDB_list, stats = ranks, maxSize = 1000, nproc = 0, scoreType = "std")
  FGSEA_list[[i]] <- FGSEA %>% mutate(pert = pert[i])
  
}

FGSEA_tib <- purrr::reduce(FGSEA_list, bind_rows) %>% 
  as_tibble() %>% 
  select(pathway, NES, pert)
FGSEA_tib$pert <- factor(FGSEA_tib$pert, levels = pert_sort$pert, ordered = TRUE)

ggplot(FGSEA_tib %>% filter(pathway == "FerrDB_suppresors"), aes(x = NES, y = 0.5)) +
  annotation_custom(grob = g2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  geom_point(color = "red", size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Normalized enrichment score", y = "", color = "FerrDB genes") +
  facet_grid(pert ~ .) +
  theme_minimal(base_size = 16) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_blank(),
        axis.text = element_text(color = "black"))

ggsave("CMAP scores rankings per perturbation (FerrDB suppressors, FGSEA version).tiff", dpi = 300, width = 4, height = 5)

#save CMAP scores in supplemental tables (Supplemental table S4)
thisFig <- CMAP_analysis %>% 
  mutate(cell = sapply(str_split(ID, "_"), function(x) x[1])) %>% 
  mutate(gene = sapply(str_split(ID, "_"), function(x) x[2])) %>% 
  select(gene, cell, matches("up_20$"), matches("down_20$")) %>% 
  rename_with(.cols = everything(), .fn = ~str_remove(.x, "_20$"))


wb <- createWorkbook()
options("openxlsx.numFmt" = "0.00")
addWorksheet(wb, "CMAP", gridLines = FALSE)
writeData(wb, sheet = "CMAP", x = "Connectivity map scores", startCol = 1, startRow = 1)
writeDataTable(wb, sheet = "CMAP", x = thisFig, startCol = 1, startRow = 3, 
               tableStyle = "TableStyleLight1", bandedRows = FALSE)

saveWorkbook(wb, file = "CMAP scores.xlsx", overwrite = TRUE)

