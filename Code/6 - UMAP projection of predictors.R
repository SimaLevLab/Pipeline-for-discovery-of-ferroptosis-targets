#this script is the main compatational pipeline, using the predictors to generate the UMAP shown in figure 5

library(tidyverse)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(ggridges)
library(rstatix)
library(umap)
library(tidytext)
library(cmapR)
library(cogena)
library(rentrez)
library(gplots)

#-----------------------------------------------------------------------------------------
# Gather info from other scripts (run the script, or load the data they generate)
#-----------------------------------------------------------------------------------------

#1. run "RNAseq analysis - Main script"
source("3 - RNAseq analysis.R")
#2. Transcriptomic similarity scores - run the script "5 - Transcriptomic similarity scores.R", or load the script result here:
CMAP <- read_csv("Data/CMAP analysis max combined.csv")
#3, Correlation scores - run the script "2 - Correlations of AUCs to gene data.R", or load the script result here:
Correlations <- read_csv("Data/Correlation of gene expression and depScore to AUCs (2023).csv")
#4, Geneshot query (done in the geneshot site, results saved in the file below)
Geneshot_query <- read_csv("Data/Geneshot query for ferroptosis terms.csv")
Geneshot_terms <- c("Ferroptosis", "Lipid peroxidation", "Iron", "Glutathione", "GPX4", "SLC7A11")
#other data
Achilles_common_essential <- read_csv("Data/Achilles_common_essentials.csv") %>% 
  mutate(gene = str_extract(gene, "^[:alnum:]*")) %>% pull(gene)

#make gradient function
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

#-----------------------------------------------------------------------------------------
# Predictor selection step (figure 4D)
#-----------------------------------------------------------------------------------------

#create pairwise pearson correlations between perturbations
cor_matrix_FC <- DEG_combined %>% 
  select(-genes, -contains("pval")) %>% 
  na.omit() %>% 
  cor(method = "pearson")

weights <- cor_matrix_FC %>% 
  as_tibble(rownames = "pert") %>% 
  mutate(across(-pert, ~ifelse(.x == 1, 0, .x))) %>% 
  mutate(across(-pert, ~ifelse(.x < 0, 0, .x))) %>% 
  mutate(sum = rowSums(across(-pert))) %>% 
  mutate(weights = sum / sum(sum) * 100) %>% 
  select(pert, weights) %>% 
  arrange(weights) %>% 
  mutate(pert = fct_reorder(pert, weights))

p1 <- ggplot(weights, aes(x = weights, y = pert)) +
  geom_col(aes(fill = weights < 9)) +
  scale_fill_manual(values = c("green4", "gray")) +
  labs(x = "Degree of Similarity (%)", y = "") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 14))

matrix <- cor_matrix_FC %>% 
  as_tibble(rownames = "pert") %>% 
  gather(-pert, key = "treat", value = "R")
sorted_pert <- levels(weights$pert)
matrix$pert <- factor(matrix$pert, levels = sorted_pert)
matrix$treat <- factor(matrix$treat, levels = sorted_pert)

p2 <- ggplot(matrix, aes(x = pert, y = treat, fill = R)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(R, 2)), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "green4") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 16) +
  theme(legend.position = "left",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggarrange(p2, p1, align = "hv")
ggsave("Perturbations consensus matrix and similarity scores.tiff", dpi = 300, height = 4, width = 10)

#-----------------------------------------------------------------------------------------
# Text mining validation for predictor selection step (figure 4E)
#-----------------------------------------------------------------------------------------

#make a histogram for the distributions of relevant terms upon random selection of 100 genes
set.seed(42)
iterations <- 2000
sample_size <- 100
citation_cutoff <- 0 #genes with less citations than this will be ignored
geneshot_density <- tibble()
geneshot_query_cmap <- Geneshot_query %>% 
  filter(gene %in% CMAP$gene) %>% #pick only genes appearing in CMAP
  mutate(across(-gene, ~ifelse(.x > citation_cutoff, .x, NA)))

for(i in 1:iterations) {
  temp <- geneshot_query_cmap %>% 
    slice_sample(n = sample_size, replace = FALSE) %>% 
    summarise(across(-gene, ~sample_size - sum(is.na(.x)))) %>% 
    mutate("iteration" = i)
  geneshot_density <- bind_rows(geneshot_density, temp)
  cat("\r",i, " / ", iterations)
}

Geneshot_random_results <- geneshot_density %>% 
  mutate(Average = rowMeans(across(-iteration))) %>% 
  gather(-iteration, key = "term", value = "value") %>% 
  mutate(value = (value / sample_size) * 100) %>% 
  mutate(category = "Random") %>% 
  dplyr::select(-iteration)

Geneshot_random_means <- Geneshot_random_results %>% 
  group_by(category, term) %>% 
  get_summary_stats()
  
#Perform all permutations of CMAP scors (20 geneset sizes) to determine the set of best predictors
#helper function - this function will take a set of the transcriptomic similarity predictors (out of the 10), 
#add the 6 correlations predictors, and calculate the ratio of ferroptosis related terms among the
#50 genes with the best predictors value.
Quantify_terms_ratio_method_A <- function(selected_pert) {
  citation_cutoff <- 0
  sample_size <- 50
  selected_pert <- str_split(selected_pert, ",", simplify = TRUE)[1,]
  print(paste("Measuring for combination of", paste(selected_pert, collapse = " ")))
  
  tempDF <- Correlations %>%
    full_join(CMAP, by = "gene") %>%
    select(gene, contains("AUC"), all_of(selected_pert)) %>%
    na.omit() %>%
    left_join(Geneshot_query, by = "gene") %>%
    mutate(across(any_of(Geneshot_terms), ~ifelse(.x > citation_cutoff, cur_column(), NA))) %>%
    mutate(mean_r_expression = rowMeans(across(contains("AUC to expression"))),
           mean_r_depscore = rowMeans(across(contains("AUC to Depscore"))),
           biggest_r = pmax(mean_r_expression, mean_r_depscore),
           mean_CMAP_scores = rowMeans(across(all_of(selected_pert)))) %>%
    filter(biggest_r > 0.5) %>%
    slice_max(mean_CMAP_scores, n = sample_size) %>%
    select(gene, any_of(Geneshot_terms)) %>%
    gather(-gene, key = "term", value = "value") %>%
    group_by(value) %>%
    dplyr::count() %>%
    na.omit() %>%
    mutate(fraction = (n / sample_size) * 100) %>%
    select(-n) %>%
    spread(key = "value", value = "fraction")
  
  return(tempDF)
}  

#run the helper function on all combinations of the 10 transcriptomic similarity predictors
selected_pert <- tibble("pert" = colnames(CMAP)) %>% 
  filter(str_detect(pert, "_20$")) %>% 
  pull(pert)

combinations_result <- list()
for(t in 1:10) {
  pert_combinations <- combn(selected_pert, m = t) 
  Tb_result <- pert_combinations %>% 
    t() %>% 
    as_tibble() %>% 
    unite(col = "combination", everything(), sep = ",") %>% 
    mutate(result <- purrr::map(combination, Quantify_terms_ratio_method_A))
  combinations_result[[t]] <- Tb_result %>% 
    unnest()
}
combinations_result_final <- Reduce(x = combinations_result, f = bind_rows)
pert_picker <- combinations_result_final %>% mutate(mean = rowMeans(across(any_of(Geneshot_terms)))) 

#Selected best predictors set
best_predictors_set <- c("BT549_E_20", "BT549_R_20", "HCC70_E_20", "HCC70_R_20", "MDA468_R_20", "SUM159_E_20", "SUM159_R_20")
best_predictors_set_collapsed <- paste0(best_predictors_set, collapse = ",")

all_predictors_set <- c("BT549_E_20", "BT549_R_20", "HCC70_E_20", "HCC70_R_20", "HS578_E_20", "HS578_R_20", 
                        "MDA468_E_20", "MDA468_R_20", "SUM159_E_20", "SUM159_R_20")
all_predictors_set_collapsed <- paste0(all_predictors_set, collapse = ",")

#generate histograms
Geneshot_pert_results <- combinations_result_final %>% 
  gather(-combination, key = "term", value = "value") %>% 
  mutate(category = "Perturbation") %>% 
  dplyr::select(-combination)

Geneshot_pert_means <- Geneshot_pert_results %>% 
  group_by(category, term) %>% 
  get_summary_stats()

total_results <- bind_rows(Geneshot_random_results %>% filter(term != "Average"), Geneshot_pert_results)
total_means <- bind_rows(Geneshot_random_means %>% filter(term != "Average"), Geneshot_pert_means)

selected_pert_summs <- combinations_result_final %>% 
  filter(combination == best_predictors_set_collapsed) %>% 
  select(-combination) %>% 
  gather(key = "term", value = "value") %>% 
  mutate(categoty = "Selected_pert")

all_pert_summs <- combinations_result_final %>% 
  filter(combination == all_predictors_set_collapsed) %>% 
  select(-combination) %>% 
  gather(key = "term", value = "value") %>% 
  mutate(categoty = "All_pert")

summary_table <- total_means %>% 
  mutate(value = str_c(round(mean, 1), "% \u00B1 ", round(sd, 1))) %>% 
  select(term, category, mean, sd, value) %>% 
  pivot_wider(names_from = "category", values_from = c("mean", "sd", "value")) %>% 
  left_join(selected_pert_summs %>% select(term, "selected_pert" = value), by = "term") %>% 
  mutate(z_score_rand = ((selected_pert - mean_Random) / sd_Random)) %>% 
  mutate(p_value_rand = 2 * pnorm(-abs(z_score_rand))) %>% 
  mutate(z_score_pert = ((selected_pert - mean_Perturbation) / sd_Perturbation)) %>% 
  mutate(p_value_pert = 2 * pnorm(-abs(z_score_pert))) %>% 
  mutate(selected_pert = str_c(round(selected_pert, 1), "%")) %>% 
  mutate(sig_level = case_when(p_value_rand < 0.01 ~ "**",
                               p_value_rand < 0.05 ~ "*",
                               TRUE ~ ""))

#plot ridgeplot - using the distribution of all perturbation and the random distribution
thisFig <- total_results %>% 
  mutate(term = fct_relevel(term, "GPX4", "Ferroptosis", "Lipid peroxidation", "SLC7A11", "Iron", "Glutathione")) %>% 
  mutate(category = ifelse(category == "Perturbation", "All predictor combinations", category))

ggplot(thisFig, aes(x = value, y = term)) + 
  geom_density_ridges2(aes(fill = category), calc_ecdf = TRUE, alpha = 0.75) +
  geom_point(data = selected_pert_summs, color = "red", size = 2) +
  geom_point(data = all_pert_summs, color = "black", size = 2) +
  geom_text(data = summary_table, aes(label = value_Random, color = p_value_rand < 0.05), x = 75, hjust = 0, size = 4) +
  geom_text(data = summary_table, aes(label = str_c(selected_pert, sig_level), color = p_value_rand < 0.05), 
            x = 95, hjust = 0, size = 4) +
  scale_fill_manual(values = c("pink", "skyblue")) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  labs(x = "% of genes positive to the term", y = "", fill = "") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top",
        plot.margin = margin(0,150,0,0), 
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, color = "black"))

ggsave("Density of perturbation combinations vs random and perturbations (ridge plot).tiff", 
       dpi = 300, height = 4, width = 8)  

#-----------------------------------------------------------------------------------------
#combining files for clustering
#-----------------------------------------------------------------------------------------
selected_pert <- c("HCC70_E_20", "MDA468_R_20", "BT549_E_20", "SUM159_E_20", "BT549_R_20", "SUM159_R_20", "HCC70_R_20")

#combining lists
GeneList <- Correlations %>% 
  full_join(CMAP, by = "gene") %>% 
  dplyr::select(gene, contains("AUC"), all_of(selected_pert)) %>%
  left_join(Geneshot_query, by = "gene")

#calculate ranking of genes
best_predictors_rank <- GeneList %>% 
  mutate(mean_r_expression = rowMeans(across(contains("AUC to expression"))),
         mean_r_depscore = rowMeans(across(contains("AUC to Depscore"))),
         biggest_r = pmax(mean_r_expression, mean_r_depscore, na.rm = TRUE),
         mean_CMAP_scores = rowMeans(across(all_of(selected_pert)))) %>% 
  mutate(across(any_of(Geneshot_terms), ~ifelse(.x > 0, cur_column(), NA))) %>% 
  mutate(priority_term = coalesce(Ferroptosis, `Lipid peroxidation`, Glutathione, Iron, GPX4, SLC7A11)) %>% 
  dplyr::select(gene, mean_r_expression, mean_r_depscore, biggest_r, mean_CMAP_scores, priority_term) %>% 
  mutate(priority_term = ifelse(is.na(priority_term), "none", priority_term),
         priority_term = fct_relevel(priority_term, "none", after = Inf)) %>% 
  filter(!is.na(mean_CMAP_scores)) %>% 
  mutate(rank_corr_expression = min_rank(mean_r_expression),
         rank_corr_depscoree = min_rank(mean_r_depscore),
         rank_corr = pmax(rank_corr_expression, rank_corr_depscoree, na.rm = TRUE),
         rank_CMAP_score = min_rank(mean_CMAP_scores), 
         mean_rank = (rank_corr + rank_CMAP_score) / 2,
         rank_total = min_rank(mean_rank))

#get the top ranking genes from above
rank_cutoff <- 150
top_ranking_genes <- best_predictors_rank %>% 
  slice_max(rank_total, n = rank_cutoff) %>% 
  pull(gene)

#clustering by UMAP
set.seed(42)

thisDF.prep <- GeneList %>% 
  arrange(gene) %>% 
  select(-any_of(Geneshot_terms)) %>% 
  na.omit() %>% 
  mutate(across(where(is.numeric), scale))

thisDF.data <- thisDF.prep %>% select(where(is.numeric))
thisDF.labels <- thisDF.prep %>% select(gene)
thisDF.map <- umap(thisDF.data, preserve.seed = TRUE, random_state = 5023, min_dist = 0.05)
km.clusters <- kmeans(thisDF.map$layout, 5, nstart = 25)
result <- cbind(thisDF.labels, thisDF.map$layout, km.clusters$cluster) %>% 
  mutate(zoom = NA) %>% 
  dplyr::rename("cluster" = "km.clusters$cluster") %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  left_join(best_predictors_rank %>% select(gene, rank_total))

#plot the UMAP
selected_labels <- c("GPX4", "SLC7A11", "HMOX1", "ACSL4", "CBS", "ZEB1", "GCH1")

ggplot(result, aes(x = `1`, y = `2`)) + 
  geom_point(aes(fill = rank_total), size = 2, shape = 21, alpha = 0.8, color = "gray") +
  geom_label_repel(data = result %>% filter(gene %in% selected_labels), aes(label = gene), 
                   size = 4, color = "black", fill = "gray90", nudge_y = 0.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "skyblue3") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank())

ggsave("UMAP of best predictors (colored by best predictor ranks).tiff", dpi = 300, width = 6, height = 6)

#-----------------------------------------------------------------------------------------
# plot UMAP with zoom-in on nodes
#-----------------------------------------------------------------------------------------

validated_genes <- c("DCXR", "DECR1", "SORD", "NDUFV1", "CANT1", "FRAT1", "TNFRSF18", "AMD1", "CTBP1", "PDAP1")

get_pubmed_info <- function(genes, Pubmed_terms) {
  pubmed_search <- list()
  for(i in seq_along(genes)) {
    print(genes[i])
    
    for(j in seq_along(Pubmed_terms)) {
      print(str_c("---", Pubmed_terms[j]))
      search_term <- paste0(genes[i], " AND ", Pubmed_terms[j])
      pubmed_result <- entrez_search(db = "pubmed", term = search_term)
      pubmed_search[[i*10 + j]] <- tibble("gene" = genes[i], "term" = Pubmed_terms[j], "count" = pubmed_result$count)
    }
    Sys.sleep(0.1)
  }
  pubmed_search_tib <- purrr::reduce(pubmed_search, bind_rows)
}

nodes_nearest_neighbors <- function(zoom_on_gene, distance_plot = FALSE, table = FALSE, excel = FALSE) {
  coord.x <- result[result$gene == zoom_on_gene, ]$`1`
  coord.y <- result[result$gene == zoom_on_gene, ]$`2`
  neighborhood <- 0.4
  
  result.text <- result %>% 
    filter(`1` > coord.x - neighborhood & `1` < coord.x + neighborhood) %>% 
    filter(`2` > coord.y - neighborhood & `2` < coord.y + neighborhood) %>% 
    mutate(zoom = TRUE)
  
  thisFig <- result %>% 
    as_tibble() %>% 
    left_join(best_predictors_rank, by = "gene")
  
  ggplot(thisFig, aes(x = `1`, y = `2`)) + 
    geom_point(aes(fill = rank_total.y), size = 2, shape = 21, alpha = 0.8, color = "gray") +
    geom_text_repel(data = result.text %>% filter(!(gene %in% selected_labels)), 
                    aes(label = gene, color = gene %in% validated_genes), size = 4, max.overlaps = 30) +
    geom_label_repel(data = result.text %>% filter(gene %in% selected_labels), 
                     aes(label = gene), fill = "white",  color = "red", size = 4, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "skyblue3") +
    scale_color_manual(values = c("gray40", "red")) +
    facet_zoom(xlim = c(coord.x - neighborhood, coord.x + neighborhood), 
               ylim = c(coord.y - neighborhood, coord.y + neighborhood), 
               zoom.size = 1, zoom.data = zoom) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "gray95"),
          panel.grid = element_blank())
  
  ggsave(paste0("UMAP of ", zoom_on_gene, " node predictors, with zoom-in on node.tiff"), 
        dpi = 200, width = 12, height = 6)
  
  #generate list of closest neighbors to the node
  node_neighbors <- result %>% 
    mutate(distance = sqrt((`1` - coord.x)^2 + (`2` - coord.y)^2)) %>% 
    select(gene, distance) %>% 
    arrange(distance) %>% 
    slice_min(distance, n = 75) %>% 
    left_join(GeneList, by = "gene") 
  
  #search pubmed for ferroptosis publications
  if(distance_plot) {
    set_entrez_key("") #write here the pubmed access key
    pubmed_filename <- str_c("Data/Pubmed search for ", zoom_on_gene, " node neighbors.RData")
    if(file.exists(pubmed_filename)) {
      load(pubmed_filename)
    } else {
      pubmed_data <- get_pubmed_info(as.vector(node_neighbors$gene), Geneshot_terms)
      save(pubmed_data, file = pubmed_filename)
    }
    
    #plot the pubmed plot
    thisFig <- node_neighbors %>%
      as_tibble() %>% 
      dplyr::select(gene, distance) %>% 
      slice_head(n = 50) 
    
    thisFig2 <- thisFig %>% 
      left_join(pubmed_data, by = "gene") %>% 
      mutate(count = ifelse(count == 0, NA, count)) %>% 
      arrange(distance) %>% 
      mutate(gene = fct_inorder(gene))
    
    p1 <- ggplot(thisFig2, aes()) +
      geom_segment(aes(x = gene, xend = gene, y = 0, yend = distance), color = "lightblue", size = 4) +
      labs(x = "", y = paste("Distance from\n", zoom_on_gene, "node")) +
      theme_bw(base_size = 14) +
      theme(legend.position = "none",
            plot.margin = margin(0,1,0,1),
            panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_blank())
  
    p2 <- ggplot(thisFig2, aes(x = gene, y = term, fill = count)) +
      geom_tile(color = "black") +
      scale_fill_gradient(low = "salmon", high = "red", na.value = "white") +
      labs(x = "", y = "") +
      theme_bw(base_size = 16) +
      theme(legend.position = "none",
            plot.margin = margin(0,1,0,1),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    thisFig3 <- thisFig2 %>% 
      select(gene) %>% 
      mutate(Validation = gene %in% validated_genes,
            `Non pan-essential` = !(gene %in% Achilles_common_essential)) %>% #citation: doi: https://doi.org/10.1101/720243
      gather(-gene, key = "parameter", value = "value")
    
    p3 <- ggplot(thisFig3, aes(x = gene, y = parameter)) +
      geom_point(aes(color = value), size = 4) + 
      scale_color_manual(values = c("white", "yellowgreen")) +
      labs(x = "", y = "") +
      theme_bw(base_size = 16) +
      theme(legend.position = "none",
            plot.margin = margin(0,1,0,1),
            panel.border = element_blank(),
            axis.text = element_text(color = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
  
    ggarrange(p1 + rremove("x.text"), ggplot() + theme_void(),
              p3, ggplot() + theme_void(),
              p2, ncol = 1, align = "v", heights = c(0.7, -0.08, 0.3, -0.08, 1))
    ggsave(paste0("Distance from ", zoom_on_gene," on UMAP.tiff"), dpi = 300, height = 4, width = 10)
  }
  
  #create excel table of the targets
  if(excel) {
    library(rentrez)
    library(openxlsx)
    
    Enrichr_result1 <- EnrichrFun(genes = node_neighbors$gene, "GO_Molecular_Function_2018", "test")
    Enrichr_result2 <- EnrichrFun(genes = node_neighbors$gene, "GO_Biological_Process_2018", "test")
    Enrichr_result3 <- EnrichrFun(genes = node_neighbors$gene, "KEGG_2019_Human", "test")
    Enrichr_result4 <- EnrichrFun(genes = node_neighbors$gene, "WikiPathway_2023_Human", "test")
    
    set_entrez_key("8f89599695f96f6a5ff245741dfe7e02d808")
    get_entrez_info <- function(gene) {
      cat("\r", "getting info for: ", gene, "                  ")
      Sys.sleep(0.1)
      search_term <- paste0(gene, "[GENE] AND Homo sapiens[ORGN]")
      search <- entrez_search(db = "gene", term = search_term)
      retrieve <- entrez_summary(db = "gene", id = search$ids[1])
    }
    
    this_Genelist <- node_neighbors %>% 
      left_join(Enrichr_result1$`long results`, by = c("gene" = "Gene")) %>% 
      left_join(Enrichr_result2$`long results`, by = c("gene" = "Gene")) %>% 
      left_join(Enrichr_result3$`long results`, by = c("gene" = "Gene")) %>% 
      dplyr::rename("KEGG" = "Terms", "GO-Molecular Function" = "Terms.x", "GO-Biological process" = "Terms.y") %>% 
      mutate(info = purrr::map(gene, get_entrez_info)) %>% 
      mutate(description = map_chr(info, "description"),
             summary = map_chr(info, "summary")) 
    
    this_Genelist_1 <- as_tibble(this_Genelist) %>% 
      select(-any_of(Geneshot_terms)) %>% 
      left_join(pubmed_data %>% spread(key = "term", value = "count"), by = "gene") %>% 
      mutate(across(any_of(Geneshot_terms), ~ifelse(.x == 0, NA, .x))) %>% 
      mutate(In_top_100 = ifelse(gene %in% top_ranking_genes, "yes", ""), 
             Is_common_essential = ifelse(gene %in% Achilles_common_essential, "yes", ""),
             Validated = ifelse(gene %in% validated_genes, "yes", "")) %>% 
      select(gene, contains("distance"), description, summary, 
             "Best Predictors" = In_top_100, Is_common_essential, Validated, 
             any_of(Geneshot_terms), KEGG, `GO-Molecular Function`,
             `GO-Biological process`, contains("expression"), contains("Depscore"), any_of(selected_pert))
    
    #export the target list to excel
    wb <- createWorkbook()
    addWorksheet(wb, zoom_on_gene, gridLines = FALSE)
    writeData(wb, sheet = zoom_on_gene, x = str_c("Gene list - ", zoom_on_gene, " nearest neighbors"), startCol = 1, startRow = 1)
    writeDataTable(wb, sheet = zoom_on_gene, x = this_Genelist_1, startCol = 1, startRow = 3, 
                   tableStyle = "TableStyleLight1", bandedRows = FALSE)
    saveWorkbook(wb, file = str_c("UMAP gene list - ", zoom_on_gene, " Best predictors.xlsx"), overwrite = TRUE)
  }
  
  return(node_neighbors)
}

GCH1_node <- nodes_nearest_neighbors("GCH1", distance_plot = TRUE)
GPX4_node <- nodes_nearest_neighbors("GPX4", distance_plot = TRUE)

#-----------------------------------------------------------------------------------------
# validate independence of UMAP random seed
#-----------------------------------------------------------------------------------------

set.seed(123)
random_seeds <- round(runif(n = 500, min = 1, max = 10000), 0)

collect_results <- list()
for(i in seq_along(random_seeds)) {
  print(str_c("Iteration ", i, ", seed ", random_seeds[i]))
  thisDF.map <- umap(thisDF.data, preserve.seed = TRUE, random_state = random_seeds[i], min_dist = 0.05)
  result <- cbind(thisDF.labels, thisDF.map$layout) %>% mutate(zoom = NA)
  
  zoom_on_gene <- "GCH1"
  coord.x <- result[result$gene == zoom_on_gene, ]$`1`
  coord.y <- result[result$gene == zoom_on_gene, ]$`2`
  GCH1_gene <- result %>% 
    mutate(distance = sqrt((`1` - coord.x)^2 + (`2` - coord.y)^2)) %>% 
    select(gene, distance) %>% 
    mutate(rank = min_rank(distance)) %>% 
    filter(gene %in% c("CANT1", "FRAT1", "DCXR", "DECR1", "EBP", "NDUFV1", "TNFRSF18", "SORD")) %>% 
    mutate(iteration = i, seed = random_seeds[i], node = "GCH1")
  
  zoom_on_gene <- "GPX4"
  coord.x <- result[result$gene == zoom_on_gene, ]$`1`
  coord.y <- result[result$gene == zoom_on_gene, ]$`2`
  GPX4_gene <- result %>% 
    mutate(distance = sqrt((`1` - coord.x)^2 + (`2` - coord.y)^2)) %>% 
    select(gene, distance) %>% 
    mutate(rank = min_rank(distance)) %>% 
    filter(gene %in% c("CTBP1", "PDAP1")) %>% 
    mutate(iteration = i, seed = random_seeds[i], node = "GPX4")
  
  collect_results[[i]] <- bind_rows(GCH1_gene, GPX4_gene)
}

collect_results_tib <- purrr::reduce(collect_results, bind_rows)

#plots
thisFig <- collect_results_tib %>% 
  as_tibble() %>% 
  filter(!(gene %in% c("PDAP1", "EBP"))) %>% 
  mutate(mean_per_gene = mean(rank), .by = gene) %>% 
  mutate(mean_per_seed = mean(rank), .by = seed) %>% 
  mutate(iteration1 = reorder_within(x = iteration, by = rank, within = gene))

ggplot(thisFig, aes(x = iteration1, y = rank)) + 
  geom_col(aes(fill = rank < 50)) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  geom_text(data = thisFig %>% distinct(gene, .keep_all = TRUE),
            aes(label = str_c("Major node: ", node)), x = 50, y = 110) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  facet_wrap(~gene, scales = "free_x", ncol = 3) +
  labs(x = "500 Random UMAPs", y = "Ranked distance from major node") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

ggsave("Validate independece of UMAP random seed plot.tiff", dpi = 300, width = 8, height = 7)

#-----------------------------------------------------------------------------------------
# UMAP distances vs. predictor matrix distances (Figure S4F)
#-----------------------------------------------------------------------------------------

set.seed(42)

#create predictor matrix
heatmap <- GeneList %>% 
  select(-any_of(Geneshot_terms)) %>% 
  na.omit() %>% 
  mutate(across(where(is.numeric), scale)) %>% 
  column_to_rownames(var = "gene")

#create predictor distance matrix
distance_mat <- heatmap %>% 
  dist(method = "euclidean") %>% 
  as.matrix() %>% 
  as_tibble(rownames = "ID") %>% 
  select(ID, GCH1) %>% 
  gather(-ID, key = "node", value = "distance") %>% 
  mutate(dist_fun = "euclidean")

#compare the distance of the predictor matrix to the UMAP matrix
UMAP_distances <- result %>% 
  as_tibble() %>% 
  mutate(GCH1 = sqrt((`1` - result %>% filter(gene == "GCH1") %>% pull(`1`)) ^ 2 +
                       (`2` - result %>% filter(gene == "GCH1") %>% pull(`2`)) ^ 2)) %>% 
  select("ID" = gene, GCH1) %>% 
  gather(-ID, key = "node", value = "distance") %>% 
  mutate(dist_fun = "UMAP")

#figure for paper
thisFig <- distance_mat %>% 
  bind_rows(UMAP_distances) %>% 
  filter(node == "GCH1") %>% 
  filter(dist_fun %in% c("euclidean", "UMAP")) %>% 
  spread(key = "dist_fun", value = "distance")

ggplot(thisFig, aes(x = UMAP, y = euclidean)) +
  geom_point(color = "gray", alpha = 0.5) + 
  geom_point(data = thisFig %>% filter(ID %in% c("GCH1", "DCXR", "DECR1", "SORD", "NDUFV1", "CANT1", "FRAT1", "TNFRSF18")), 
             color = "red") +
  geom_text_repel(data = thisFig %>% filter(ID == "GCH1"), 
                  aes(label = ID), color = "red") +
  geom_smooth(method = "lm") +
  geom_smooth(data = thisFig %>% slice_min(UMAP, n = 200), method = "lm", color = "red") +
  stat_cor(color = "blue") +
  labs(x = "Distance from GCH1\nin the UMAP", 
       y = "Distance from GCH1\nin the predictors matrix") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 16, color = "black"))

ggsave("Compare predictor distance to UMAP distance.tiff", dpi = 300, width = 5, height = 5)
