#this script takes the dose response curve data, outputs the DRC curves and caluclate the AUCs, as shown in figure 2.

library(tidyverse)
library(drc)
library(ggsignif)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(PharmacoGx)

#------------------------------------------------------------------------------------------------------
# Plot DRC plots and calculate AUCs
#------------------------------------------------------------------------------------------------------

#function definition
Plot_DRC <- function(DRC_data, FIN) {
  
  #DRC plots
  DRC_plot <- DRC_data %>% 
    dplyr::select(-conc) %>% 
    gather(-cell, -logConc, key = "rep", value = "viability") %>% 
    group_by(cell, logConc) %>% 
    get_summary_stats() %>% 
    filter(!(cell %in% c("HCC1143", "BT474"))) %>% 
    mutate(cell = factor(cell, levels = c("MCF7", "T47D", "JIMT1", "SKBR3", "MDAMB453", "MDAMB468", "BT20", 
                                          "HCC1937", "HCC70", "MDAMB231", "SUM159PT", "HS578T", "MDAMB436", "BT549")))

  ggplot(DRC_plot, aes(x = logConc, y = mean, color = cell)) + 
    geom_point() +
    geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd), size = 0.5, width = 0.1) +
    stat_smooth(se = FALSE, method = "drm", 
                method.args = list(fct=L.4(fixed = c(NA, NA, 100, NA)))) +
    scale_color_manual(values = c("grey80", "grey60", "grey40", "grey25", "grey10", "red", "tan3", 
                                        "green", "pink", "steelblue1", "blue", "yellowgreen",
                                        "darkviolet", "aquamarine4")) +
    scale_y_continuous(breaks = seq(0, 120, 20), limits = c(0, 120)) +
    theme_classic(base_size = 16) +
    theme(axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(color = "black", size = 14), 
          legend.title = element_text(size = 11),
          legend.text = element_text(color = "black", size = 12)) +
    labs(x = paste0("Log [", FIN, "], M"), y = "Cell viability (%)", title = "", color = "")
  
  ggsave(paste0("DRC plot for ", FIN, " (Varsha).tiff"), dpi = 300, height = 4, width = 6)

  #calculate parameters
  DRC_measures <- DRC_data %>% 
    dplyr::select(-conc) %>% 
    gather(-cell, -logConc, key = "rep", value = "viability")

  cells <- DRC_data %>% distinct(cell) %>% pull(cell)

  parameters <- list()
  for(i in seq_along(cells)) {
    thisDF <- DRC_data %>% 
      filter(cell == cells[i])
    n_reps <- thisDF %>% summarise(across(-c(cell, conc, logConc), ~ sum(is.na(.)))) %>% 
      gather() %>% filter(value < 7) %>% pull(key) %>% as.numeric()
    thisDFlong <- thisDF %>% 
      gather(-cell, -conc, -logConc, key = "rep", value = "viability") %>% 
      na.omit()
    drmModel <- drm(viability ~ logConc, data = thisDFlong, curveid = rep, 
                    fct = L.4(names = c("slope", "lower limit", "upper limit", "EC50"),
                              fixed = c(NA, 0, 100, NA)))
    result <- drmModel[["parmMat"]]
    auc <- c()
    for(j in seq_along(n_reps)) {
      auc[j] <- 100 - computeAUC(thisDF$logConc, thisDF[[j+3]], conc_as_log = TRUE, viability_as_pct = TRUE, area.type = "Actual")
    }
    parameters[[i]] <- tibble(
      "inducer" = FIN,
      "cell" = cells[i],
      "repeat" = colnames(result),
      auc = auc,
      logEC50 = (result[2,]))
  }

  parameters_tb <- purrr::reduce(parameters, bind_rows)
  return(parameters_tb)
}

DRC_data <- read_csv("Data/drc curve - FIN56.csv")
FIN56_AUCs <- Plot_DRC(DRC_data, FIN = "FIN56")

DRC_data <- read_csv("Data/drc curve - Erastin.csv")
Erastin_AUCs <- Plot_DRC(DRC_data, FIN = "Erastin")

DRC_data <- read_csv("Data/drc curve - RSL3.csv")
RSL3_AUCs <- Plot_DRC(DRC_data, FIN = "RSL3")

#gather AUCS
FIN_AUCs <- bind_rows(FIN56_AUCs, Erastin_AUCs, RSL3_AUCs) %>% 
  dplyr::select(-`repeat`, -logEC50) %>% 
  group_by(inducer, cell) %>% 
  get_summary_stats()

#------------------------------------------------------------------------------------------------------
# AUCs analysis
#------------------------------------------------------------------------------------------------------

#plot AUCs by TNBC status
PAM50 <- read_csv("Data/PAM50 info.csv")

thisFig <- FIN_AUCs %>% 
  left_join(PAM50, by = c("cell" = "ID"))
  
ggplot(thisFig, aes(x = `TNBC status`, y = mean, fill = `TNBC status`)) +
  geom_boxplot(outlier.shape = NA, position=position_dodge(1)) +
  geom_dotplot(dotsize = 0.4, binaxis = 'y', stackdir = 'center',
               position = position_dodge(1), show.legend = FALSE) +
  geom_signif(comparisons = list(c("TNBC", "Non-TNBC")),
              test = "t.test", test.args = list(var.equal = TRUE),
              map_signif_level = TRUE, vjust = 0.3) +
  scale_fill_manual(values = c("gray", "salmon")) +
  scale_y_continuous(limits = c(30, 105), breaks = c(40, 60, 80, 100)) +
  #coord_flip() +
  facet_wrap(~ inducer, nrow = 1, strip.position = "bottom") +
  labs(x = "", y = "AUC", fill = "") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        strip.placement = "outside",
        legend.position = "right",
        legend.text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave("AUCs summary by TNBC.tiff", dpi = 300, width = 5, height = 3)

#plot AUCs heatmap
thisFig <- FIN_AUCs %>% 
  left_join(PAM50, by = c("cell" = "ID")) %>% 
  arrange(mean) %>% 
  mutate(cell = fct_inorder(cell)) %>% 
  complete(inducer, cell)

cols <- c("Non-TNBC" = "black", "Mesenchymal" = "red", "Basal-like" = "blue", "UCL" = "blue")

ggplot(thisFig, aes(x = cell, y = inducer, fill = mean)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(mean, 0))) +
  scale_fill_gradient2(low = "red", high = "black", na.value = "gray", midpoint = 77) +
  labs(x = "", y = "", fill = "AUC") +
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 3.8, ymax = 4.5, fill = "salmon") + 
  annotate("rect", xmin = 9.5, xmax = 16.5, ymin = 3.8, ymax = 4.5, fill = "gray") + 
  annotate("text", x = 4.5, y = 4.2, label = "TNBC", size = 5, color = "black") + 
  annotate("text", x = 13, y = 4.2, label = "Non-TNBC", size = 5, color = "black") + 
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("AUCs summary - boxplot.tiff", dpi = 300, height = 3, width = 8)
