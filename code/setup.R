suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(ggbeeswarm)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(splines)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(broom)
  library(broom.mixed)
  library(survival)
  library(survminer)
})

data_dir <- "../demo_data"
out_dir  <- "../demo_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

biomarker_vars   <- c("ptau217", "ptau181", "GFAP", "NfL")
biomarker_labels <- c(ptau217 = "p-tau217", ptau181 = "p-tau181", GFAP = "GFAP", NfL = "NfL")

theme_manuscript <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(size = base_size - 1, color = "black"),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.4, "cm"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "gray40"),
      panel.spacing = unit(0.3, "lines")
    )
}
theme_set(theme_manuscript())

colors_adnc <- c("Low ADNC" = "#FEE5D9", "Int ADNC" = "#F08870", "Sev ADNC" = "#CB181D")
colors_biomarkers <- c("p-tau217" = "#B2182B", "p-tau181" = "#E66101", "GFAP" = "#1B7837", "NfL" = "#762A83")
colors_lbd  <- c("Absent" = "white", "Brainstem" = "#A4D8E0", "Limbic" = "#5BBEC9", "Neocortical" = "#2BA7B8")
colors_late <- c("Absent" = "white", "LATE Amygdala" = "#C8B3EF", "LATE Limbic" = "#B79CE8", "LATE Neocortical" = "#6F42C1")
colors_cerad <- c("None" = "#F0F0F0", "Sparse" = "#FDB863", "Moderate" = "#E08214", "Frequent" = "#B35806")
colors_late_copathology <- c("Low Path" = "#E8E8E8", "LATE" = "#E9E0F7", "ADNC" = "#F08870", "ADNC + LATE" = "#B79CE8")
colors_lbd_copathology  <- c("Low Path" = "#E8E8E8", "LBD" = "#9AD4D6", "ADNC" = "#F08870", "ADNC + LBD" = "#2BA7B8")
colors_ftld_copathology <- c("Low Path" = "#E8E8E8", "FTLD" = "#8C564B", "ADNC" = "#F08870", "ADNC + FTLD" = "#6F4E37")
colors_overall_pathology <- c("Low Path" = "#E8E8E8", "Int ADNC" = "#FFB366", "Sev ADNC" = "#FF8533",
                              "FTLD" = "#8C564B", "LBD" = "#2BA7B8", "LATE" = "#9B59B6",
                              "Mixed ADNC" = "#BCBD22", "Other" = "#969696", "Unclassified" = "white")
colors_overall_group2 <- c("Low Path" = "#E8E8E8", "Other" = "#B8B8B8", "ADNC" = "#F08870", "Mixed ADNC" = "#B85450")
colors_part_group <- c("Control" = "#E8E8E8", "PART" = "#7B9FB8", "Int ADNC" = "#FFB366", "Sev ADNC" = "#FF8533")
colors_age_group <- c("<70" = "#016c59", "70-80" = "#1c9099", "80-90" = "#67a9cf", "90+" = "#bdc9e1")
colors_time_to_death_data <- c("0-5 years" = "#B35806", "5-10 years" = "#E08214", "10-15 years" = "#FDB863", "15+ years" = "#fdead2")
colors_mmse_gradient <- c("#67000D", "#A50F15", "#CB181D", "#FB6A4A", "#FCBBA1")
colors_apoe <- c("APOE ε4 -" = "#91BFDB", "APOE ε4 +" = "#FC8D59")
colors_sex <- c("Male" = "#a6dba0", "Female" = "#c2a5cf")
colors_binary <- c("Absent" = "white", "Present" = "#2E5090", "Missing" = "#E8E8E8")
colors_ftld <- c("FTLD-tau" = "#8C564B", "FTLD-TDP" = "#C49C94", "FTLD-FUS" = "#E7969C", "FTLD-UPS" = "#7B4173", "FTLD-NOS" = "#A55194")
colors_tauopathy <- c("Pick's Disease" = "#8C6D62", "CBD" = "#A0826D", "PSP" = "#B39A7C", "AGD" = "#C9B18A", "Other Tauopathy" = "#6F5F5A")
colors_copathology_roc <- c("No co-pathology" = "#F08870", "With LATE" = "#B79CE8", "With HS" = "#FDAE6B",
                            "With LBD" = "#2BA7B8", "With FTLD" = "#8C564B", "With Other" = "#7F7F7F")

format_pvalue <- function(p, threshold = 0.001) {
  if (is.na(p)) return("NA")
  if (p < threshold) return(paste0("<", threshold))
  sprintf("%.3f", p)
}

get_stars <- function(p) {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

calc_roc_metrics <- function(data, outcome, predictor, direction = "<") {
  d <- data %>% filter(!is.na(.data[[outcome]]), !is.na(.data[[predictor]]))
  if (nrow(d) < 10) return(NULL)
  ro <- roc(d[[outcome]], d[[predictor]], direction = direction, quiet = TRUE)
  ci <- ci.auc(ro, conf.level = 0.95)
  yc <- coords(ro, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity", "accuracy", "ppv", "npv"))
  data.frame(N = nrow(d), AUC = as.numeric(auc(ro)), AUC_Lower = ci[1], AUC_Upper = ci[3],
             Cutoff = yc$threshold, Sensitivity = yc$sensitivity, Specificity = yc$specificity,
             Accuracy = yc$accuracy, PPV = yc$ppv, NPV = yc$npv,
             Youden = yc$sensitivity + yc$specificity - 1)
}

add_overall_group2 <- function(data) {
  data %>% mutate(
    Overall_Group2 = case_when(
      Overall_Group == "Low Path" ~ "Low Path",
      Overall_Group %in% c("FTLD", "LBD", "Other", "LATE") ~ "Other",
      Overall_Group %in% c("Int ADNC", "Sev ADNC") ~ "ADNC",
      Overall_Group == "Mixed ADNC" ~ "Mixed ADNC",
      TRUE ~ NA_character_),
    Overall_Group2 = factor(Overall_Group2, levels = c("Low Path", "Other", "ADNC", "Mixed ADNC")))
}

add_adnc_binary <- function(data) {
  data %>% mutate(ad_binary = case_when(ad == "Low ADNC" ~ 0, ad %in% c("Int ADNC", "Sev ADNC") ~ 1, TRUE ~ NA_real_))
}
