# Extended Data Figure 2 - Plasma biomarkers across FTLD classes and subtypes.
#   a  Biomarkers in ADNC w/o FTLD, FTLD-tau, FTLD-TDP, FTLD-FUS, Mixed FTLD (KW + Dunn FDR)
#   b  GFAP/NfL ratio in pure FTLD-tau vs pure FTLD-TDP (Wilcoxon rank-sum)
#   c  Biomarkers across specific FTLD subtypes + Low pathology + ADNC (exploratory)

source("setup.R")

db <- readRDS(file.path(data_dir, "demo_cohort.rds"))

dunn_brackets <- function(d, group) {
  groups_present <- d %>% group_by(biomarker) %>%
    summarise(n_groups = n_distinct(.data[[group]]), .groups = "drop") %>%
    filter(n_groups >= 2) %>% pull(biomarker)
  if (length(groups_present) == 0) return(d[0, ] %>% mutate(group1 = character(), group2 = character(),
                                                            p.adj = numeric(), p.adj.label = character(), y.position = numeric()))
  dn <- d %>% filter(biomarker %in% groups_present) %>% group_by(biomarker) %>%
    dunn_test(as.formula(paste("value ~", group)), p.adjust.method = "fdr") %>%
    filter(p.adj < 0.05) %>% mutate(p.adj.label = get_stars(p.adj)) %>% ungroup()
  d %>% group_by(biomarker) %>% summarise(max_val = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(dn, by = "biomarker") %>% group_by(biomarker) %>%
    mutate(y.position = max_val * (1.05 + (row_number() - 1) * 0.12)) %>%
    ungroup() %>% filter(!is.na(p.adj))
}

#Panel a: biomarkers by FTLD molecular class 
db_a <- db %>%
  mutate(ftld_class = case_when(
    ftld_any == 1 & (ftld_tau + ftld_tdp + ftld_fus) > 1 ~ "Mixed FTLD",
    ftld_tau == 1 ~ "FTLD-tau",
    ftld_tdp == 1 ~ "FTLD-TDP",
    ftld_fus == 1 ~ "FTLD-FUS",
    ad %in% c("Int ADNC", "Sev ADNC") & !FTLD_present ~ "ADNC w/o FTLD",
    TRUE ~ NA_character_),
    ftld_class = factor(ftld_class,
                        levels = c("ADNC w/o FTLD", "FTLD-tau", "FTLD-TDP", "FTLD-FUS", "Mixed FTLD")))

colors_a <- c("ADNC w/o FTLD" = "#F08870", "FTLD-tau" = colors_ftld[["FTLD-tau"]],
              "FTLD-TDP" = colors_ftld[["FTLD-TDP"]], "FTLD-FUS" = colors_ftld[["FTLD-FUS"]],
              "Mixed FTLD" = "#A0522D")

box_a <- db_a %>% filter(!is.na(ftld_class)) %>%
  select(ID, ftld_class, all_of(biomarker_vars)) %>%
  pivot_longer(all_of(biomarker_vars), names_to = "biomarker", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(biomarker = factor(biomarker, biomarker_vars, biomarker_labels),
         ftld_class = droplevels(ftld_class))

brackets_a <- dunn_brackets(box_a, "ftld_class")

p_a <- ggplot(box_a) +
  geom_boxplot(aes(ftld_class, value, fill = ftld_class), outlier.shape = NA, alpha = 0.75, linewidth = 0.4, width = 0.7) +
  geom_quasirandom(aes(ftld_class, value), size = 0.5, shape = 21, color = "black", stroke = 0.2, alpha = 0.7, width = 0.2) +
  scale_fill_manual(values = colors_a) +
  facet_wrap(~ biomarker, nrow = 1, scales = "free_y") +
  theme_manuscript(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Concentration (pg/mL)", title = "Ext Data Fig 2a: Plasma biomarkers by FTLD class")
if (nrow(brackets_a) > 0)
  p_a <- p_a + stat_pvalue_manual(brackets_a, label = "p.adj.label", tip.length = 0.01, size = 2.5, bracket.size = 0.3)
ggsave(file.path(out_dir, "ExtData2a_biomarkers_by_ftld_class.pdf"), p_a, width = 7, height = 3, dpi = 300)

# Panel b: GFAP/NfL ratio, pure FTLD-tau vs pure FTLD-TDP
db_b <- db %>%
  filter(!is.na(GFAP), !is.na(NfL), NfL > 0) %>%
  mutate(GFAP_NfL_ratio = GFAP / NfL,
         ftld_molecular = case_when(
           ftld_tau == 1 & ftld_tdp == 0 ~ "FTLD-tau",
           ftld_tdp == 1 & ftld_tau == 0 ~ "FTLD-TDP",
           TRUE ~ NA_character_),
         ftld_molecular = factor(ftld_molecular, levels = c("FTLD-tau", "FTLD-TDP"))) %>%
  filter(!is.na(ftld_molecular)) %>% mutate(ftld_molecular = droplevels(ftld_molecular))

colors_b <- c("FTLD-tau" = colors_ftld[["FTLD-tau"]], "FTLD-TDP" = colors_ftld[["FTLD-TDP"]])

p_b <- ggplot(db_b) +
  geom_boxplot(aes(ftld_molecular, GFAP_NfL_ratio, fill = ftld_molecular), outlier.shape = NA, alpha = 0.75, linewidth = 0.4, width = 0.6) +
  geom_quasirandom(aes(ftld_molecular, GFAP_NfL_ratio), size = 0.5, shape = 21, color = "black", stroke = 0.2, alpha = 0.7, width = 0.2) +
  scale_fill_manual(values = colors_b) +
  theme_manuscript(9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = NULL, y = "GFAP / NfL Ratio", title = "Ext Data Fig 2b:\nGFAP/NfL ratio by FTLD subtype")

if (n_distinct(db_b$ftld_molecular) == 2 &&
    all(table(db_b$ftld_molecular) >= 2)) {
  wilcox_bracket <- db_b %>% wilcox_test(GFAP_NfL_ratio ~ ftld_molecular) %>%
    mutate(p.label = get_stars(p),
           y.position = max(db_b$GFAP_NfL_ratio, na.rm = TRUE) * 1.10,
           p.adj.label = ifelse(p < 0.05, paste0("p = ", format_pvalue(p), " ", p.label),
                                paste0("p = ", format_pvalue(p), " ns")))
  p_b <- p_b + stat_pvalue_manual(wilcox_bracket, label = "p.adj.label",
                                  tip.length = 0.02, size = 3, bracket.size = 0.4)
}
ggsave(file.path(out_dir, "ExtData2b_gfap_nfl_ratio.pdf"), p_b, width = 2.5, height = 3, dpi = 300)

# Panel c: specific FTLD subtypes 
db_c <- db %>%
  mutate(ftld_subtype_detail = case_when(
    ftld_tau == 1 & ftld_tdp == 1 ~ "Mixed FTLD",
    ftld_tau_specific == "PSP" ~ "PSP",
    ftld_tau_specific == "CBD" ~ "CBD",
    ftld_tau_specific == "Pick's Disease" ~ "Pick's Disease",
    ftld_tau_specific == "AGD" ~ "AGD",
    ftld_tau == 1 ~ "FTLD-tau NOS",
    ftld_tdp == 1 ~ "FTLD-TDP",
    ftld_fus == 1 ~ "FTLD-FUS",
    Overall_Group == "Low Path" ~ "Low pathology",
    Overall_Group %in% c("Int ADNC", "Sev ADNC") ~ "ADNC",
    TRUE ~ NA_character_),
    ftld_subtype_detail = factor(ftld_subtype_detail,
                                 levels = c("Low pathology", "PSP", "CBD", "Pick's Disease", "AGD",
                                            "FTLD-tau NOS", "FTLD-FUS", "FTLD-TDP", "Mixed FTLD", "ADNC")))

group_counts <- table(db_c$ftld_subtype_detail)
sufficient_groups <- names(group_counts[group_counts >= 3])

colors_c <- c("Low pathology" = "#E8E8E8", "PSP" = colors_tauopathy[["PSP"]], "CBD" = colors_tauopathy[["CBD"]],
              "Pick's Disease" = colors_tauopathy[["Pick's Disease"]], "AGD" = colors_tauopathy[["AGD"]],
              "FTLD-tau NOS" = colors_tauopathy[["Other Tauopathy"]], "FTLD-FUS" = colors_ftld[["FTLD-FUS"]],
              "FTLD-TDP" = colors_ftld[["FTLD-TDP"]], "Mixed FTLD" = "#A0522D", "ADNC" = "#F08870")

box_c <- db_c %>% filter(ftld_subtype_detail %in% sufficient_groups) %>%
  select(ID, ftld_subtype_detail, all_of(biomarker_vars)) %>%
  pivot_longer(all_of(biomarker_vars), names_to = "biomarker", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(biomarker = factor(biomarker, biomarker_vars, biomarker_labels),
         ftld_subtype_detail = droplevels(ftld_subtype_detail))

p_c <- ggplot(box_c) +
  geom_boxplot(aes(ftld_subtype_detail, value, fill = ftld_subtype_detail), outlier.shape = NA, alpha = 0.75, linewidth = 0.4, width = 0.7) +
  geom_quasirandom(aes(ftld_subtype_detail, value), size = 0.5, shape = 21, color = "black", stroke = 0.2, alpha = 0.7, width = 0.2) +
  scale_fill_manual(values = colors_c) +
  facet_wrap(~ biomarker, nrow = 1, scales = "free_y") +
  theme_manuscript(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Concentration (pg/mL)", title = "Ext Data Fig 2c: Plasma biomarkers by FTLD subtype")
ggsave(file.path(out_dir, "ExtData2c_biomarkers_by_subtype.pdf"), p_c, width = 9, height = 3, dpi = 300)

cat("Extended Data Figure 2 panels written to", out_dir, "\n")
