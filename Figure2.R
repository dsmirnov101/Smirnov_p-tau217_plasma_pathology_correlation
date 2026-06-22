# Figure 1 - Plasma biomarker performance for autopsy-confirmed ADNC.
#   a  Participant-level co-pathology / biomarker heatmap
#   b  Biomarker concentrations by ADNC severity
#   c  ROC curves per biomarker (Intermediate/High vs Low ADNC)
#   d  ROC curves for multi-marker logistic models
#   e  Biomarkers across Braak 0-II / PART / Intermediate / High ADNC groups

source("setup.R")

db <- readRDS(file.path(data_dir, "demo_cohort.rds")) %>% add_overall_group2() %>% add_adnc_binary()

# Panel b: biomarker concentrations by ADNC severity
box_data <- db %>%
  filter(!is.na(ad)) %>%
  select(ID, ad, all_of(biomarker_vars)) %>%
  pivot_longer(all_of(biomarker_vars), names_to = "biomarker", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(biomarker = factor(biomarker, biomarker_vars, biomarker_labels),
         ad = factor(ad, c("Low ADNC", "Int ADNC", "Sev ADNC")))

dunn_brackets <- function(d, group) {
  dn <- d %>% group_by(biomarker) %>%
    dunn_test(as.formula(paste("value ~", group)), p.adjust.method = "fdr") %>%
    filter(p.adj < 0.05) %>% mutate(p.adj.label = get_stars(p.adj)) %>% ungroup()
  d %>% group_by(biomarker) %>% summarise(max_val = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(dn, by = "biomarker") %>% group_by(biomarker) %>%
    mutate(y.position = max_val * (1.05 + (row_number() - 1) * 0.12)) %>%
    ungroup() %>% filter(!is.na(p.adj))
}
box_brackets <- dunn_brackets(box_data, "ad")

p_b <- ggplot(box_data) +
  geom_boxplot(aes(ad, value, fill = ad), outlier.shape = NA, alpha = 0.75, linewidth = 0.4, width = 0.7) +
  geom_quasirandom(aes(ad, value), size = 0.5, shape = 21, color = "black", stroke = 0.2, alpha = 0.7, width = 0.2) +
  scale_fill_manual(values = colors_adnc) +
  facet_wrap(~ biomarker, nrow = 1, scales = "free_y") +
  theme_manuscript(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Concentration (pg/mL)", title = "Figure 1b: Plasma biomarkers by ADNC severity") +
  stat_pvalue_manual(box_brackets, label = "p.adj.label", tip.length = 0.01, size = 2.5, bracket.size = 0.3)
ggsave(file.path(out_dir, "Figure1b_biomarkers_by_adnc.pdf"), p_b, width = 7, height = 3, dpi = 300)

# Panel c: ROC curves per biomarker
roc_data <- db %>% filter(!is.na(ad_binary))
roc_objs <- list()
roc_curves <- bind_rows(lapply(seq_along(biomarker_vars), function(i) {
  b <- biomarker_vars[i]; d <- roc_data %>% filter(!is.na(.data[[b]]))
  ro <- roc(d$ad_binary, d[[b]], direction = "<", quiet = TRUE)
  roc_objs[[unname(biomarker_labels[i])]] <<- ro
  data.frame(sensitivity = ro$sensitivities, specificity = ro$specificities,
             biomarker = unname(biomarker_labels[i]), auc = as.numeric(auc(ro)))
})) %>%
  mutate(fpr = 1 - specificity,
         label = sprintf("%s (AUC=%.3f)", biomarker, auc),
         label = factor(label, levels = unique(label[order(match(biomarker, biomarker_labels))])))

p_c <- ggplot(roc_curves, aes(fpr, sensitivity, color = label)) +
  geom_line(linewidth = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_color_manual(values = unname(colors_biomarkers[biomarker_labels])) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_manuscript(9) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.position = c(0.70, 0.22), legend.background = element_blank()) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Figure 1c: Biomarker ROC for ADNC")
ggsave(file.path(out_dir, "Figure1c_roc_all_biomarkers.pdf"), p_c, width = 4, height = 4, dpi = 300)

# Panel d: multi-marker logistic models
model_data <- db %>%
  mutate(sex_numeric = as.numeric(factor(Sex)) - 1, apoe4_carrier = as.integer(apoe4 >= 1)) %>%
  select(ad_binary, ptau217, ptau181, GFAP, NfL, NPDAGE, sex_numeric, apoe4_carrier) %>%
  filter(complete.cases(.))
models <- list(
  "p-tau217 + NfL" = ad_binary ~ ptau217 + NfL,
  "All 4 biomarkers" = ad_binary ~ ptau217 + ptau181 + GFAP + NfL,
  "All 4 + Age + Sex" = ad_binary ~ ptau217 + ptau181 + GFAP + NfL + NPDAGE + sex_numeric,
  "All 4 + Age + Sex + APOE4" = ad_binary ~ ptau217 + ptau181 + GFAP + NfL + NPDAGE + sex_numeric + apoe4_carrier)
ptau217_alone <- roc(model_data$ad_binary, model_data$ptau217, direction = "<", quiet = TRUE)
model_curves <- bind_rows(
  data.frame(sensitivity = ptau217_alone$sensitivities, specificity = ptau217_alone$specificities,
             model = "p-tau217 alone", auc = as.numeric(auc(ptau217_alone))),
  bind_rows(lapply(names(models), function(nm) {
    fit <- glm(models[[nm]], data = model_data, family = binomial)
    ro <- roc(model_data$ad_binary, predict(fit, type = "response"), direction = "<", quiet = TRUE)
    data.frame(sensitivity = ro$sensitivities, specificity = ro$specificities,
               model = nm, auc = as.numeric(auc(ro))) }))) %>%
  mutate(fpr = 1 - specificity,
         label = sprintf("%s (AUC=%.3f)", model, auc),
         label = factor(label, levels = unique(label[order(match(model, c("p-tau217 alone", names(models))))])))
colors_model <- c("#B2182B", "#762A83", "#FFA600", "#01557b", "#1A9850")
p_d <- ggplot(model_curves, aes(fpr, sensitivity, color = label)) +
  geom_line(linewidth = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_color_manual(values = colors_model) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_manuscript(9) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7), legend.position = "right") +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Figure 1d: Multi-marker models")
ggsave(file.path(out_dir, "Figure1d_multimarker_models.pdf"), p_d, width = 6, height = 3, dpi = 300)

#  Panel e: Braak 0-II / PART / Intermediate / High ADNC

part_data <- db %>%
  mutate(PART_group = factor(case_when(
    ad == "Low ADNC" & Braak_grouped %in% c("0", "I-II") ~ "Control",
    Braak_grouped == "III-IV" & NEUR %in% c("None", "Sparse") ~ "PART",
    Braak_grouped == "III-IV" & NEUR %in% c("Moderate", "Frequent") ~ "Int ADNC",
    ad == "Sev ADNC" ~ "Sev ADNC", TRUE ~ NA_character_),
    levels = c("Control", "PART", "Int ADNC", "Sev ADNC"))) %>%
  filter(!is.na(PART_group)) %>%
  select(ID, PART_group, all_of(biomarker_vars)) %>%
  pivot_longer(all_of(biomarker_vars), names_to = "biomarker", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(biomarker = factor(biomarker, biomarker_vars, biomarker_labels))
part_brackets <- dunn_brackets(part_data, "PART_group")
p_e <- ggplot(part_data) +
  geom_boxplot(aes(PART_group, value, fill = PART_group), outlier.shape = NA, alpha = 0.75, linewidth = 0.4, width = 0.7) +
  geom_quasirandom(aes(PART_group, value), size = 0.5, shape = 21, color = "black", stroke = 0.2, alpha = 0.7, width = 0.2) +
  scale_fill_manual(values = colors_part_group) +
  facet_wrap(~ biomarker, nrow = 1, scales = "free_y") +
  theme_manuscript(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Concentration (pg/mL)", title = "Figure 1e: Biomarkers by Braak 0-II / PART / Int / High ADNC")
if (nrow(part_brackets) > 0)
  p_e <- p_e + stat_pvalue_manual(part_brackets, label = "p.adj.label", tip.length = 0.01, size = 2.5, bracket.size = 0.3)
ggsave(file.path(out_dir, "Figure1e_biomarkers_by_part.pdf"), p_e, width = 7, height = 3, dpi = 300)

#Panel a: participant-level heatmap 
dat <- db %>% filter(!is.na(NPDAGE), !is.na(ptau181) | !is.na(ptau217) | !is.na(GFAP) | !is.na(NfL)) %>%
  mutate(
    flag_ftld_tau_nos = as.integer(ftld_tau %in% 1 & !flag_picks %in% 1 & !flag_psp %in% 1 & !flag_cbd %in% 1 & !flag_agd %in% 1),
    flag_other_nontau = as.integer(flag_ftld_fus %in% 1 | flag_nifid %in% 1 | flag_bibd %in% 1),
    flag_als_combined = pmax(flag_als, flag_tdp43_spinal, na.rm = TRUE),
    flag_tbi = as.integer(flag_tbi_chronic %in% 1 | flag_tbi_acute %in% 1))
for (b in biomarker_vars)
  dat[[paste0(b, "_rank")]] <- rank(dat[[b]], na.last = "keep", ties.method = "average") / sum(!is.na(dat[[b]]))

config <- tribble(
  ~Category, ~Pathology, ~Pathology_label, ~row_order,
  "1. AD Pathology", "adnc_pathology", "ADNC", 1, "1. AD Pathology", "part", "PART", 2,
  "2. Synucleinopathies", "lbd_stage", "Lewy Body Disease", 3,
  "3. LATE-NC", "late_stage", "LATE-NC", 4, "3. LATE-NC", "hippo_sclerosis", "Hippocampal Sclerosis", 5,
  "4. FTLD/ALS", "psp", "PSP", 6, "4. FTLD/ALS", "cbd", "CBD", 7, "4. FTLD/ALS", "picks", "Pick's Disease", 8,
  "4. FTLD/ALS", "agd", "AGD", 9, "4. FTLD/ALS", "ftld_tau_nos", "FTLD-tau NOS", 10,
  "4. FTLD/ALS", "other_nontau", "FTLD-FUS", 11, "4. FTLD/ALS", "ftld_tdp", "FTLD-TDP", 12,
  "4. FTLD/ALS", "als_combined", "ALS/MND", 13,
  "5. Biomarkers", "ptau181", "p-tau181", 14, "5. Biomarkers", "ptau217", "p-tau217", 15,
  "5. Biomarkers", "GFAP", "GFAP", 16, "5. Biomarkers", "NfL", "NfL", 17)

ord <- dat %>%
  mutate(.ad = case_when(ad == "Sev ADNC" ~ 4, ad == "Int ADNC" ~ 3, ad == "Low ADNC" ~ 1, TRUE ~ 0),
         .lbd = as.integer(lbd_stage), .late = as.integer(late_stage)) %>%
  arrange(desc(.ad), desc(.lbd), desc(.late)) %>% transmute(ID, participant_order = row_number())

special <- bind_rows(
  dat %>% transmute(ID, Category = "1. AD Pathology", Pathology = "adnc_pathology", Present = as.character(ad)),
  dat %>% transmute(ID, Category = "1. AD Pathology", Pathology = "part",
                    Present = if_else(Braak_grouped == "III-IV" & NEUR %in% c("None", "Sparse"), "Present", "Absent")),
  dat %>% transmute(ID, Category = "2. Synucleinopathies", Pathology = "lbd_stage",
                    Present = case_when(lbd_stage == "Neocortical" ~ "LBD Neocortical", lbd_stage == "Limbic" ~ "LBD Limbic",
                                        lbd_stage == "Brainstem" ~ "LBD Brainstem", TRUE ~ "Absent")),
  dat %>% transmute(ID, Category = "3. LATE-NC", Pathology = "late_stage", Present = as.character(late_stage)))
binary_flags <- dat %>% select(ID, starts_with("flag_")) %>%
  pivot_longer(-ID, names_to = "Pathology", values_to = "val") %>%
  mutate(Pathology = str_remove(Pathology, "^flag_")) %>%
  inner_join(config %>% select(Category, Pathology), by = "Pathology") %>%
  mutate(Present = if_else(val == 1, "Present", "Absent")) %>%
  group_by(Pathology) %>% filter(sum(Present == "Present") > 0) %>% ungroup() %>%
  select(ID, Category, Pathology, Present)
bio_rows <- dat %>% select(ID, ptau181_rank, ptau217_rank, GFAP_rank, NfL_rank) %>%
  pivot_longer(-ID, names_to = "Pathology", values_to = "r") %>%
  mutate(Pathology = str_remove(Pathology, "_rank"), Category = "5. Biomarkers",
         Present = if_else(is.na(r), "Missing", as.character(r))) %>%
  select(ID, Category, Pathology, Present)

facet_levels <- c("1. AD Pathology", "2. Synucleinopathies", "3. LATE-NC", "4. FTLD/ALS", "5. Biomarkers")
plot_df <- bind_rows(special, binary_flags, bio_rows) %>%
  inner_join(config, by = c("Category", "Pathology")) %>%
  left_join(ord, by = "ID") %>%
  mutate(Category = factor(Category, facet_levels),
         fill_value = case_when(
           is.na(Present) ~ "Missing",
           Pathology %in% c("adnc_pathology", "lbd_stage", "late_stage") ~ Present,
           Pathology == "part" & Present == "Present" ~ "PART",
           Pathology == "hippo_sclerosis" & Present == "Present" ~ "Hippo Sclerosis",
           Category == "5. Biomarkers" & Present != "Missing" ~
             paste0("Rank_", sprintf("%03d", pmax(1L, pmin(100L, as.integer(ceiling(suppressWarnings(as.numeric(Present)) * 100)))))),
           Category == "5. Biomarkers" ~ "Missing",
           Present == "Present" ~ "Present", TRUE ~ "Absent"),
         Pathology_label = factor(Pathology_label, levels = config$Pathology_label[order(config$row_order)]))

rank_cols <- setNames(colorRampPalette(c("#5E9DC9", "#E0EFF5", "#F5F5F5", "#F49887", "#DF5649"))(100),
                      paste0("Rank_", sprintf("%03d", 1:100)))
base_cols <- c("Absent" = "white", "Missing" = "grey90",
               "Low ADNC" = unname(colors_adnc["Low ADNC"]), "Int ADNC" = unname(colors_adnc["Int ADNC"]),
               "Sev ADNC" = unname(colors_adnc["Sev ADNC"]), "PART" = "#7B9FB8",
               "LBD Brainstem" = unname(colors_lbd["Brainstem"]), "LBD Limbic" = unname(colors_lbd["Limbic"]),
               "LBD Neocortical" = unname(colors_lbd["Neocortical"]),
               "LATE Amygdala" = unname(colors_late["LATE Amygdala"]), "LATE Limbic" = unname(colors_late["LATE Limbic"]),
               "LATE Neocortical" = unname(colors_late["LATE Neocortical"]), "Hippo Sclerosis" = "#9B4D96",
               "Present" = unname(colors_binary["Present"]))
legend_breaks <- names(base_cols)[names(base_cols) != "Missing"]

p_a <- ggplot(plot_df, aes(participant_order, fct_rev(Pathology_label))) +
  geom_tile(aes(fill = fill_value), color = NA) +
  scale_fill_manual(name = "Status", values = c(base_cols, rank_cols),
                    breaks = legend_breaks, na.value = "white") +
  facet_grid(rows = vars(Category), scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Figure 1a: Co-pathology and biomarker patterns", x = paste0("Participants (N = ", nrow(dat), ")"), y = NULL) +
  theme_minimal(8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 7, hjust = 1),
        strip.text = element_blank(), panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.4),
        legend.position = "bottom", legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 7, face = "bold"), legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
ggsave(file.path(out_dir, "Figure1a_copathology_heatmap.pdf"), p_a, width = 7, height = 4, dpi = 300)

cat("Figure 1 panels written to", out_dir, "\n")
