# Figure 3 - Robustness of p-tau217 discrimination (Int/High vs Low ADNC) across strata.
#   a  Demographic strata: full cohort, age at death, sex
#   b  Clinical strata: clinical AD vs non-AD (NACCALZD), MMSE quintiles
#   c  Genetic stratum: APOE e4 carrier status
#   d  Co-pathology strata: none / LATE / LBD / FTLD; vascular subgroups
# Per-stratum AUC + 95% CI via pROC. Cross-stratum delta-AUC via subject-level
# cluster bootstrap (B = 2000) with BH-FDR within each stratification family.

source("setup.R")

db <- readRDS(file.path(data_dir, "demo_cohort.rds")) %>% add_adnc_binary()

db <- db %>%
  mutate(
    ad = factor(ad, levels = c("Low ADNC", "Int ADNC", "Sev ADNC")),
    age_group = cut(NPDAGE, breaks = c(-Inf, 70, 80, 90, Inf),
                    labels = c("<70", "70-80", "80-90", "90+"), right = FALSE),
    Sex = case_when(Sex == 1 ~ "Male", Sex == 2 ~ "Female", TRUE ~ NA_character_),
    APOE4_status = case_when(apoe4 >= 1 ~ "APOE e4 +", apoe4 == 0 ~ "APOE e4 -",
                             TRUE ~ NA_character_),
    etiology_group = case_when(NACCALZD == 1 ~ "Clinical AD etiology",
                               NACCALZD == 0 ~ "Clinical non-AD etiology",
                               TRUE ~ NA_character_),
    ADNC_present  = ad %in% c("Int ADNC", "Sev ADNC"),
    LATE_present  = late_stage %in% c("LATE Limbic", "LATE Neocortical"),
    LBD_present   = lbd_stage %in% c("Limbic", "Neocortical"),
    FTLD_present  = ftld_any == 1,
    grossinf_flag  = NACCINF  %in% 1L,
    microinf_flag  = NACCMICR %in% 1L,
    athero_flag    = NACCAVAS %in% c(2L, 3L),
    arteriolo_flag = NACCARTE %in% c(2L, 3L),
    caa_flag       = NACCAMY  %in% c(2L, 3L)
  ) %>%
  filter(!is.na(ptau217), !is.na(ad))

db <- db %>% mutate(mmse_quintile_num = ntile(MMSE, 5))
mmse_ranges <- db %>%
  filter(!is.na(mmse_quintile_num), !is.na(MMSE)) %>%
  group_by(mmse_quintile_num) %>%
  summarize(min_mmse = floor(min(MMSE, na.rm = TRUE)),
            max_mmse = ceiling(max(MMSE, na.rm = TRUE)), .groups = "drop") %>%
  mutate(mmse_label = sprintf("MMSE %d-%d", min_mmse, max_mmse))
db <- db %>%
  left_join(mmse_ranges %>% dplyr::select(mmse_quintile_num, mmse_label),
            by = "mmse_quintile_num") %>%
  mutate(mmse_quintile = factor(mmse_label,
           levels = mmse_ranges$mmse_label[order(mmse_ranges$mmse_quintile_num)]))

age_groups      <- c("<70", "70-80", "80-90", "90+")
sexes           <- c("Male", "Female")
etiology_levels <- c("Clinical AD etiology", "Clinical non-AD etiology")
apoe_status     <- c("APOE e4 -", "APOE e4 +")

colors_etiology <- c("Clinical AD etiology" = "#D7301F", "Clinical non-AD etiology" = "#FDCC8A")
colors_apoe_data <- setNames(unname(colors_apoe), c("APOE e4 -", "APOE e4 +"))
colors_mmse_data <- setNames(colors_mmse_gradient, mmse_ranges$mmse_label)
colors_vascular <- c(
  "No infarcts" = "#E5F5E0", "Gross infarcts" = "#31A354", "Microinfarcts" = "#A1D99B",
  "No microvascular disease" = "#EFEDF5", "Atherosclerosis (mod-sev)" = "#9E9AC8",
  "Arteriolosclerosis (mod-sev)" = "#807DBA", "CAA (mod-sev)" = "#54278F"
)

# Per-stratum AUC (drops strata with <5 of either ADNC class) 
empty_auc <- data.frame(Group = character(), N = integer(), AUC = numeric(),
                        AUC_Lower = numeric(), AUC_Upper = numeric())

stratum_auc <- function(d) {
  d <- d %>% filter(!is.na(ad_binary), !is.na(ptau217))
  if (sum(d$ad_binary == 1) < 5 || sum(d$ad_binary == 0) < 5) return(NULL)
  calc_roc_metrics(d, "ad_binary", "ptau217", direction = "<")
}

auc_by_group <- function(data, var, levels_vec) {
  res <- lapply(levels_vec, function(lv) {
    m <- stratum_auc(data %>% filter(.data[[var]] == lv))
    if (is.null(m)) return(NULL)
    data.frame(Group = lv, N = m$N, AUC = m$AUC,
               AUC_Lower = m$AUC_Lower, AUC_Upper = m$AUC_Upper)
  })
  bind_rows(c(list(empty_auc), res))
}

auc_overall <- {
  m <- stratum_auc(db)
  data.frame(Group = "Full cohort", N = m$N, AUC = m$AUC,
             AUC_Lower = m$AUC_Lower, AUC_Upper = m$AUC_Upper)
}
auc_age      <- auc_by_group(db, "age_group", age_groups)
auc_sex      <- auc_by_group(db, "Sex", sexes)
auc_etiology <- auc_by_group(db, "etiology_group", etiology_levels)
auc_mmse     <- auc_by_group(db, "mmse_quintile", levels(db$mmse_quintile))
auc_apoe     <- auc_by_group(db, "APOE4_status", apoe_status)

copat_strata <- list(
  "No co-pathology" = db %>%
    filter((ADNC_present  & !LATE_present & !LBD_present & !FTLD_present) |
           (!ADNC_present & !LATE_present & !LBD_present & !FTLD_present)) %>%
    mutate(ad_binary = ifelse(ADNC_present, 1, 0)),
  "With LATE" = db %>% filter(LATE_present),
  "With LBD"  = db %>% filter(LBD_present),
  "With FTLD" = db %>% filter(FTLD_present)
)
auc_copat <- bind_rows(c(list(empty_auc), lapply(names(copat_strata), function(sn) {
  m <- stratum_auc(copat_strata[[sn]])
  if (is.null(m)) return(NULL)
  data.frame(Group = dplyr::recode(sn, "No co-pathology" = "None"),
             N = m$N, AUC = m$AUC, AUC_Lower = m$AUC_Lower, AUC_Upper = m$AUC_Upper)
})))

vasc_strata <- list(
  "No infarcts"    = db %>% filter(NACCINF %in% 0L & NACCMICR %in% 0L),
  "Gross infarcts" = db %>% filter(grossinf_flag),
  "Microinfarcts"  = db %>% filter(microinf_flag),
  "No microvascular disease"     = db %>% filter(!athero_flag & !arteriolo_flag & !caa_flag),
  "Atherosclerosis (mod-sev)"    = db %>% filter(athero_flag),
  "Arteriolosclerosis (mod-sev)" = db %>% filter(arteriolo_flag),
  "CAA (mod-sev)"                = db %>% filter(caa_flag)
)
auc_vasc <- bind_rows(c(list(empty_auc), lapply(names(vasc_strata), function(sn) {
  m <- stratum_auc(vasc_strata[[sn]])
  if (is.null(m)) return(NULL)
  data.frame(Group = sn, N = m$N, AUC = m$AUC,
             AUC_Lower = m$AUC_Lower, AUC_Upper = m$AUC_Upper)
})))

# Subject-level cluster bootstrap of delta-AUC
bootstrap_auc_diff <- function(db_full, strata_fns, pred_col = "ptau217",
                               B = 2000, min_n = 10, min_pos = 5, seed = 42) {
  group_nms <- names(strata_fns)
  ids   <- unique(db_full$ID)
  n_ids <- length(ids)
  auc_one <- function(d) {
    d <- d[!is.na(d[[pred_col]]) & !is.na(d$outcome), , drop = FALSE]
    npos <- sum(d$outcome == 1); nneg <- sum(d$outcome == 0)
    if (nrow(d) < min_n || npos < min_pos || nneg < min_pos) return(NA_real_)
    tryCatch(suppressWarnings(as.numeric(pROC::auc(
      pROC::roc(d$outcome, d[[pred_col]], direction = "<", quiet = TRUE)))),
      error = function(e) NA_real_)
  }
  obs_auc <- sapply(group_nms, function(g) auc_one(strata_fns[[g]](db_full)))
  id_to_row <- match(ids, db_full$ID)
  set.seed(seed)
  boot_aucs <- matrix(NA_real_, nrow = B, ncol = length(group_nms),
                      dimnames = list(NULL, group_nms))
  for (b in seq_len(B)) {
    sampled <- sample.int(n_ids, n_ids, replace = TRUE)
    boot_db <- db_full[id_to_row[sampled], , drop = FALSE]
    for (g in group_nms) boot_aucs[b, g] <- auc_one(strata_fns[[g]](boot_db))
  }
  pairs <- combn(group_nms, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    if (is.na(obs_auc[p[1]]) || is.na(obs_auc[p[2]])) return(NULL)
    diffs <- boot_aucs[, p[1]] - boot_aucs[, p[2]]
    diffs <- diffs[!is.na(diffs)]
    if (length(diffs) < B / 2) return(NULL)
    obs_diff <- obs_auc[p[1]] - obs_auc[p[2]]
    ci <- quantile(diffs, c(0.025, 0.975), na.rm = TRUE)
    p_two <- max(2 * min(mean(diffs <= 0), mean(diffs >= 0)), 1 / length(diffs))
    data.frame(Stratum_1 = p[1], AUC_1 = round(unname(obs_auc[p[1]]), 3),
               Stratum_2 = p[2], AUC_2 = round(unname(obs_auc[p[2]]), 3),
               Delta_AUC = round(unname(obs_diff), 3),
               Delta_CI_Lower = round(unname(ci[1]), 3),
               Delta_CI_Upper = round(unname(ci[2]), 3),
               p_raw = p_two, stringsAsFactors = FALSE)
  })
  res <- bind_rows(Filter(Negate(is.null), rows))
  if (nrow(res) == 0) return(NULL)
  res$p_adj_BH <- p.adjust(res$p_raw, method = "BH")
  res
}

make_fns <- function(var, levels_vec) {
  setNames(lapply(levels_vec, function(lv)
    eval(bquote(function(d) d %>%
      filter(.data[[.(var)]] == .(lv), !is.na(ad_binary)) %>%
      mutate(outcome = ad_binary)))), levels_vec)
}

copat_fns <- list(
  "No co-pathology" = function(d) d %>%
    filter((ADNC_present  & !LATE_present & !LBD_present & !FTLD_present) |
           (!ADNC_present & !LATE_present & !LBD_present & !FTLD_present)) %>%
    mutate(outcome = ifelse(ADNC_present, 1L, 0L)),
  "With LATE" = function(d) d %>% filter(LATE_present, !is.na(ad_binary)) %>% mutate(outcome = ad_binary),
  "With LBD"  = function(d) d %>% filter(LBD_present,  !is.na(ad_binary)) %>% mutate(outcome = ad_binary),
  "With FTLD" = function(d) d %>% filter(FTLD_present, !is.na(ad_binary)) %>% mutate(outcome = ad_binary)
)
vasc_fns <- list(
  "No infarcts"    = function(d) d %>% filter(NACCINF %in% 0L & NACCMICR %in% 0L, !is.na(ad_binary)) %>% mutate(outcome = ad_binary),
  "Gross infarcts" = function(d) d %>% filter(grossinf_flag, !is.na(ad_binary)) %>% mutate(outcome = ad_binary),
  "Microinfarcts"  = function(d) d %>% filter(microinf_flag, !is.na(ad_binary)) %>% mutate(outcome = ad_binary)
)
microvasc_fns <- list(
  "No microvascular disease"     = function(d) d %>% filter(!athero_flag & !arteriolo_flag & !caa_flag, !is.na(ad_binary)) %>% mutate(outcome = ad_binary),
  "Atherosclerosis (mod-sev)"    = function(d) d %>% filter(athero_flag,    !is.na(ad_binary)) %>% mutate(outcome = ad_binary),
  "Arteriolosclerosis (mod-sev)" = function(d) d %>% filter(arteriolo_flag, !is.na(ad_binary)) %>% mutate(outcome = ad_binary),
  "CAA (mod-sev)"                = function(d) d %>% filter(caa_flag,        !is.na(ad_binary)) %>% mutate(outcome = ad_binary)
)

B_BOOT <- 2000
boot_tables <- list(
  age       = bootstrap_auc_diff(db, make_fns("age_group", age_groups),      B = B_BOOT, seed = 101),
  sex       = bootstrap_auc_diff(db, make_fns("Sex", sexes),                 B = B_BOOT, seed = 102),
  etiology  = bootstrap_auc_diff(db, make_fns("etiology_group", etiology_levels), B = B_BOOT, seed = 103),
  mmse      = bootstrap_auc_diff(db, make_fns("mmse_quintile", levels(db$mmse_quintile)), B = B_BOOT, seed = 104),
  apoe      = bootstrap_auc_diff(db, make_fns("APOE4_status", apoe_status),  B = B_BOOT, seed = 105),
  copat     = bootstrap_auc_diff(db, copat_fns,     B = B_BOOT, seed = 106),
  vasc      = bootstrap_auc_diff(db, vasc_fns,      B = B_BOOT, seed = 107),
  microvasc = bootstrap_auc_diff(db, microvasc_fns, B = B_BOOT, seed = 108)
)
bind_rows(boot_tables, .id = "Family") %>%
  write_csv(file.path(out_dir, "Figure3_AUC_comparisons.csv"))

# Forest plot assembly 
prep_forest <- function(df, strat) {
  df %>%
    mutate(Stratification = strat,
           label   = sprintf("%.3f (%.3f-%.3f)", AUC, AUC_Lower, AUC_Upper),
           n_label = sprintf("%s (n=%d)", Group, N))
}

make_forest_plot <- function(data, color_vec, title_str, x_lo = 0.495, x_hi = 1.18) {
  data <- data %>% mutate(n_label = factor(n_label, levels = rev(unique(n_label))))
  ggplot(data, aes(x = AUC, y = n_label)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = AUC_Lower, xmax = AUC_Upper), height = 0, linewidth = 0.5, color = "black") +
    geom_point(aes(fill = color_key), shape = 21, size = 2.5, color = "black", stroke = 0.4) +
    geom_text(aes(label = label), hjust = 0.5, vjust = -1.1, size = 2.3, color = "black") +
    scale_fill_manual(values = color_vec, guide = "none") +
    scale_x_continuous(limits = c(x_lo, x_hi), breaks = seq(0.5, 1.0, 0.1),
                       expand = expansion(mult = c(0, 0.03))) +
    facet_grid(Stratification ~ ., scales = "free_y", space = "free_y") +
    theme_manuscript(base_size = 9) +
    theme(panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.major.y = element_blank(), panel.border = element_blank(),
          axis.line.y = element_blank(), strip.text = element_text(face = "bold", size = 8),
          strip.background = element_blank(), panel.spacing = unit(0.5, "lines"), plot.margin = margin(5, 10, 5, 5)) +
    labs(x = "AUC (95% CI)", y = NULL, title = title_str)
}

# Panel A: Demographic
forest_a <- bind_rows(
  prep_forest(auc_overall, "Full cohort"),
  prep_forest(auc_age, "Age at death"),
  prep_forest(auc_sex, "Sex")
) %>%
  mutate(Stratification = factor(Stratification, levels = c("Full cohort", "Age at death", "Sex")),
         color_key = Group) %>%
  arrange(Stratification)
colors_a <- c("Full cohort" = "gray40", colors_age_group, colors_sex)
p_a <- make_forest_plot(forest_a, colors_a, "Figure 3a: Demographic strata")
ggsave(file.path(out_dir, "Figure3a_demographic.pdf"), p_a, width = 5.6, height = 2.6, dpi = 300)

# Panel B: Clinical
forest_b <- bind_rows(
  prep_forest(auc_etiology %>% mutate(Group = dplyr::recode(Group,
                "Clinical AD etiology" = "Clinical AD",
                "Clinical non-AD etiology" = "Clinical non-AD")) %>%
              mutate(color_key = dplyr::recode(Group, "Clinical AD" = "Clinical AD etiology",
                "Clinical non-AD" = "Clinical non-AD etiology")), "Clinical etiology"),
  prep_forest(auc_mmse %>% mutate(color_key = Group), "MMSE quintile")
) %>%
  mutate(Stratification = factor(Stratification, levels = c("Clinical etiology", "MMSE quintile")))
colors_b <- c("Clinical AD etiology" = unname(colors_etiology[["Clinical AD etiology"]]),
              "Clinical non-AD etiology" = unname(colors_etiology[["Clinical non-AD etiology"]]),
              colors_mmse_data)
p_b <- make_forest_plot(forest_b, colors_b, "Figure 3b: Clinical strata")
ggsave(file.path(out_dir, "Figure3b_clinical.pdf"), p_b, width = 5.6, height = 2.6, dpi = 300)

# Panel C: Genetic
forest_c <- prep_forest(auc_apoe %>% mutate(color_key = Group), "APOE e4")
p_c <- make_forest_plot(forest_c, colors_apoe_data, "Figure 3c: Genetic stratum")
ggsave(file.path(out_dir, "Figure3c_genetic.pdf"), p_c, width = 5.6, height = 1.4, dpi = 300)

# Panel D: Co-pathology
forest_d <- bind_rows(
  prep_forest(auc_copat %>% mutate(color_key = Group), "Neurodeg."),
  prep_forest(auc_vasc %>% filter(Group %in% c("No infarcts", "Gross infarcts", "Microinfarcts")) %>%
                mutate(color_key = Group), "Vascular"),
  prep_forest(auc_vasc %>% filter(Group %in% c("No microvascular disease", "Atherosclerosis (mod-sev)",
                "Arteriolosclerosis (mod-sev)", "CAA (mod-sev)")) %>%
                mutate(color_key = Group), "Microvascular")
) %>%
  mutate(Stratification = factor(Stratification,
           levels = c("Neurodeg.", "Vascular", "Microvascular")))
colors_d <- c(colors_copathology_roc[c("No co-pathology", "With LATE", "With LBD", "With FTLD")],
              None = unname(colors_copathology_roc[["No co-pathology"]]),
              colors_vascular)
p_d <- make_forest_plot(forest_d, colors_d, "Figure 3d: Co-pathology strata")
ggsave(file.path(out_dir, "Figure3d_copathology.pdf"), p_d, width = 5.8, height = 3.6, dpi = 300)

cat("Figure 3 written to", out_dir, "\n")
