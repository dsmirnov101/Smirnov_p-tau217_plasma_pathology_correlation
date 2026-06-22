# Figure 4 - Biomarker-defined cognitive decline and clinical progression.
#   a  Annual rate of MMSE change (slope) by biomarker status, from lmer
#      biomarker x time interaction models
#   b  Modeled MMSE trajectories by p-tau217 status, stratified by co-pathology
#   c  Modeled MMSE trajectories by p-tau217 status, stratified by baseline dementia
#   d  Kaplan-Meier dementia-free survival among non-demented at baseline, by biomarker
# Binary biomarker status from autopsy-derived Youden cutpoints (Fig 1c); NfL median split.

source("setup.R")

db_long <- readRDS(file.path(data_dir, "demo_cohort_long.rds"))
db <- readRDS(file.path(data_dir, "demo_cohort.rds")) %>% add_adnc_binary()

TIME_AFTER_BASELINE <- 10
biomarker_vars_raw <- c("ptau217", "ptau181", "GFAP", "NfL")
biomarker_status_labels <- c("p-tau217", "p-tau181", "GFAP", "NfL")
biomarker_status_factor_vars <- c("ptau217_status", "ptau181_status", "GFAP_status", "NfL_status")
biomarker_pos_vars <- c("ptau217_pos", "ptau181_pos", "GFAP_pos", "NfL_pos")

# Youden cutpoints from cross-sectional ADNC ROC (NfL = median split)
youden_cutpoints <- setNames(rep(NA_real_, length(biomarker_vars_raw)), biomarker_vars_raw)
for (bio in biomarker_vars_raw) {
  m <- calc_roc_metrics(db, "ad_binary", bio, direction = "<")
  if (!is.null(m)) youden_cutpoints[bio] <- m$Cutoff
}
nfl_median <- median(db$NfL, na.rm = TRUE)
youden_cutpoints["NfL"] <- nfl_median

apply_status <- function(d) {
  d %>% mutate(
    ptau217_pos = as.numeric(ptau217_src >= youden_cutpoints["ptau217"]),
    ptau181_pos = as.numeric(ptau181_src >= youden_cutpoints["ptau181"]),
    GFAP_pos    = as.numeric(GFAP_src    >= youden_cutpoints["GFAP"]),
    NfL_pos     = as.numeric(NfL_src     >= youden_cutpoints["NfL"]),
    ptau217_status = factor(ifelse(ptau217_pos == 1, "p-tau217 +", "p-tau217 -"),
                            levels = c("p-tau217 -", "p-tau217 +")),
    ptau181_status = factor(ifelse(ptau181_pos == 1, "p-tau181 +", "p-tau181 -"),
                            levels = c("p-tau181 -", "p-tau181 +")),
    GFAP_status = factor(ifelse(GFAP_pos == 1, "GFAP +", "GFAP -"),
                         levels = c("GFAP -", "GFAP +")),
    NfL_status = factor(ifelse(NfL_pos == 1, "NfL +", "NfL -"),
                        levels = c("NfL -", "NfL +"))
  )
}

# Baseline
baseline_biomarkers <- db_long %>%
  filter(!is.na(ptau217)) %>%
  group_by(ID) %>%
  arrange(VISEQ) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(ID, first_viseq = VISEQ, baseline_age = AgeAtVisit,
                ptau217_src = ptau217, ptau181_src = ptau181,
                GFAP_src = GFAP, NfL_src = NfL)

#Longitudinal trajectory data
db_traj <- db_long %>%
  inner_join(baseline_biomarkers, by = "ID") %>%
  add_overall_group2() %>%
  mutate(
    time_from_baseline = AgeAtVisit - baseline_age,
    age_centered = baseline_age - mean(baseline_age, na.rm = TRUE),
    educ_centered = Educ - mean(Educ, na.rm = TRUE),
    sex_factor = factor(Sex),
    copathology_group = case_when(
      Overall_Group2 %in% c("Other", "Mixed ADNC") ~ "Co-pathology Present",
      Overall_Group2 %in% c("Low Path", "ADNC") ~ "No Co-pathology",
      TRUE ~ NA_character_),
    copathology_group = factor(copathology_group,
                               levels = c("No Co-pathology", "Co-pathology Present"))
  ) %>%
  filter(time_from_baseline >= 0, time_from_baseline <= TIME_AFTER_BASELINE) %>%
  apply_status() %>%
  filter(!is.na(age_centered), !is.na(sex_factor), !is.na(educ_centered))

baseline_dx_for_dem <- db_traj %>%
  group_by(ID) %>%
  arrange(VISEQ) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(ID, baseline_dx_dem = dx2) %>%
  mutate(dementia_status = case_when(
    baseline_dx_dem %in% c("Cognitively Normal", "MCI") ~ "Non-Dementia",
    baseline_dx_dem == "Dementia" ~ "Dementia",
    TRUE ~ NA_character_))

db_traj_dem <- db_traj %>%
  left_join(baseline_dx_for_dem, by = "ID") %>%
  filter(!is.na(dementia_status)) %>%
  mutate(dementia_status = factor(dementia_status, levels = c("Non-Dementia", "Dementia")))

# lmer
fit_traj <- function(data, status_var, min_subj = 20) {
  temp <- data %>% filter(!is.na(.data[[status_var]]), !is.na(MMSE))
  if (length(unique(temp$ID)) < min_subj) return(NULL)
  n_subj <- length(unique(temp$ID))
  f_rs <- as.formula(paste0("MMSE ~ ", status_var, " * time_from_baseline + ",
    "time_from_baseline * (age_centered + sex_factor + educ_centered) + (1 + time_from_baseline | ID)"))
  f_ri <- as.formula(paste0("MMSE ~ ", status_var, " * time_from_baseline + ",
    "time_from_baseline * (age_centered + sex_factor + educ_centered) + (1 | ID)"))
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  model <- tryCatch(lmer(f_rs, data = temp, REML = TRUE, control = ctrl),
    error = function(e) tryCatch(lmer(f_ri, data = temp, REML = TRUE, control = ctrl),
      error = function(e2) NULL))
  if (is.null(model)) return(NULL)

  coef_tidy <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  it <- coef_tidy %>% filter(grepl(paste0(status_var, ".*time_from_baseline|time_from_baseline.*", status_var), term))
  slope_diff_beta    <- if (nrow(it) > 0) it$estimate[1]  else NA
  slope_diff_p       <- if (nrow(it) > 0) it$p.value[1]   else NA
  slope_diff_ci_low  <- if (nrow(it) > 0) it$conf.low[1]  else NA
  slope_diff_ci_high <- if (nrow(it) > 0) it$conf.high[1] else NA

  time_seq <- seq(0, TIME_AFTER_BASELINE, by = 0.5)
  mode_sex <- names(sort(table(temp$sex_factor), decreasing = TRUE))[1]
  pred_grid <- expand_grid(
    !!status_var := levels(temp[[status_var]]),
    time_from_baseline = time_seq,
    age_centered = mean(temp$age_centered, na.rm = TRUE),
    sex_factor = factor(mode_sex, levels = levels(temp$sex_factor)),
    educ_centered = mean(temp$educ_centered, na.rm = TRUE))
  pred_grid[[status_var]] <- factor(pred_grid[[status_var]], levels = levels(temp[[status_var]]))
  pred_grid$predicted_MMSE <- predict(model, newdata = pred_grid, re.form = NA)
  mm <- model.matrix(delete.response(terms(model)), pred_grid)
  mm <- mm[, colnames(vcov(model)), drop = FALSE]
  pred_grid$se <- sqrt(diag(mm %*% tcrossprod(vcov(model), mm)))
  pred_grid$ci_low  <- pred_grid$predicted_MMSE - 1.96 * pred_grid$se
  pred_grid$ci_high <- pred_grid$predicted_MMSE + 1.96 * pred_grid$se
  pred_grid$Group <- pred_grid[[status_var]]

  raw <- temp %>% dplyr::select(ID, time_from_baseline, MMSE, !!status_var) %>%
    rename(Group = !!status_var)

  emm <- tryCatch(emtrends(model, as.formula(paste0("~ ", status_var)),
    var = "time_from_baseline", lmer.df = "satterthwaite"), error = function(e) NULL)
  group_slopes <- NULL
  if (!is.null(emm)) {
    es <- as.data.frame(summary(emm))
    group_slopes <- data.frame(Group = as.character(es[[1]]),
      slope = es$time_from_baseline.trend, se = es$SE,
      ci_low = es$lower.CL, ci_high = es$upper.CL, stringsAsFactors = FALSE)
  }
  list(predictions = pred_grid, raw_data = raw, group_slopes = group_slopes, n_subj = n_subj,
       slope_diff_beta = slope_diff_beta, slope_diff_p = slope_diff_p,
       slope_diff_ci_low = slope_diff_ci_low, slope_diff_ci_high = slope_diff_ci_high)
}

# Slopes across strata for forest plot 
group_slopes_summary <- data.frame()
slope_summary <- data.frame()

add_slopes <- function(res, bio_label, stratum, condition) {
  if (is.null(res)) return(invisible(NULL))
  slope_summary <<- rbind(slope_summary, data.frame(
    Biomarker = bio_label, Stratum = stratum, Condition = condition,
    beta = res$slope_diff_beta, ci_low = res$slope_diff_ci_low,
    ci_high = res$slope_diff_ci_high, p_value = res$slope_diff_p,
    n_subj = res$n_subj, stringsAsFactors = FALSE))
  if (!is.null(res$group_slopes)) {
    group_slopes_summary <<- rbind(group_slopes_summary, data.frame(
      Biomarker = bio_label, Stratum = stratum, Condition = condition,
      Group = res$group_slopes$Group, slope = res$group_slopes$slope,
      se = res$group_slopes$se, ci_low = res$group_slopes$ci_low,
      ci_high = res$group_slopes$ci_high, n_subj = res$n_subj, stringsAsFactors = FALSE))
  }
}

traj_results <- list()
for (i in seq_along(biomarker_status_factor_vars)) {
  res <- fit_traj(db_traj, biomarker_status_factor_vars[i], min_subj = 30)
  traj_results[[biomarker_status_labels[i]]] <- res
  add_slopes(res, biomarker_status_labels[i], "Full Cohort", "Full Cohort")
}

for (i in seq_along(biomarker_status_factor_vars)) {
  for (lv in c("No Co-pathology", "Co-pathology Present")) {
    res <- fit_traj(db_traj %>% filter(copathology_group == lv),
                    biomarker_status_factor_vars[i])
    add_slopes(res, biomarker_status_labels[i], "Co-pathology", lv)
  }
}

dem_traj_results <- list()
for (i in seq_along(biomarker_status_factor_vars)) {
  for (lv in c("Non-Dementia", "Dementia")) {
    res <- fit_traj(db_traj_dem %>% filter(dementia_status == lv),
                    biomarker_status_factor_vars[i])
    add_slopes(res, biomarker_status_labels[i], "Dementia Status", lv)
    if (biomarker_status_labels[i] == "p-tau217") dem_traj_results[[lv]] <- res
  }
}


# PANEL A: forest plot of per-group annual MMSE slopes

condition_order <- c("Full Cohort", "No Co-pathology", "Co-pathology Present",
                     "Non-Dementia", "Dementia")
biomarker_colors <- c("p-tau217" = "#B2182B", "p-tau181" = "#E66101",
                      "GFAP" = "#1B7837", "NfL" = "#762A83")

if (nrow(group_slopes_summary) > 0) {
  slope_summary <- slope_summary %>%
    mutate(sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                           p_value < 0.05 ~ "*", TRUE ~ ""),
           Condition = factor(Condition, levels = condition_order),
           Biomarker = factor(Biomarker, levels = biomarker_status_labels))

  forest_data <- group_slopes_summary %>%
    mutate(pos_neg = ifelse(grepl("\\+$", Group), "Positive", "Negative"),
           pos_neg = factor(pos_neg, levels = c("Negative", "Positive")),
           Condition = factor(Condition, levels = condition_order),
           Biomarker = factor(Biomarker, levels = biomarker_status_labels))

  sig_df <- slope_summary %>% group_by(Biomarker, Condition) %>% slice(1) %>% ungroup()
  star_df <- forest_data %>%
    group_by(Biomarker, Condition) %>%
    summarize(y = max(ci_high, na.rm = TRUE), .groups = "drop") %>%
    left_join(sig_df %>% dplyr::select(Biomarker, Condition, sig), by = c("Biomarker", "Condition"))
  dodge_width <- 0.55

  p_a <- ggplot(forest_data,
                aes(x = Biomarker, y = slope, color = Biomarker, shape = pos_neg, group = pos_neg)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.7,
                  position = position_dodge(width = dodge_width)) +
    geom_point(size = 3, position = position_dodge(width = dodge_width)) +
    geom_text(data = star_df, aes(x = Biomarker, y = y, label = sig), inherit.aes = FALSE,
              vjust = -0.4, size = 4, color = "gray25") +
    facet_wrap(~ Condition, nrow = 1) +
    scale_color_manual(values = biomarker_colors, name = "Biomarker") +
    scale_shape_manual(values = c("Negative" = 1, "Positive" = 16),
                       labels = c("Negative" = "Biomarker -", "Positive" = "Biomarker +"),
                       name = "Status") +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3)),
           shape = guide_legend(override.aes = list(size = 3, color = "black"))) +
    labs(x = NULL, y = "Estimated MMSE change (points/year)",
         title = "Annual Rate of MMSE Decline by Biomarker Status") +
    theme_manuscript(base_size = 9) +
    theme(legend.position = "bottom", strip.text = element_text(size = 9, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3))

  ggsave(file.path(out_dir, "Figure4a_slope_forest.pdf"), p_a, width = 11, height = 4.2, dpi = 300)
}

# PANELS B / C: p-tau217 MMSE trajectories stratified by co-pathology / dementia

plot_strat_trajectories <- function(strat_data, status_var, strat_var, bio_label, levels_vec, title_str) {
  preds <- list(); raws <- list()
  annot <- data.frame()
  for (lv in levels_vec) {
    res <- fit_traj(strat_data %>% filter(.data[[strat_var]] == lv), status_var)
    if (is.null(res)) next
    pg <- res$predictions; pg[[strat_var]] <- lv; preds[[lv]] <- pg
    rw <- res$raw_data; rw[[strat_var]] <- lv; raws[[lv]] <- rw
    annot <- rbind(annot, data.frame(
      lvl = lv,
      label = sprintf("Slope diff: b=%.3f\np=%s, N=%d",
                      res$slope_diff_beta, format_pvalue(res$slope_diff_p), res$n_subj),
      stringsAsFactors = FALSE))
  }
  if (length(preds) == 0) return(NULL)
  all_preds <- bind_rows(preds); all_raw <- bind_rows(raws)
  all_preds[[strat_var]] <- factor(all_preds[[strat_var]], levels = levels_vec)
  all_raw[[strat_var]]   <- factor(all_raw[[strat_var]],   levels = levels_vec)
  neg_label <- paste0(bio_label, " -"); pos_label <- paste0(bio_label, " +")
  fill_colors <- setNames(c("#2166AC", "#B2182B"), c(neg_label, pos_label))

  p <- ggplot() +
    geom_line(data = all_raw, aes(x = time_from_baseline, y = MMSE, group = ID, color = Group),
              alpha = 0.1, linewidth = 0.3) +
    geom_ribbon(data = all_preds, aes(x = time_from_baseline, ymin = ci_low, ymax = ci_high, fill = Group),
                alpha = 0.2) +
    geom_line(data = all_preds, aes(x = time_from_baseline, y = predicted_MMSE, color = Group),
              linewidth = 1.3) +
    facet_wrap(as.formula(paste0("~ ", strat_var)), ncol = 2) +
    scale_color_manual(values = fill_colors) +
    scale_fill_manual(values = fill_colors) +
    scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
    coord_cartesian(ylim = c(0, 30)) +
    theme_manuscript(base_size = 9) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(size = 12, face = "bold")) +
    labs(x = "Time from Baseline (years)", y = "MMSE Score", title = title_str)

  if (nrow(annot) > 0) {
    annot[[strat_var]] <- factor(annot$lvl, levels = levels_vec)
    p <- p + geom_text(data = annot, aes(x = 0.5, y = 3, label = label),
                       hjust = 0, vjust = 0, size = 3, color = "gray30")
  }
  p
}

p_b <- plot_strat_trajectories(db_traj, "ptau217_status", "copathology_group", "p-tau217",
  c("No Co-pathology", "Co-pathology Present"),
  "p-tau217: MMSE Trajectories by Co-pathology Status")
if (!is.null(p_b)) ggsave(file.path(out_dir, "Figure4b_ptau217_traj_copathology.pdf"),
                          p_b, width = 8, height = 4, dpi = 300)

p_c <- plot_strat_trajectories(db_traj_dem, "ptau217_status", "dementia_status", "p-tau217",
  c("Non-Dementia", "Dementia"),
  "p-tau217: MMSE Trajectories by Baseline Dementia Status")
if (!is.null(p_c)) ggsave(file.path(out_dir, "Figure4c_ptau217_traj_dementia.pdf"),
                          p_c, width = 8, height = 4, dpi = 300)


# PANEL D: Kaplan-Meier dementia-free survival

baseline_dx_data <- db_long %>%
  filter(!is.na(ptau217)) %>%
  group_by(ID) %>% arrange(VISEQ) %>% slice(1) %>% ungroup() %>%
  dplyr::select(ID, baseline_dx = dx2, baseline_age = AgeAtVisit)

conversion_summary <- db_long %>%
  filter(!is.na(dx2)) %>%
  inner_join(baseline_dx_data, by = "ID") %>%
  mutate(reached_dementia = as.numeric(dx2 == "Dementia")) %>%
  group_by(ID) %>% arrange(VISEQ) %>%
  summarize(baseline_dx = first(baseline_dx),
            last_age = last(AgeAtVisit),
            converted_to_dementia = as.numeric(any(reached_dementia == 1, na.rm = TRUE)),
            first_dementia_age = ifelse(any(reached_dementia == 1, na.rm = TRUE),
                                        min(AgeAtVisit[reached_dementia == 1]), NA_real_),
            .groups = "drop") %>%
  left_join(baseline_dx_data %>% dplyr::select(ID, baseline_age), by = "ID") %>%
  mutate(followup_years = last_age - baseline_age,
         converted = coalesce(converted_to_dementia, 0),
         time_to_event = ifelse(converted == 1, first_dementia_age - baseline_age, followup_years)) %>%
  inner_join(baseline_biomarkers %>%
               dplyr::select(ID, ptau217_src, ptau181_src, GFAP_src, NfL_src), by = "ID") %>%
  filter(baseline_dx %in% c("Cognitively Normal", "MCI"),
         !is.na(time_to_event), time_to_event > 0) %>%
  apply_status()

km_plots <- list()
for (i in seq_along(biomarker_status_factor_vars)) {
  status_var <- biomarker_status_factor_vars[i]
  bio_label <- biomarker_status_labels[i]
  km_data <- conversion_summary %>% filter(!is.na(.data[[status_var]]))
  if (nrow(km_data) < 20 || sum(km_data$converted) < 1) next
  if (length(unique(km_data[[status_var]])) < 2) next

  surv_formula <- as.formula(paste0("Surv(time_to_event, converted) ~ ", status_var))
  fit <- tryCatch(survfit(surv_formula, data = km_data), error = function(e) NULL)
  if (is.null(fit)) next
  fit$call$formula <- surv_formula
  pval <- tryCatch({
    lr <- survdiff(surv_formula, data = km_data)
    1 - pchisq(lr$chisq, length(lr$n) - 1)
  }, error = function(e) NA_real_)

  neg_label <- paste0(bio_label, " -"); pos_label <- paste0(bio_label, " +")
  gg <- tryCatch(ggsurvplot(fit, data = km_data,
    palette = c("#2166AC", "#B2182B"),
    conf.int = TRUE, conf.int.alpha = 0.12,
    legend.labs = c(neg_label, pos_label), legend.title = "",
    ylim = c(0, 1), ylab = "Dementia-Free Probability",
    xlab = "Time from Baseline (years)",
    pval = ifelse(is.na(pval), FALSE, paste0("Log-rank p = ", format_pvalue(pval))),
    pval.size = 3, ggtheme = theme_manuscript(base_size = 9)),
    error = function(e) NULL)
  if (is.null(gg)) next
  km_plots[[bio_label]] <- gg$plot + labs(title = bio_label) +
    theme(legend.position = "bottom")
}

if (length(km_plots) > 0) {
  p_d <- wrap_plots(km_plots, nrow = 1) +
    plot_annotation(title = "Dementia-Free Survival by Biomarker Status",
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
  ggsave(file.path(out_dir, "Figure4d_KM_conversion.pdf"), p_d,
         width = 3.5 * length(km_plots), height = 3.5, dpi = 300)
}

cat("Figure 4 panels written to", out_dir, "\n")
