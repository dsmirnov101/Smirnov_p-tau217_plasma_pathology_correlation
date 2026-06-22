# Figure 2 - Longitudinal p-tau217 trajectories before death and stability of
#   diagnostic performance across blood-to-autopsy intervals.
#   a  Population LME trajectories of p-tau217 vs years before death by ADNC
#   b  Fixed-effects contrasts vs Low ADNC (delta-method 95% CIs) with divergence years
#   c  Cross-sectional p-tau217 by ADNC and time-to-death interval (Dunn FDR)
#   d  AUC (Int/High vs Low ADNC) within blood-to-death intervals (points + 95% CI)

source("setup.R")

adnc_colors <- c(
  "Low ADNC" = "#2E7D32",
  "Int ADNC" = colors_adnc[["Int ADNC"]],
  "Sev ADNC" = colors_adnc[["Sev ADNC"]]
)
contrast_colors <- c(
  "Sev - Low" = colors_adnc[["Sev ADNC"]],
  "Int - Low" = colors_adnc[["Int ADNC"]]
)

round_up5 <- function(x) ceiling(x / 5) * 5

# Panels a, b: longitudinal LME trajectories + contrasts
db_long <- readRDS(file.path(data_dir, "demo_cohort_long.rds"))

db_traj <- db_long %>%
  filter(!is.na(NPDAGE), !is.na(ad)) %>%
  mutate(
    ad              = factor(ad, levels = c("Low ADNC", "Int ADNC", "Sev ADNC")),
    time_from_death = -as.numeric(time_to_death),
    NACCID          = factor(ID),
    Sex             = factor(Sex)
  ) %>%
  filter(!is.na(time_from_death), !is.na(NPDAGE), !is.na(Sex))

make_pred_grid <- function(model, db) {
  tfd_range <- range(db$time_from_death, na.rm = TRUE)
  tfd_seq   <- seq(tfd_range[1], tfd_range[2], length.out = 200)

  pg <- expand_grid(
    ad              = factor(levels(db$ad), levels = levels(db$ad)),
    time_from_death = tfd_seq,
    NPDAGE          = mean(db$NPDAGE, na.rm = TRUE),
    Sex             = factor(levels(db$Sex)[1], levels = levels(db$Sex))
  )

  mm       <- model.matrix(delete.response(terms(model)), pg)
  pg$pred  <- as.numeric(mm %*% fixef(model))
  pvar     <- diag(mm %*% tcrossprod(as.matrix(vcov(model)), mm))
  pg$se    <- sqrt(pvar)
  pg$lower <- pg$pred - 1.96 * pg$se
  pg$upper <- pg$pred + 1.96 * pg$se
  pg
}

compute_group_contrasts <- function(model, data, biomarker_label) {
  time_grid <- seq(floor(min(data$time_from_death, na.rm = TRUE)), 0, by = 0.5)
  at_vals   <- list(
    time_from_death = time_grid,
    NPDAGE          = mean(data$NPDAGE, na.rm = TRUE),
    Sex             = levels(data$Sex)[1]
  )

  emm <- emmeans(model, ~ ad | time_from_death, at = at_vals,
                 pbkrtest.limit = Inf, lmerTest.limit = Inf)

  ct <- contrast(emm,
    method = list("Sev - Low" = c(-1, 0, 1), "Int - Low" = c(-1, 1, 0))) |>
    confint(level = 0.95) |>
    as.data.frame() |>
    mutate(sig = lower.CL > 0, biomarker = biomarker_label)

  div_sev <- ct |> filter(contrast == "Sev - Low", sig) |>
    slice_min(time_from_death) |> pull(time_from_death)
  div_int <- ct |> filter(contrast == "Int - Low", sig) |>
    slice_min(time_from_death) |> pull(time_from_death)
  if (!length(div_sev)) div_sev <- NA_real_
  if (!length(div_int)) div_int <- NA_real_

  list(contrasts = ct, divergence_sev = div_sev, divergence_int = div_int)
}

plot_trajectory_final <- function(db, col, pred_grid, biomarker_label, ylim = NULL) {
  db2   <- db %>% mutate(.ybd = -time_from_death)
  pg2   <- pred_grid %>% mutate(.ybd = -time_from_death)
  x_max <- round_up5(max(db2$.ybd, na.rm = TRUE))

  p <- ggplot() +
    geom_point(data = db2, aes(.ybd, .data[[col]], color = ad),
               alpha = 0.10, size = 0.55, shape = 16) +
    geom_line(data = db2, aes(.ybd, .data[[col]], color = ad, group = ID),
              alpha = 0.07, linewidth = 0.22) +
    geom_ribbon(data = pg2, aes(.ybd, ymin = lower, ymax = upper, fill = ad), alpha = 0.16) +
    geom_line(data = pg2, aes(.ybd, pred, color = ad), linewidth = 1.35) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray35",
               linewidth = 0.35, alpha = 0.65) +
    scale_x_reverse(limits = c(x_max, 0), breaks = seq(0, x_max, by = 5), expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0.03, 0)) +
    scale_color_manual(values = adnc_colors, name = "ADNC Severity") +
    scale_fill_manual(values = adnc_colors, name = "ADNC Severity") +
    theme_manuscript(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.key.size  = unit(0.35, "cm"),
      legend.text      = element_text(size = 7.5),
      legend.title     = element_text(size = 8, face = "bold"),
      legend.margin    = margin(t = -2, b = 2)
    ) +
    labs(x = "Years Before Death", y = sprintf("%s (pg/mL)", biomarker_label),
         title = "Figure 2a: p-tau217 trajectories by ADNC")

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

plot_contrast_final <- function(ct_result, y_label, ylim = NULL) {
  ct_df       <- ct_result$contrasts
  div_sev     <- ct_result$divergence_sev
  div_int     <- ct_result$divergence_int
  ct_df2      <- ct_df %>% mutate(.ybd = -time_from_death)
  div_sev_pos <- if (!is.na(div_sev)) -div_sev else NA_real_
  div_int_pos <- if (!is.na(div_int)) -div_int else NA_real_
  x_max       <- round_up5(max(ct_df2$.ybd, na.rm = TRUE))

  p <- ggplot(ct_df2, aes(x = .ybd, color = contrast, fill = contrast))

  if (!is.na(div_sev_pos))
    p <- p + annotate("rect", xmin = 0, xmax = div_sev_pos, ymin = -Inf, ymax = Inf,
                      fill = contrast_colors["Sev - Low"], alpha = 0.055)
  if (!is.na(div_int_pos))
    p <- p + annotate("rect", xmin = 0, xmax = div_int_pos, ymin = -Inf, ymax = Inf,
                      fill = contrast_colors["Int - Low"], alpha = 0.055)

  p <- p +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "gray20") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55",
               linewidth = 0.30, alpha = 0.65) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.22, color = NA) +
    geom_line(aes(y = estimate), linewidth = 1.15) +
    scale_x_reverse(limits = c(x_max, 0), breaks = seq(0, x_max, by = 5), expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0.05, 0)) +
    scale_color_manual(values = contrast_colors, name = NULL) +
    scale_fill_manual(values = contrast_colors, name = NULL) +
    theme_manuscript(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.key.size  = unit(0.35, "cm"),
      legend.text      = element_text(size = 7.5),
      legend.margin    = margin(t = -2, b = 2)
    ) +
    labs(x = "Years Before Death", y = sprintf("Δ %s (pg/mL)", y_label),
         title = "Figure 2b: ADNC vs Low ADNC contrasts")

  if (!is.na(div_sev_pos))
    p <- p +
      geom_vline(xintercept = div_sev_pos, linetype = "dotted", linewidth = 0.65,
                 color = contrast_colors["Sev - Low"]) +
      annotate("text", x = div_sev_pos, y = Inf,
               label = sprintf("Sev−Low:\n%.0f yr", div_sev_pos),
               color = contrast_colors["Sev - Low"], vjust = 1.4, hjust = 1.1,
               size = 2.5, lineheight = 0.9)
  if (!is.na(div_int_pos))
    p <- p +
      geom_vline(xintercept = div_int_pos, linetype = "dotted", linewidth = 0.65,
                 color = contrast_colors["Int - Low"]) +
      annotate("text", x = div_int_pos, y = Inf,
               label = sprintf("Int−Low:\n%.0f yr", div_int_pos),
               color = contrast_colors["Int - Low"], vjust = 3.8, hjust = 1.1,
               size = 2.5, lineheight = 0.9)

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

db_ptau217_long <- db_traj %>% filter(!is.na(ptau217)) %>% droplevels()

model_ptau217 <- lmer(
  ptau217 ~ ad * ns(time_from_death, df = 2) + NPDAGE + Sex + (1 | NACCID),
  data    = db_ptau217_long,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)

pg_ptau217 <- make_pred_grid(model_ptau217, db_ptau217_long)
ct_ptau217 <- compute_group_contrasts(model_ptau217, db_ptau217_long, "p-tau217")

p_a <- plot_trajectory_final(db_ptau217_long, "ptau217", pg_ptau217, "p-tau217", ylim = c(0, 3))
ggsave(file.path(out_dir, "Figure2a_ptau217_trajectory.pdf"), p_a, width = 4, height = 3.6, dpi = 300)

p_b <- plot_contrast_final(ct_ptau217, "p-tau217")
ggsave(file.path(out_dir, "Figure2b_ptau217_contrast.pdf"), p_b, width = 4, height = 3.6, dpi = 300)

# Panels c, d: cross-sectional stratification by time-to-death interval
db <- readRDS(file.path(data_dir, "demo_cohort.rds")) %>%
  filter(!is.na(ptau217), !is.na(ad)) %>%
  add_adnc_binary() %>%
  mutate(
    ad         = factor(ad, levels = c("Low ADNC", "Int ADNC", "Sev ADNC")),
    time_group = cut(as.numeric(time_to_death),
                     breaks = c(-Inf, 5, 10, 15, Inf),
                     labels = c("0-5 years", "5-10 years", "10-15 years", "15+ years"),
                     right  = FALSE)
  )

time_groups <- c("0-5 years", "5-10 years", "10-15 years", "15+ years")

# Panel c: boxplots with Dunn post-hoc (FDR)
plot_data <- db %>% filter(!is.na(time_group))

dunn_groups <- plot_data %>%
  group_by(time_group) %>%
  filter(n_distinct(ad) > 1) %>%
  ungroup()

dunn_results <- dunn_groups %>%
  group_by(time_group) %>%
  dunn_test(ptau217 ~ ad, p.adjust.method = "fdr") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.label = get_stars(p.adj)) %>%
  ungroup()

dunn_with_positions <- plot_data %>%
  group_by(time_group) %>%
  summarize(max_val = max(ptau217, na.rm = TRUE), .groups = "drop") %>%
  left_join(dunn_results, by = "time_group") %>%
  group_by(time_group) %>%
  mutate(bracket_level = row_number(),
         y.position = max_val * (1.05 + (bracket_level - 1) * 0.12)) %>%
  ungroup() %>%
  filter(!is.na(p.adj))

p_c <- ggplot(plot_data, aes(ad, ptau217)) +
  geom_boxplot(aes(fill = ad), outlier.shape = NA, alpha = 0.75, linewidth = 0.4, width = 0.7) +
  geom_quasirandom(aes(fill = ad), size = 0.5, shape = 21, color = "black", stroke = 0.2, alpha = 0.7, width = 0.2) +
  scale_fill_manual(values = colors_adnc) +
  facet_wrap(~ time_group, nrow = 1) +
  theme_manuscript(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none", panel.spacing = unit(0.2, "lines")) +
  labs(x = NULL, y = "p-tau217 (pg/mL)",
       title = "Figure 2c: p-tau217 by ADNC and time to death")

if (nrow(dunn_with_positions) > 0)
  p_c <- p_c + stat_pvalue_manual(dunn_with_positions, label = "p.adj.label",
                                  tip.length = 0.01, size = 2.5, bracket.size = 0.3)
ggsave(file.path(out_dir, "Figure2c_ptau217_by_time.pdf"), p_c, width = 6, height = 3, dpi = 300)

# Panel d: AUC within intervals (95% CI)
auc_time <- bind_rows(lapply(time_groups, function(tg) {
  m <- calc_roc_metrics(db %>% filter(time_group == tg), "ad_binary", "ptau217", direction = "<")
  if (is.null(m)) return(NULL)
  m$Time_Group <- tg
  m
}))

auc_time <- auc_time %>%
  mutate(Group = factor(Time_Group, levels = rev(time_groups)),
         label = sprintf("%.3f (%.3f-%.3f)", AUC, AUC_Lower, AUC_Upper))

p_d <- ggplot(auc_time, aes(x = AUC, y = Group)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = AUC_Lower, xmax = AUC_Upper), height = 0.2, linewidth = 0.5, color = "gray30") +
  geom_point(aes(fill = Time_Group), shape = 21, size = 3, color = "black", stroke = 0.4) +
  geom_text(aes(label = label), hjust = -0.15, size = 2.5, color = "gray20") +
  scale_fill_manual(values = colors_time_to_death_data, guide = "none") +
  scale_x_continuous(limits = c(0.4, 1.35), breaks = seq(0.5, 1.0, 0.1)) +
  theme_manuscript(base_size = 9) +
  theme(panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.major.y = element_blank()) +
  labs(x = "AUC (95% CI)", y = NULL, title = "Figure 2d: p-tau217 AUC by time to death")
ggsave(file.path(out_dir, "Figure2d_auc_by_time.pdf"), p_d, width = 5.6, height = 2.5, dpi = 300)

cat("Figure 2 panels written to", out_dir, "\n")
