# Extended Data Figure 1 - Plasma biomarkers across pathology groupings.
#   a  Biomarker concentrations across overall neuropathologic groupings (8 groups)
#   b  p-tau by ADNC and Lewy body disease (2x2, excluding LATE/FTLD/Other)
#   c  p-tau by ADNC and LATE (2x2, excluding LBD/FTLD/Other)

source("setup.R")

db <- readRDS(file.path(data_dir, "demo_cohort.rds")) %>% add_overall_group2()

db <- db %>%
  mutate(
    ADNC_present  = ad %in% c("Int ADNC", "Sev ADNC"),
    LBD_present   = lbd_stage %in% c("Limbic", "Neocortical"),
    LATE_present  = late_stage %in% c("LATE Limbic", "LATE Neocortical"),
    FTLD_present  = FTLD_present,
    Other_present = Other_present,
    LBD_group = case_when(
      !ADNC_present & !LBD_present & !LATE_present & !FTLD_present & !Other_present ~ "Low Path",
      !ADNC_present &  LBD_present & !LATE_present & !FTLD_present & !Other_present ~ "LBD",
       ADNC_present & !LBD_present & !LATE_present & !FTLD_present & !Other_present ~ "ADNC",
       ADNC_present &  LBD_present & !LATE_present & !FTLD_present & !Other_present ~ "ADNC + LBD",
      TRUE ~ NA_character_),
    LBD_group = factor(LBD_group, levels = c("Low Path", "LBD", "ADNC", "ADNC + LBD")),
    LATE_group = case_when(
      !ADNC_present & !LATE_present & !LBD_present & !FTLD_present & !Other_present ~ "Low Path",
      !ADNC_present &  LATE_present & !LBD_present & !FTLD_present & !Other_present ~ "LATE",
       ADNC_present & !LATE_present & !LBD_present & !FTLD_present & !Other_present ~ "ADNC",
       ADNC_present &  LATE_present & !LBD_present & !FTLD_present & !Other_present ~ "ADNC + LATE",
      TRUE ~ NA_character_),
    LATE_group = factor(LATE_group, levels = c("Low Path", "LATE", "ADNC", "ADNC + LATE"))
  )

dunn_positions <- function(plot_data, dunn_results) {
  plot_data %>%
    group_by(biomarker) %>%
    summarize(max_val = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(dunn_results, by = "biomarker") %>%
    group_by(biomarker) %>%
    mutate(bracket_level = row_number(),
           y.position = max_val * (1.05 + (bracket_level - 1) * 0.12)) %>%
    ungroup() %>%
    filter(!is.na(p.adj))
}

biomarker_box <- function(plot_data, color_vec, title_str, annotate = TRUE) {
  kw <- plot_data %>%
    group_by(biomarker) %>%
    kruskal_test(value ~ grp) %>%
    mutate(significant = p < 0.05)
  sig_bio <- kw %>% filter(significant) %>% pull(biomarker)

  dunn_pos <- NULL
  if (annotate && length(sig_bio) > 0) {
    dunn_res <- plot_data %>%
      filter(biomarker %in% sig_bio) %>%
      group_by(biomarker) %>%
      dunn_test(value ~ grp, p.adjust.method = "fdr") %>%
      filter(p.adj < 0.05) %>%
      mutate(p.adj.label = get_stars(p.adj)) %>%
      ungroup()
    if (nrow(dunn_res) > 0) dunn_pos <- dunn_positions(plot_data, dunn_res)
  }

  p <- ggplot(plot_data, aes(x = grp, y = value, fill = grp)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, linewidth = 0.4, width = 0.7) +
    geom_quasirandom(size = 0.5, shape = 21, color = "black", stroke = 0.2,
                     alpha = 0.7, width = 0.2) +
    scale_fill_manual(values = color_vec) +
    facet_wrap(~ biomarker, nrow = 1, scales = "free_y") +
    theme_manuscript(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none",
          panel.spacing = unit(0.3, "lines")) +
    labs(x = NULL, y = "Concentration (pg/mL)", title = title_str)

  if (!is.null(dunn_pos) && nrow(dunn_pos) > 0) {
    p <- p + stat_pvalue_manual(dunn_pos, label = "p.adj.label",
                                tip.length = 0.01, size = 2.5, bracket.size = 0.3,
                                inherit.aes = FALSE)
  }
  p
}

long_biomarkers <- function(data, group_var, markers) {
  data %>%
    filter(!is.na(.data[[group_var]])) %>%
    dplyr::select(ID, grp = all_of(group_var), all_of(markers)) %>%
    pivot_longer(cols = all_of(markers), names_to = "biomarker", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(biomarker = factor(biomarker, levels = markers,
                              labels = unname(biomarker_labels[markers])))
}

# Panel a: all biomarkers across the overall pathology groups
plot_a <- long_biomarkers(db, "Overall_Group", c("ptau217", "ptau181", "GFAP", "NfL"))
p_a <- biomarker_box(plot_a, colors_overall_pathology,
                     "Biomarkers by detailed pathology classification", annotate = FALSE)
ggsave(file.path(out_dir, "ExtData1a_biomarkers_by_overall_group.pdf"),
       p_a, width = 8, height = 2.8, dpi = 300)

# Panel b: biomarkers by ADNC and LBD (2x2)
plot_b <- long_biomarkers(db, "LBD_group", biomarker_vars)
colors_lbd_group <- c("Low Path" = colors_lbd_copathology[["Low Path"]],
                      "LBD" = colors_lbd_copathology[["LBD"]],
                      "ADNC" = colors_lbd_copathology[["ADNC"]],
                      "ADNC + LBD" = colors_lbd_copathology[["ADNC + LBD"]])
p_b <- biomarker_box(plot_b, colors_lbd_group,
                     "Biomarkers by ADNC and Lewy body disease")
ggsave(file.path(out_dir, "ExtData1b_biomarkers_by_lbd.pdf"),
       p_b, width = 7, height = 3, dpi = 300)

# Panel c: biomarkers by ADNC and LATE (2x2)
plot_c <- long_biomarkers(db, "LATE_group", biomarker_vars)
colors_late_group <- c("Low Path" = colors_late_copathology[["Low Path"]],
                       "LATE" = colors_late_copathology[["LATE"]],
                       "ADNC" = colors_late_copathology[["ADNC"]],
                       "ADNC + LATE" = colors_late_copathology[["ADNC + LATE"]])
p_c <- biomarker_box(plot_c, colors_late_group,
                     "Biomarkers by ADNC and LATE")
ggsave(file.path(out_dir, "ExtData1c_biomarkers_by_late.pdf"),
       p_c, width = 7, height = 3, dpi = 300)

cat("Extended Data Figure 1 panels written to", out_dir, "\n")
