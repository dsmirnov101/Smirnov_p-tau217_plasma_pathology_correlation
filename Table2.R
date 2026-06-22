# Generates the synthetic demonstration cohort used by all figure and table
# scripts. Values are random draws with a fixed seed and contain no real
# participant information. Produces ../demo_data/demo_cohort.rds (one row per
# participant) and ../demo_data/demo_cohort_long.rds (multiple visits).

suppressPackageStartupMessages(library(tidyverse))
set.seed(42)

out_dir <- "../demo_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
N <- 200

ad_levels <- c("Low ADNC", "Int ADNC", "Sev ADNC")
ad  <- factor(sample(ad_levels, N, replace = TRUE, prob = c(0.22, 0.20, 0.58)), levels = ad_levels)
sev <- as.integer(ad)

braak_by_ad <- c("Low ADNC" = "I-II", "Int ADNC" = "III-IV", "Sev ADNC" = "V-VI")
Braak_grouped <- factor(ifelse(ad == "Low ADNC" & runif(N) < 0.4, "0", braak_by_ad[as.character(ad)]),
                        levels = c("0", "I-II", "III-IV", "V-VI"))
neur_levels <- c("None", "Sparse", "Moderate", "Frequent")
NEUR <- factor(neur_levels[pmin(4L, pmax(1L, sev + sample(-1:1, N, replace = TRUE)))], levels = neur_levels)

# Plasma biomarkers: lognormal, group medians near the manuscript's Youden  cutpoints
lnorm <- function(m, s) rlnorm(N, log(m[sev]), s)
ptau217 <- lnorm(c(0.30, 0.62, 1.25), 0.55)
ptau181 <- lnorm(c(28, 44, 60), 0.45)
GFAP    <- lnorm(c(175, 300, 405), 0.45)
NfL     <- rlnorm(N, log(48), 0.5)
ptau181[sample(N, 4)] <- NA; GFAP[sample(N, 2)] <- NA; NfL[sample(N, 2)] <- NA
log_ptau217 <- log(ptau217); log_ptau181 <- log(ptau181); log_GFAP <- log(GFAP); log_NfL <- log(NfL)

AgeAtVisit    <- pmin(99, pmax(55, round(rnorm(N, 78, 9), 1)))
time_to_death <- round(pmin(25, rgamma(N, shape = 2.2, scale = 4.2)), 2)  # wide spread; ~15% beyond 15 y
NPDAGE        <- pmin(104, round(AgeAtVisit + time_to_death))
AGEOFONSET    <- pmax(45, round(AgeAtVisit - rgamma(N, 2, 0.4)))
Sex   <- sample(c(1, 2), N, replace = TRUE)
Race  <- sample(c(1, 2, 3), N, replace = TRUE, prob = c(0.85, 0.10, 0.05))
Educ  <- pmin(22, pmax(6, round(rnorm(N, 15, 3))))
apoe4 <- pmin(2, sample(0:2, N, replace = TRUE, prob = c(0.55, 0.35, 0.10)) +
                (sev == 3) * sample(0:1, N, replace = TRUE, prob = c(0.7, 0.3)))
APOE  <- c(33, 34, 44)[apoe4 + 1]
VISEQ <- sample(1:6, N, replace = TRUE)

dx2_levels <- c("Cognitively Normal", "MCI", "Dementia", "Excluded/Unknown")
dx2 <- factor(dx2_levels[pmin(3L, pmax(1L, sev + sample(-1:1, N, replace = TRUE)))], levels = dx2_levels)
dx_levels <- c("Control", "MCI", "AD Dementia", "LBD Dementia", "FTD Dementia", "Vascular Dementia",
               "Parkinsonian Dementia", "Other Dementia", "Excluded/Unknown", "Other/Unknown")
dx <- factor(case_when(
  dx2 == "Cognitively Normal" ~ "Control", dx2 == "MCI" ~ "MCI",
  TRUE ~ sample(c("AD Dementia", "LBD Dementia", "FTD Dementia", "Vascular Dementia"),
                N, replace = TRUE, prob = c(0.6, 0.15, 0.15, 0.10))), levels = dx_levels)
MMSE <- round(pmin(30, pmax(0, case_when(
  dx2 == "Cognitively Normal" ~ rnorm(N, 28, 1.5), dx2 == "MCI" ~ rnorm(N, 24, 2.5), TRUE ~ rnorm(N, 16, 5)))))
CDRSum <- round(pmin(18, pmax(0, case_when(
  dx2 == "Cognitively Normal" ~ rnorm(N, 0.3, 0.4), dx2 == "MCI" ~ rnorm(N, 2, 1), TRUE ~ rnorm(N, 9, 4)))), 1)
PersCare <- as.character(pmin(3, pmax(0, round(CDRSum / 5))))

lbd_stage  <- factor(sample(c("Absent", "Brainstem", "Limbic", "Neocortical"), N, TRUE, c(0.70, 0.10, 0.10, 0.10)),
                     levels = c("Absent", "Brainstem", "Limbic", "Neocortical"))
late_stage <- factor(sample(c("Absent", "LATE Amygdala", "LATE Limbic", "LATE Neocortical"), N, TRUE, c(0.68, 0.12, 0.12, 0.08)),
                     levels = c("Absent", "LATE Amygdala", "LATE Limbic", "LATE Neocortical"))

ftld_any   <- rbinom(N, 1, 0.18)
ftld_tau   <- ftld_any * rbinom(N, 1, 0.55)
ftld_tdp   <- ftld_any * (1 - ftld_tau) * rbinom(N, 1, 0.7)
ftld_fus   <- ftld_any * (1 - ftld_tau) * (1 - ftld_tdp) * rbinom(N, 1, 0.3)
ftld_other <- ftld_any * (1 - ftld_tau) * (1 - ftld_tdp) * (1 - ftld_fus)
ftld_mnd   <- ftld_tdp * rbinom(N, 1, 0.3)
Other_present <- rbinom(N, 1, 0.12) == 1
FTLD_present  <- ftld_any == 1
na_if0 <- function(flag, value) ifelse(flag == 1, value, NA_character_)
ftld_subtype <- case_when(ftld_tau == 1 ~ "FTLD-Tau", ftld_tdp == 1 ~ "FTLD-TDP",
                          ftld_fus == 1 | ftld_other == 1 ~ "FTLD-Other/NOS", TRUE ~ NA_character_)
ftld_tau_specific  <- na_if0(ftld_tau, sample(c("PSP", "CBD", "AGD", "Pick's Disease", "Other Tauopathy"), N, TRUE))
ftld_tdp_mnd_status <- case_when(ftld_mnd == 1 ~ "FTLD-MND", ftld_tdp == 1 ~ "FTLD-TDP without MND", TRUE ~ NA_character_)
ftld_other_specific <- na_if0(ftld_fus, "FTLD-FUS")
ftld_adnc_group <- case_when(ftld_any == 1 & ad == "Sev ADNC" ~ "FTLD with Sev ADNC",
                             ftld_any == 1 & ad == "Int ADNC" ~ "FTLD with Int ADNC",
                             ftld_any == 1 ~ "FTLD without ADNC", TRUE ~ NA_character_)
ftld_subtype_adnc <- case_when(ftld_tau == 1 & sev >= 2 ~ "FTLD-Tau with ADNC", ftld_tau == 1 ~ "FTLD-Tau without ADNC",
                               ftld_tdp == 1 & sev >= 2 ~ "FTLD-TDP with ADNC", ftld_tdp == 1 ~ "FTLD-TDP without ADNC",
                               TRUE ~ NA_character_)

og_levels <- c("Low Path", "Int ADNC", "Sev ADNC", "FTLD", "LBD", "LATE", "Other", "Mixed ADNC", "Unclassified")
has_copath <- (lbd_stage != "Absent") | (late_stage != "Absent") | FTLD_present | Other_present
Overall_Group <- factor(case_when(
  ad == "Low ADNC" & FTLD_present ~ "FTLD", ad == "Low ADNC" & lbd_stage != "Absent" ~ "LBD",
  ad == "Low ADNC" & late_stage != "Absent" ~ "LATE", ad == "Low ADNC" & Other_present ~ "Other",
  ad == "Low ADNC" ~ "Low Path",
  ad == "Int ADNC" & has_copath ~ "Mixed ADNC", ad == "Int ADNC" ~ "Int ADNC",
  ad == "Sev ADNC" & has_copath ~ "Mixed ADNC", ad == "Sev ADNC" ~ "Sev ADNC",
  TRUE ~ "Unclassified"), levels = og_levels)

# Clinical AD etiology
NACCALZD <- rbinom(N, 1, ifelse(ad == "Low ADNC", 0.20, 0.85))
int <- as.integer(ad == "Int ADNC")
ID <- sprintf("DEMO%04d", seq_len(N))

analysis <- tibble(
  ID, VISEQ, AGEOFONSET, APOE, Sex, Race, Educ, AgeAtVisit, CDRSum, PersCare, MMSE, NPDAGE,
  ad, Braak_grouped, NEUR, lbd_stage, late_stage, dx, dx2,
  ptau217, ptau181, GFAP, NfL, log_ptau217, log_ptau181, log_GFAP, log_NfL,
  time_to_death, apoe4, NACCALZD, int, ftld_tau, ftld_tdp, ftld_fus, ftld_mnd, ftld_other, ftld_any,
  FTLD_present, Other_present, ftld_subtype, ftld_tau_specific, ftld_tdp_mnd_status,
  ftld_other_specific, ftld_adnc_group, ftld_subtype_adnc, Overall_Group)

# NACC ordinal pathology codes and FTLD
np_codes <- c("NPALSMND", "NACCAVAS", "NACCINF", "NACCMICR", "NACCHEM", "NACCARTE", "NACCAMY")
for (nm in np_codes) analysis[[nm]] <- sample(0:3, N, replace = TRUE, prob = c(0.7, 0.15, 0.1, 0.05))
flag_cols <- c("flag_picks", "flag_psp", "flag_cbd", "flag_agd", "flag_ftld_fus", "flag_nifid",
               "flag_bibd", "flag_als", "flag_tdp43_spinal", "flag_tbi_acute", "flag_tbi_chronic",
               "flag_hippo_sclerosis")
for (nm in flag_cols) analysis[[nm]] <- rbinom(N, 1, 0.05)

saveRDS(analysis, file.path(out_dir, "demo_cohort.rds"))
write_csv(analysis, file.path(out_dir, "demo_cohort.csv"))

# Longitudinal: expand each participant to multiple visits
long <- analysis %>%
  rowwise() %>% mutate(n_visits = sample(1:4, 1, prob = c(0.40, 0.35, 0.15, 0.10))) %>%  # ~2 visits/subject
  uncount(n_visits, .id = "visit") %>%
  group_by(ID) %>%
  mutate(
    nv = n(), years_back = (nv - visit) * runif(1, 1.5, 4),
    VISEQ = visit, AgeAtVisit = round(AgeAtVisit - years_back, 1),
    time_to_death = round(time_to_death + years_back, 2),
    MMSE = round(pmin(30, MMSE + years_back * runif(1, 0.4, 1.2))),
    dx2_sev = pmax(1L, as.integer(dx2) - (years_back > 2.5) - (years_back > 6)),
    dx2 = factor(dx2_levels[dx2_sev], levels = dx2_levels),
    dx = factor(ifelse(dx2 == "Cognitively Normal", "Control",
                ifelse(dx2 == "MCI", "MCI", as.character(dx))), levels = dx_levels),
    ptau217 = ptau217 * exp(-years_back * 0.05), ptau181 = ptau181 * exp(-years_back * 0.04),
    GFAP = GFAP * exp(-years_back * 0.04), NfL = NfL * exp(-years_back * 0.05),
    log_ptau217 = log(ptau217), log_ptau181 = log(ptau181), log_GFAP = log(GFAP), log_NfL = log(NfL),
    death_date = as.Date("2024-01-01"), IDATE = death_date - round(time_to_death * 365.25)) %>%
  ungroup() %>% select(-years_back, -nv, -visit, -dx2_sev)

saveRDS(long, file.path(out_dir, "demo_cohort_long.rds"))
cat(sprintf("Wrote %d participants and %d longitudinal visits to %s\n",
            nrow(analysis), nrow(long), normalizePath(out_dir)))
