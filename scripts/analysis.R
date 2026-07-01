########################################
## Air pollution and birthweight
## Generation C cohort analysis
##
## This script estimates associations
## between prenatal PM2.5 exposure and
## BW-for-GA z-scores, including effect
## modification by neighborhood-level
## structural disadvantage.
##
## Raw participant-level data are not
## included in this repository.
########################################

# Required external package:
# install.packages("remotes")
# remotes::install_github("ki-tools/growthstandards")
#
# Additional packages:
# install.packages(c("zoo", "EValue", "tableone", "WeightIt", "here"))

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(fst)
  library(stringr)
  library(mice)
  library(bdlim)
  library(ggplot2)
  library(janitor)
  library(forcats)
  library(growthstandards)
  library(zoo)
  library(scales)
  library(flextable)
  library(officer)
  library(ggtext)
  library(gt)
  library(grid)
  library(EValue)
  library(tableone)
  library(WeightIt)
})

options(dplyr.summarise.inform = FALSE)
set.seed(500)

########################################
## Helper functions
########################################

first_nonNA <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

season_from_date <- function(x) {
  m <- month(x)
  case_when(
    m %in% 3:5       ~ "Spring",
    m %in% 6:8       ~ "Summer",
    m %in% 9:11      ~ "Autumn",
    m %in% c(12, 1, 2) ~ "Winter",
    TRUE             ~ NA_character_
  )
}

med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  sprintf("%.1f (%.1f, %.1f)",
          median(x), quantile(x, 0.25), quantile(x, 0.75))
}

n_pct <- function(x) {
  n   <- sum(x, na.rm = TRUE)
  pct <- round(100 * mean(x, na.rm = TRUE), 0)
  sprintf("%s (%s%%)", comma(n), pct)
}

wilcox_p <- function(var, group, data) {
  d <- data %>%
    filter(!is.na(.data[[var]]), !is.na(.data[[group]])) %>%
    mutate(grp = droplevels(as.factor(.data[[group]])))
  if (nrow(d) < 10 || nlevels(d$grp) < 2) return(NA_real_)
  wilcox.test(as.formula(paste(var, "~ grp")), data = d, exact = FALSE)$p.value
}

chisq_p <- function(var, group, data) {
  tbl <- table(data[[var]], data[[group]])
  if (any(dim(tbl) < 2)) return(NA_real_)
  suppressWarnings(chisq.test(tbl)$p.value)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

make_row <- function(label, overall, high, low, p = "",
                     is_header = FALSE, indent = FALSE) {
  data.frame(
    Characteristic = label,
    Overall        = overall,
    Low_SREI       = low,
    High_SREI      = high,
    p_value        = as.character(p),
    is_header      = is_header,
    indent         = indent,
    stringsAsFactors = FALSE
  )
}

########################################
## Paths
########################################

dir_data    <- here::here("data")
dir_outputs <- here::here("outputs")
dir_tables  <- file.path(dir_outputs, "tables")
dir_figures <- file.path(dir_outputs, "figures")

dir.create(dir_outputs, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

path_predictions <- file.path(dir_data, "predictions_2025-10-10.fst")
path_addr_csv    <- file.path(dir_data, "merged_address_data_2025-10-22.csv")
path_birth_csv   <- file.path(dir_data, "birth_data_2025-10-14.csv")
path_hdp_csv     <- file.path(dir_data, "hypertensive_disorders_2024-10-16.csv")
path_srei_csv    <- file.path(dir_data, "mean_srei_genc.csv")

########################################
## Read data
########################################

# Address data
data_addr <- read_csv(
  path_addr_csv,
  col_types = cols(
    .default      = col_guess(),
    id            = col_character(),
    inputrow      = col_character(),
    birthdate     = col_date(format = ""),
    lmp           = col_date(format = ""),
    start_poll    = col_date(format = ""),
    end_poll      = col_date(format = ""),
    conception90d = col_date(format = "")
  )
) %>%
  clean_names()

# Predictions: normalize units from 10 ng/m³ to µg/m³
predictions <- read_fst(path_predictions, as.data.table = FALSE) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(
    inputrow   = as.character(inputrow),
    date       = as_date(date),
    pm25_ug_m3 = pm25_10ng_m3 * 0.01
  )

# Birth and covariate data
birth_data <- read_csv(
  path_birth_csv,
  col_types = cols(
    .default           = col_guess(),
    subject_id         = col_character(),
    birthdate          = col_character(),
    lmp                = col_character(),
    tri_1_end          = col_character(),
    tri_2_end          = col_character(),
    tri_3_end          = col_character(),
    gestationalagedays = col_double(),
    birthweight_grams  = col_double()
  )
) %>%
  clean_names() %>%
  mutate(
    birthweight_grams = dplyr::na_if(birthweight_grams, -99),
    birthdate         = suppressWarnings(mdy(birthdate)),
    lmp               = suppressWarnings(mdy(lmp)),
    tri_1_end         = suppressWarnings(mdy(tri_1_end)),
    tri_2_end         = suppressWarnings(mdy(tri_2_end)),
    tri_3_end         = suppressWarnings(mdy(tri_3_end))
  )

# Hypertensive disorder data
hyp_raw <- suppressMessages(
  read_csv(path_hdp_csv, show_col_types = FALSE)
) %>%
  clean_names() %>%
  mutate(subject_id = as.character(subject_id))

needed_hdp <- c(
  "gestationalhypertension", "gestationaldiabetes", "preeclampsia",
  "emr_ghtn_current", "emr_gdm_current", "emr_pec_current",
  "gest_htn_preg", "dm_preg", "gest_dm_preg", "gest_dm_pp", "pec_preg"
)

missing_hdp_cols <- setdiff(needed_hdp, names(hyp_raw))
if (length(missing_hdp_cols) > 0) {
  stop("Missing HDP columns: ", paste(missing_hdp_cols, collapse = ", "))
}

hyp_raw <- hyp_raw %>%
  mutate(
    across(
      all_of(needed_hdp),
      ~ as.integer(as.character(.) %in% c("1", "Yes", "TRUE", "True"))
    )
  ) %>%
  mutate(
    any_hypertensive_disorder_preg = if_else(
      rowSums(across(all_of(needed_hdp)), na.rm = TRUE) > 0, 1L, 0L
    ),
    diabetes_preg = if_else(
      rowSums(across(c(
        gestationaldiabetes, emr_gdm_current, dm_preg,
        gest_dm_preg, gest_dm_pp
      )), na.rm = TRUE) > 0, 1L, 0L
    ),
    preeclampsia_preg = if_else(
      rowSums(across(c(preeclampsia, emr_pec_current, pec_preg)), na.rm = TRUE) > 0,
      1L, 0L
    ),
    hypertension_preg = if_else(
      rowSums(across(c(
        gestationalhypertension, emr_ghtn_current, gest_htn_preg
      )), na.rm = TRUE) > 0,
      1L, 0L
    )
  ) %>%
  mutate(
    any_hypertensive_disorder_preg_fac = factor(
      any_hypertensive_disorder_preg, levels = c(1, 0), labels = c("Yes", "No")
    ),
    diabetes_preg_fac = factor(
      diabetes_preg, levels = c(1, 0), labels = c("Yes", "No")
    ),
    preeclampsia_preg_fac = factor(
      preeclampsia_preg, levels = c(1, 0), labels = c("Yes", "No")
    ),
    hypertension_preg_fac = factor(
      hypertension_preg, levels = c(1, 0), labels = c("Yes", "No")
    )
  )

hyp_raw2 <- hyp_raw %>%
  dplyr::select(
    subject_id,
    any_hypertensive_disorder_preg_fac,
    diabetes_preg_fac,
    preeclampsia_preg_fac,
    hypertension_preg_fac
  )

# SREI data
srei_raw <- read_csv(path_srei_csv, show_col_types = FALSE) %>%
  clean_names() %>%
  dplyr::rename(subject_id = id) %>%
  mutate(subject_id = as.character(subject_id))

########################################
## Quick sanity checks
########################################

cat("birth_data IDs:        ", n_distinct(birth_data$subject_id), "\n")
cat("address IDs:           ", n_distinct(data_addr$id), "\n")
cat("address inputrows:     ", n_distinct(data_addr$inputrow), "\n")
cat("predictions inputrows: ", n_distinct(predictions$inputrow), "\n")
cat(
  "predictions dates:     ",
  paste(range(predictions$date, na.rm = TRUE), collapse = " to "),
  "\n\n"
)

########################################
## Build long analytic file
########################################

# Join births to address history
df_birth_addr <- birth_data %>%
  mutate(subject_id = as.character(subject_id)) %>%
  left_join(
    data_addr %>% mutate(id = as.character(id)),
    by = c("subject_id" = "id")
  )

# Join predictions; restrict PM2.5 to residential address window
df_long <- df_birth_addr %>%
  left_join(
    predictions %>% dplyr::select(inputrow, date, pm25_ug_m3),
    by = "inputrow"
  ) %>%
  mutate(
    lmp_calc = as_date(birthdate.x) - as.numeric(gestationalagedays),
    lmp      = coalesce(as_date(lmp.x), lmp_calc),
    use_air  = case_when(
      !is.na(start_poll) & !is.na(end_poll) &
        !is.na(date) & date >= start_poll & date <= end_poll ~ pm25_ug_m3,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(subject_id = as.character(subject_id)) %>%
  left_join(hyp_raw2, by = "subject_id") %>%
  left_join(srei_raw, by = "subject_id")

# Gestational timing variables
df_long <- df_long %>%
  mutate(
    gaweek        = as.numeric(difftime(date, lmp, units = "weeks")),
    week          = floor(gaweek),
    trimester_num = case_when(
      gaweek < 13                ~ 1,
      gaweek >= 13 & gaweek < 27 ~ 2,
      gaweek >= 27               ~ 3,
      TRUE                       ~ NA_real_
    ),
    trimester = factor(
      trimester_num,
      levels = c(1, 2, 3),
      labels = c("First (0-13w)", "Second (14-27w)", "Third (28-42w)")
    ),
    birth_year = year(birthdate.x),
    year_fac   = factor(
      case_when(
        birthdate.x >= ymd("2020-01-01") & birthdate.x <= ymd("2020-12-31") ~ "2020",
        birthdate.x >= ymd("2021-01-01") & birthdate.x <= ymd("2021-12-31") ~ "2021",
        birthdate.x >= ymd("2022-01-01") & birthdate.x <= ymd("2022-12-31") ~ "2022",
        TRUE ~ NA_character_
      ),
      levels = c("2020", "2021", "2022")
    ),
    season = factor(
      season_from_date(lmp),
      levels = c("Summer", "Spring", "Autumn", "Winter")
    )
  )

# Restrict to rows with exposure data
df_overlap <- df_long %>% filter(!is.na(use_air))

cat("IDs with ≥1 overlapping exposure row:", n_distinct(df_overlap$subject_id), "\n")

########################################
## Exposure summaries
########################################

# Pregnancy-average PM2.5 (used in sensitivity analysis)
preg_avg <- df_overlap %>%
  group_by(subject_id) %>%
  summarise(preg_avg_pm25 = mean(use_air, na.rm = TRUE), .groups = "drop")

# Weekly means: restrict to gestational weeks 1–37 (term births only)
week_avg <- df_overlap %>%
  mutate(week = pmin(pmax(floor(gaweek), 1L), 37L)) %>%
  filter(!is.na(week)) %>%
  group_by(subject_id, week) %>%
  summarise(week_avg = mean(use_air, na.rm = TRUE), .groups = "drop")

########################################
## Define final analytic sample
########################################

# Step 1: count observed weeks per subject
week_counts <- week_avg %>%
  group_by(subject_id) %>%
  summarise(n_weeks_nonNA = sum(!is.na(week_avg)), .groups = "drop")

# Step 2: subjects with all 37 weeks observed
ids_full37 <- week_counts %>%
  filter(n_weeks_nonNA == 37) %>%
  pull(subject_id)

# Step 3: subjects with 35–36 weeks — interpolate gaps of ≤2 weeks
ids_near_complete <- week_counts %>%
  filter(n_weeks_nonNA >= 35, n_weeks_nonNA < 37) %>%
  pull(subject_id)

cat("Subjects with all 37 weeks:          ", length(ids_full37), "\n")
cat("Subjects recovered via interpolation:", length(ids_near_complete), "\n")

# Keep complete subjects as-is
week_avg_complete <- week_avg %>%
  filter(subject_id %in% ids_full37)

week_avg_interpolated <- week_avg %>%
  filter(subject_id %in% ids_near_complete) %>%
  group_by(subject_id) %>%
  tidyr::complete(week = 1:37) %>%
  arrange(subject_id, week) %>%
  mutate(
    week_avg = zoo::na.approx(week_avg, na.rm = FALSE, maxgap = 2),
    week_avg = zoo::na.locf(week_avg, na.rm = FALSE, maxgap = 2),
    week_avg = zoo::na.locf(week_avg, na.rm = FALSE, maxgap = 2, fromLast = TRUE)
  ) %>%
  ungroup()

week_avg_imputed <- bind_rows(week_avg_complete, week_avg_interpolated)

week_counts_final <- week_avg_imputed %>%
  group_by(subject_id) %>%
  summarise(n_weeks_nonNA = sum(!is.na(week_avg)), .groups = "drop")

ids_full37_final <- week_counts_final %>%
  filter(n_weeks_nonNA == 37) %>%
  pull(subject_id)

cat("After interpolation:", length(ids_full37_final), "\n")

# Confirm near-complete subjects are now all at 37
week_counts_final %>%
  filter(subject_id %in% ids_near_complete) %>%
  pull(n_weeks_nonNA) %>%
  table()

# Step 4: eligibility criteria
ids_bwt_ok <- birth_data %>%
  filter(!is.na(birthweight_grams)) %>%
  distinct(subject_id) %>%
  pull(subject_id)

# Term births: ≥37 completed weeks (≥259 days)
ids_term <- birth_data %>%
  filter(!is.na(gestationalagedays), gestationalagedays >= 259) %>%
  distinct(subject_id) %>%
  pull(subject_id)

# Singletons: exclude confirmed twins; retain NA
ids_singleton <- birth_data %>%
  mutate(twin_chr = as.character(twin)) %>%
  filter(is.na(twin_chr) | twin_chr %in% c("0", "No", "NO", "no")) %>%
  distinct(subject_id) %>%
  pull(subject_id)

# Final analytic sample
keep_ids <- Reduce(intersect, list(
  ids_full37_final,
  ids_bwt_ok,
  ids_term,
  ids_singleton
))

cat("\nFinal sample counts:\n")
cat("  Complete exposure (37 weeks): ", length(ids_full37_final), "\n")
cat("  Non-missing birthweight:      ", length(ids_bwt_ok), "\n")
cat("  Term births (≥259 days):      ", length(ids_term), "\n")
cat("  Singletons:                   ", length(ids_singleton), "\n")
cat("  Final kept IDs:               ", length(keep_ids), "\n")

########################################
## Weekly exposure matrix (wide format)
########################################

pm25_week_wide <- week_avg_imputed %>%
  tidyr::pivot_wider(
    names_from   = week,
    values_from  = week_avg,
    names_prefix = "week_avg."
  )

pm25_week_37 <- pm25_week_wide %>%
  filter(subject_id %in% keep_ids) %>%
  arrange(subject_id)

dim(pm25_week_37)

########################################
## Outcome
########################################

bwt <- birth_data %>%
  distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(subject_id, birthweight_grams, gestationalagedays) %>%
  filter(subject_id %in% keep_ids) %>%
  arrange(subject_id)

summary(bwt$birthweight_grams)

########################################
## Subject-level covariates
########################################

covar <- df_overlap %>%
  filter(subject_id %in% keep_ids) %>%
  group_by(subject_id) %>%
  summarise(
    preg_pos                           = first_nonNA(preg_pos),
    lmp                                = first_nonNA(lmp),
    gender                             = first_nonNA(gender),
    raceethnicitycombined              = first_nonNA(raceethnicitycombined),
    season                             = first_nonNA(season),
    year                               = first_nonNA(year_fac),
    parity_cat                         = first_nonNA(parity_cat),
    prepregnancybmi                    = first_nonNA(prepregnancybmi),
    maternalage                        = first_nonNA(maternalage),
    maternalage_cat                    = first_nonNA(maternalage_cat),
    maritalstatus                      = first_nonNA(maritalstatus),
    insurance_cat                      = first_nonNA(insurance_cat),
    ageatdelivery                      = first_nonNA(ageatdelivery),
    postalcode                         = first_nonNA(postalcode),
    nyc_borough                        = first_nonNA(nyc_borough),
    any_hypertensive_disorder_preg_fac = first_nonNA(any_hypertensive_disorder_preg_fac),
    diabetes_preg_fac                  = first_nonNA(diabetes_preg_fac),
    preeclampsia_preg_fac              = first_nonNA(preeclampsia_preg_fac),
    hypertension_preg_fac              = first_nonNA(hypertension_preg_fac),
    .groups = "drop"
  ) %>%
  left_join(
    srei_raw %>%
      dplyr::select(subject_id, srei, srei_quartile, srei_75_split, srei_25_split),
    by = "subject_id"
  ) %>%
  mutate(preg_pos = replace_na(as.character(preg_pos), "No"))

########################################
## Multiple imputation (MICE)
########################################

covar_mice <- covar %>%
  dplyr::select(
    subject_id, preg_pos, gender, season, year,
    parity_cat, prepregnancybmi, maternalage,
    maritalstatus, insurance_cat,
    any_hypertensive_disorder_preg_fac,
    srei, srei_quartile, srei_75_split, srei_25_split
  )

missing_counts <- sapply(covar_mice, function(x) sum(is.na(x)))
print(missing_counts)

meth_subj <- make.method(covar_mice)
pred_subj <- make.predictorMatrix(covar_mice)

# Never impute subject_id
meth_subj["subject_id"] <- ""
pred_subj["subject_id", ] <- 0
pred_subj[, "subject_id"] <- 0

missing_var <- names(covar_mice)[colSums(is.na(covar_mice)) > 0]
meth_subj[missing_var] <- "rf"
meth_subj[setdiff(names(covar_mice), c("subject_id", missing_var))] <- ""

# Generate 5 imputed datasets; primary analysis uses draw 1
set.seed(500)
imp_subj <- mice(
  covar_mice,
  m = 5,
  method = meth_subj,
  predictorMatrix = pred_subj,
  maxit = 5,
  seed = 500,
  printFlag = FALSE
)

covar_imputed <- complete(imp_subj, 1) %>%
  arrange(subject_id)

missing_counts_after <- sapply(covar_imputed, function(x) sum(is.na(x)))
print(missing_counts_after)

########################################
## Assemble final analytic dataset
########################################

dat <- pm25_week_37 %>%
  left_join(bwt,           by = "subject_id") %>%
  left_join(covar_imputed, by = "subject_id") %>%
  arrange(subject_id)

cat("\nDimension checks (all should match):\n")
cat("  pm25_week_37 rows: ", nrow(pm25_week_37), "\n")
cat("  bwt rows:          ", nrow(bwt), "\n")
cat("  covar_imputed rows:", nrow(covar_imputed), "\n")
cat("  dat rows:          ", nrow(dat), "\n")

########################################
## Recode analysis variables
########################################

# SREI grouping (primary modifier: 75th-percentile split)
dat$srei1 <- factor(
  dat$srei_75_split,
  levels = c("75th_High", "Below 75th"),
  labels = c("≥75th percentile", "<75th percentile")
)

table(dat$srei1, useNA = "always")

# Parity
dat$parity <- factor(
  dat$parity_cat,
  levels = c("0", "1"),
  labels = c("Nulliparous", "Multiparous")
)

# Sex (used only for z-score computation, NOT as a model covariate)
dat$sex <- factor(
  as.character(dat$gender),
  levels = c("1", "2"),
  labels = c("Female", "Male")
)

dat$sex_clean <- case_when(
  dat$sex == "Male"   ~ "Male",
  dat$sex == "Female" ~ "Female",
  TRUE ~ NA_character_
)

# Insurance: collapse Self-pay into Private; reference = Private
dat$insurance <- factor(
  as.character(dat$insurance_cat),
  levels = c("1", "2", "3"),
  labels = c("Private", "Public", "Self-pay")
) %>%
  fct_collapse(Private = c("Private", "Self-pay"), Public = "Public") %>%
  fct_relevel("Private")

# Marital status: collapse to four categories; reference = Significant Other
dat$maritalstatus <- factor(
  as.character(dat$maritalstatus),
  levels = as.character(1:7),
  labels = c(
    "Divorced", "Legally Separated/Separated",
    "Married/Civil Union", "Significant Other/Life Partner",
    "Single", "Unknown/Other", "Widowed"
  )
) %>%
  fct_collapse(
    `Legally Separated` = c("Divorced", "Legally Separated/Separated"),
    `Significant Other` = c("Married/Civil Union", "Significant Other/Life Partner"),
    Other               = c("Widowed", "Unknown/Other")
  ) %>%
  fct_relevel("Significant Other")

# Season and year
dat$season <- factor(dat$season, levels = c("Summer", "Spring", "Autumn", "Winter"))
dat$year   <- factor(as.character(dat$year), levels = c("2020", "2021", "2022"))

# Continuous covariates
dat$bmi <- as.numeric(dat$prepregnancybmi)
dat$age <- as.numeric(dat$maternalage)

# Other
dat$anyhypdis <- factor(dat$any_hypertensive_disorder_preg_fac)
dat$pos       <- factor(as.character(dat$preg_pos), levels = c("Yes", "No"))

########################################
## BW-for-GA z-score (INTERGROWTH-21st)
########################################

dat <- dat %>%
  mutate(
    bw_kg     = birthweight_grams / 1000,
    ga_days   = gestationalagedays,
    bw_z_intergrowth = igb_wtkg2zscore(
      gagebrth = ga_days,
      wtkg     = bw_kg,
      sex      = sex_clean
    ),
    bw_centile_intergrowth = igb_wtkg2centile(
      gagebrth = ga_days,
      wtkg     = bw_kg,
      sex      = sex_clean
    )
  )

summary(dat$bw_z_intergrowth)
table(is.na(dat$bw_z_intergrowth))

# Remove any rows with missing primary outcome
dat <- dat %>% filter(!is.na(bw_z_intergrowth))

########################################
## Prepare BDLIM inputs
########################################

wk_cols <- paste0("week_avg.", 1:37)

X  <- dat %>% dplyr::select(all_of(wk_cols)) %>% as.matrix()
Yz <- dat$bw_z_intergrowth
G  <- as.factor(dat$srei1)

stopifnot(all(wk_cols %in% names(dat)))
stopifnot(!anyDuplicated(dat$subject_id))
stopifnot(length(Yz) == nrow(X))
stopifnot(length(G) == nrow(X))

cat("\nBDLIM input dimensions:\n")
cat("  Outcome (Yz):  ", length(Yz), "\n")
cat("  Exposure (X):  ", nrow(X), "x", ncol(X), "\n")
cat("  Group (G):     ", length(G), "\n")

########################################
## Primary analyses: BDLIM
########################################

## Unadjusted model
set.seed(123)
fit <- bdlim4(
  y        = Yz,
  exposure = X,
  covars   = dat[, 0, drop = FALSE],
  group    = G,
  df       = 4,
  nits     = 50000,
  nthin    = 10,
  parallel = TRUE
)

fit
print(summary(fit))
plot(summary(fit))

## Adjusted model
set.seed(123)
fit_adj <- bdlim4(
  y        = Yz,
  exposure = X,
  covars   = dat[, c(
    "bmi", "age", "parity", "maritalstatus",
    "insurance", "season", "year"
  ), drop = FALSE],
  group    = G,
  df       = 4,
  nits     = 50000,
  nthin    = 10,
  parallel = TRUE
)

fit_adj
print(summary(fit_adj))
plot(summary(fit_adj))

########################################
## Sensitivity: results across imputed datasets
## (Supplementary eTable 1)
########################################

cum_list <- vector("list", 5)

dat_base <- dat %>%
  dplyr::select(
    subject_id,
    all_of(wk_cols),
    birthweight_grams,
    gestationalagedays,
    bw_z_intergrowth
  )

for (i in 1:5) {
  cat("Running imputation", i, "\n")
  
  covar_imp_i <- complete(imp_subj, i) %>%
    arrange(subject_id) %>%
    mutate(
      srei1 = factor(
        srei_75_split,
        levels = c("75th_High", "Below 75th"),
        labels = c("≥75th percentile", "<75th percentile")
      ),
      parity = factor(
        parity_cat,
        levels = c("0", "1"),
        labels = c("Nulliparous", "Multiparous")
      ),
      insurance = factor(
        as.character(insurance_cat),
        levels = c("1", "2", "3"),
        labels = c("Private", "Public", "Self-pay")
      ) %>%
        fct_collapse(
          Private = c("Private", "Self-pay"),
          Public  = "Public"
        ) %>%
        fct_relevel("Private"),
      maritalstatus = factor(
        as.character(maritalstatus),
        levels = as.character(1:7),
        labels = c(
          "Divorced", "Legally Separated/Separated",
          "Married/Civil Union", "Significant Other/Life Partner",
          "Single", "Unknown/Other", "Widowed"
        )
      ) %>%
        fct_collapse(
          `Legally Separated` = c("Divorced", "Legally Separated/Separated"),
          `Significant Other` = c("Married/Civil Union", "Significant Other/Life Partner"),
          Other               = c("Widowed", "Unknown/Other")
        ) %>%
        fct_relevel("Significant Other"),
      season = factor(season, levels = c("Summer", "Spring", "Autumn", "Winter")),
      year   = factor(as.character(year), levels = c("2020", "2021", "2022")),
      bmi    = as.numeric(prepregnancybmi),
      age    = as.numeric(maternalage)
    ) %>%
    dplyr::select(
      subject_id, srei1, bmi, age, parity, maritalstatus, insurance, season, year
    )
  
  dat_i <- dat_base %>%
    left_join(covar_imp_i, by = "subject_id") %>%
    arrange(subject_id)
  
  Y_i <- dat_i$bw_z_intergrowth
  X_i <- dat_i %>% dplyr::select(all_of(wk_cols)) %>% as.matrix()
  G_i <- as.factor(dat_i$srei1)
  
  fit_i <- bdlim4(
    y        = Y_i,
    exposure = X_i,
    covars   = dat_i[, c(
      "bmi", "age", "parity", "maritalstatus",
      "insurance", "season", "year"
    ), drop = FALSE],
    group    = G_i,
    df       = 4,
    nits     = 20000,
    nthin    = 10,
    parallel = FALSE
  )
  
  cum_i <- summary(fit_i)$cumulative %>%
    mutate(imputation = i)
  
  cum_list[[i]] <- cum_i
}

cum_results <- bind_rows(cum_list)
print(cum_results)

supp_table_full <- cum_results %>%
  mutate(
    mean   = round(mean, 2),
    q2.5   = round(q2.5, 2),
    q97.5  = round(q97.5, 2),
    pr_gr0 = round(pr_gr0, 2)
  ) %>%
  transmute(
    Imputation   = imputation,
    `SREI group` = group,
    Mean         = mean,
    `95% CrI`    = paste0(q2.5, ", ", q97.5),
    `Pr(>0)`     = pr_gr0
  ) %>%
  arrange(Imputation, `SREI group`)

print(supp_table_full)

ft_supp <- flextable(supp_table_full) %>%
  set_caption(
    "Supplementary Table 1. Cumulative associations between prenatal PM2.5 exposure and BW-for-GA z-scores across multiply imputed datasets."
  ) %>%
  autofit()

doc_supp <- read_docx() %>%
  body_add_flextable(ft_supp)

print(
  doc_supp,
  target = file.path(dir_tables, "Supplementary_Table1_Imputed.docx")
)

########################################
## Sensitivity analysis 1: linear regression with
## pregnancy-average PM2.5, incl. PM2.5 x SREI interaction
########################################

preg_avg_sens_int <- preg_avg %>%
  filter(subject_id %in% keep_ids) %>%
  arrange(subject_id) %>%
  left_join(
    dat %>%
      dplyr::select(
        subject_id, bw_z_intergrowth, age, bmi,
        parity, maritalstatus, insurance, season, year, srei1
      ),
    by = "subject_id"
  )

# Unadjusted
lm_z_unadj <- lm(bw_z_intergrowth ~ preg_avg_pm25, data = preg_avg_sens_int)
summary(lm_z_unadj)

# Adjusted
lm_z_adj <- lm(
  bw_z_intergrowth ~ preg_avg_pm25 + bmi + parity + age +
    maritalstatus + insurance + season + year,
  data = preg_avg_sens_int
)
summary(lm_z_adj)

# Adjusted, with PM2.5 x SREI interaction
lm_z_adj_int <- lm(
  bw_z_intergrowth ~ preg_avg_pm25 * srei1 + bmi + parity + age +
    maritalstatus + insurance + season + year,
  data = preg_avg_sens_int
)
summary(lm_z_adj_int)

########################################
## Sensitivity analysis 2: 2021-2022 births only
## (excluding pandemic-onset 2020)
########################################

dat_2122 <- dat %>% filter(year %in% c("2021", "2022"))
cat("Sample size (2021-2022 births):", nrow(dat_2122), "\n")

Y_2122 <- dat_2122$bw_z_intergrowth
X_2122 <- dat_2122 %>% dplyr::select(all_of(wk_cols)) %>% as.matrix()
G_2122 <- as.factor(dat_2122$srei1)

stopifnot(
  length(Y_2122) == nrow(X_2122),
  length(G_2122) == nrow(X_2122)
)

## Unadjusted
set.seed(123)
fit_2122_unadj <- bdlim4(
  y        = Y_2122,
  exposure = X_2122,
  covars   = dat_2122[, 0, drop = FALSE],
  group    = G_2122,
  df       = 4,
  nits     = 50000,
  nthin    = 10,
  parallel = TRUE
)
print(summary(fit_2122_unadj))
plot(summary(fit_2122_unadj))

## Adjusted (year dropped: only 2021-2022 remain)
set.seed(123)
fit_2122_adj <- bdlim4(
  y        = Y_2122,
  exposure = X_2122,
  covars   = dat_2122[, c(
    "bmi", "age", "parity", "maritalstatus",
    "insurance", "season"
  ), drop = FALSE],
  group    = G_2122,
  df       = 4,
  nits     = 50000,
  nthin    = 10,
  parallel = TRUE
)
print(summary(fit_2122_adj))
plot(summary(fit_2122_adj))

supp_table_2122 <- bind_rows(
  summary(fit_2122_unadj)$cumulative %>% mutate(Model = "Unadjusted"),
  summary(fit_2122_adj)$cumulative   %>% mutate(Model = "Adjusted")
) %>%
  transmute(
    Model,
    `SREI group` = group,
    Mean     = round(mean, 4),
    Median   = round(median, 4),
    SD       = round(sd, 4),
    `2.5%`   = round(q2.5, 4),
    `97.5%`  = round(q97.5, 4),
    `Pr(>0)` = round(pr_gr0, 2)
  ) %>%
  arrange(Model, desc(`SREI group`))

print(supp_table_2122)

ft_2122 <- flextable(supp_table_2122) %>%
  set_caption(
    "Supplementary Table 2. Cumulative associations between prenatal PM2.5 exposure and BW-for-GA z-scores, restricted to 2021-2022 births, by SREI group."
  ) %>%
  autofit()

doc_2122 <- read_docx() %>%
  body_add_flextable(ft_2122)

print(doc_2122, target = file.path(dir_tables, "Supplementary_Table2_2021_2022.docx"))

########################################
## Sensitivity analysis 3: complete cases only
## (fully observed 37-week exposure series)
########################################

keep_ids_cc <- setdiff(keep_ids, ids_near_complete)
dat_cc <- dat %>% filter(subject_id %in% keep_ids_cc)
cat("Complete-case sample size:", nrow(dat_cc), "\n")

X_cc  <- dat_cc %>% dplyr::select(all_of(wk_cols)) %>% as.matrix()
Yz_cc <- dat_cc$bw_z_intergrowth
G_cc  <- as.factor(dat_cc$srei1)

set.seed(123)
fit_cc_unadj <- bdlim4(
  y        = Yz_cc,
  exposure = X_cc,
  covars   = dat_cc[, 0, drop = FALSE],
  group    = G_cc,
  df       = 4,
  nits     = 50000,
  nthin    = 10,
  parallel = TRUE
)
summary(fit_cc_unadj)$cumulative

set.seed(123)
fit_cc_adj <- bdlim4(
  y        = Yz_cc,
  exposure = X_cc,
  covars   = dat_cc[, c(
    "bmi", "age", "parity", "maritalstatus",
    "insurance", "season", "year"
  ), drop = FALSE],
  group    = G_cc,
  df       = 4,
  nits     = 50000,
  nthin    = 10,
  parallel = TRUE
)
summary(fit_cc_adj)$cumulative

supp_table_cc <- bind_rows(
  summary(fit_cc_unadj)$cumulative %>% mutate(Model = "Unadjusted"),
  summary(fit_cc_adj)$cumulative   %>% mutate(Model = "Adjusted")
) %>%
  transmute(
    Model,
    `SREI group` = group,
    Mean     = round(mean, 4),
    Median   = round(median, 4),
    SD       = round(sd, 4),
    `2.5%`   = round(q2.5, 4),
    `97.5%`  = round(q97.5, 4),
    `Pr(>0)` = round(pr_gr0, 2)
  ) %>%
  arrange(Model, desc(`SREI group`))

print(supp_table_cc)

ft_cc <- flextable(supp_table_cc) %>%
  set_caption(
    "Supplementary Table 3. Cumulative associations between prenatal PM2.5 exposure and BW-for-GA z-scores, restricted to participants with fully observed 37-week exposure series, by SREI group."
  ) %>%
  autofit()

doc_cc <- read_docx() %>%
  body_add_flextable(ft_cc)

print(doc_cc, target = file.path(dir_tables, "Supplementary_Table3_CompleteCases.docx"))

########################################
## Sensitivity analysis 4: selection bias comparison + IPW
########################################

## -- Build included vs. excluded comparison dataset --

race_lookup_full <- birth_data %>%
  distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(subject_id, raceethnicitycombined)

ids_all_enrolled <- unique(birth_data$subject_id)
ids_in_overlap   <- df_overlap %>% distinct(subject_id) %>% pull(subject_id)
ids_excluded_all <- setdiff(ids_all_enrolled, keep_ids)
ids_no_pm25      <- setdiff(ids_all_enrolled, ids_in_overlap)

cat("Total enrolled: ", length(ids_all_enrolled), "\n")
cat("Final analytic: ", length(keep_ids), "\n")
cat("Total excluded: ", length(ids_excluded_all), "\n")

covar_excl_pm25 <- df_overlap %>%
  filter(subject_id %in% ids_excluded_all) %>%
  group_by(subject_id) %>%
  summarise(
    maternalage     = first_nonNA(maternalage),
    prepregnancybmi = first_nonNA(prepregnancybmi),
    parity_cat      = first_nonNA(parity_cat),
    insurance_cat   = first_nonNA(insurance_cat),
    maritalstatus   = first_nonNA(maritalstatus),
    season          = first_nonNA(season),
    year            = first_nonNA(year_fac),
    .groups = "drop"
  ) %>%
  mutate(had_pm25 = TRUE)

covar_excl_nopm25 <- birth_data %>%
  filter(subject_id %in% ids_no_pm25) %>%
  distinct(subject_id, .keep_all = TRUE) %>%
  mutate(
    maternalage     = as.numeric(maternalage),
    prepregnancybmi = as.numeric(prepregnancybmi),
    year = factor(case_when(
      birthdate >= ymd("2020-01-01") & birthdate <= ymd("2020-12-31") ~ "2020",
      birthdate >= ymd("2021-01-01") & birthdate <= ymd("2021-12-31") ~ "2021",
      birthdate >= ymd("2022-01-01") & birthdate <= ymd("2022-12-31") ~ "2022",
      TRUE ~ NA_character_
    ), levels = c("2020", "2021", "2022")),
    season   = factor(
      season_from_date(lmp),
      levels = c("Summer", "Spring", "Autumn", "Winter")
    ),
    had_pm25 = FALSE
  ) %>%
  dplyr::select(
    subject_id, maternalage, prepregnancybmi, parity_cat,
    insurance_cat, maritalstatus, season, year, had_pm25
  )

covar_excluded_all <- bind_rows(covar_excl_pm25, covar_excl_nopm25) %>%
  left_join(
    birth_data %>%
      distinct(subject_id, .keep_all = TRUE) %>%
      dplyr::select(subject_id, birthweight_grams, gestationalagedays),
    by = "subject_id"
  ) %>%
  left_join(srei_raw %>% dplyr::select(subject_id, srei), by = "subject_id") %>%
  left_join(race_lookup_full, by = "subject_id") %>%
  mutate(group = "Excluded")

covar_included_all <- covar %>%
  filter(subject_id %in% keep_ids) %>%
  left_join(
    birth_data %>%
      distinct(subject_id, .keep_all = TRUE) %>%
      dplyr::select(subject_id, birthweight_grams, gestationalagedays),
    by = "subject_id"
  ) %>%
  left_join(race_lookup_full, by = "subject_id") %>%
  mutate(
    raceethnicitycombined = coalesce(raceethnicitycombined.x, raceethnicitycombined.y),
    maternalage     = as.numeric(maternalage),
    prepregnancybmi = as.numeric(prepregnancybmi),
    group           = "Included"
  ) %>%
  dplyr::select(-raceethnicitycombined.x, -raceethnicitycombined.y)

recode_vars <- function(df) {
  df %>%
    mutate(
      race_label = case_when(
        as.character(raceethnicitycombined) == "2" ~ "Asian",
        as.character(raceethnicitycombined) == "3" ~ "Black",
        as.character(raceethnicitycombined) == "4" ~ "Hispanic",
        as.character(raceethnicitycombined) == "6" ~ "White",
        as.character(raceethnicitycombined) == "1" ~ "Unknown",
        as.character(raceethnicitycombined) %in% c("5", "7", "8") ~ "Other",
        TRUE ~ "Unknown"
      ) %>%
        factor(levels = c("Asian", "Black", "Hispanic", "White", "Other", "Unknown")),
      parity_cat = factor(
        parity_cat, levels = c("0", "1"),
        labels = c("Nulliparous", "Multiparous")
      ),
      insurance_cat = factor(
        as.character(insurance_cat),
        levels = c("1", "2", "3"),
        labels = c("Private", "Public", "Self-pay")
      ) %>%
        fct_collapse(Private = c("Private", "Self-pay"), Public = "Public"),
      maritalstatus = factor(
        as.character(maritalstatus),
        levels = as.character(1:7),
        labels = c(
          "Divorced", "Legally Separated/Separated",
          "Married/Civil Union", "Significant Other/Life Partner",
          "Single", "Unknown/Other", "Widowed"
        )
      ) %>%
        fct_collapse(
          `Legally Separated` = c("Divorced", "Legally Separated/Separated"),
          `Significant Other` = c("Married/Civil Union", "Significant Other/Life Partner"),
          Other               = c("Widowed", "Unknown/Other")
        ) %>%
        fct_relevel("Significant Other"),
      season = factor(season, levels = c("Summer", "Spring", "Autumn", "Winter")),
      year   = factor(as.character(year), levels = c("2020", "2021", "2022"))
    )
}

compare_vars <- c(
  "maternalage", "prepregnancybmi", "parity_cat", "insurance_cat",
  "maritalstatus", "race_label", "season", "year", "srei",
  "gestationalagedays", "birthweight_grams"
)
cat_vars <- c(
  "parity_cat", "insurance_cat", "maritalstatus", "race_label", "season", "year"
)
cont_nonnormal <- c(
  "maternalage", "prepregnancybmi", "srei", "gestationalagedays", "birthweight_grams"
)

combined_all <- bind_rows(covar_included_all, covar_excluded_all) %>%
  recode_vars() %>%
  dplyr::select(subject_id, group, any_of(compare_vars))

## -- Table + SMDs --

tab_selection <- CreateTableOne(
  vars = compare_vars, strata = "group", data = combined_all,
  factorVars = cat_vars, addOverall = FALSE
)

compute_smd <- function(data, var, is_cat = FALSE) {
  inc <- data %>% filter(group == "Included") %>% pull(!!sym(var))
  exc <- data %>% filter(group == "Excluded") %>% pull(!!sym(var))
  if (is_cat) {
    t1 <- prop.table(table(inc)); t2 <- prop.table(table(exc))
    all_levels <- union(names(t1), names(t2))
    p1 <- as.numeric(t1[all_levels]); p1[is.na(p1)] <- 0
    p2 <- as.numeric(t2[all_levels]); p2[is.na(p2)] <- 0
    denom <- (p1 * (1 - p1) + p2 * (1 - p2)) / 2
    denom[denom == 0] <- NA
    sqrt(sum((p1 - p2)^2 / denom, na.rm = TRUE))
  } else {
    m1 <- mean(inc, na.rm = TRUE); s1 <- sd(inc, na.rm = TRUE)
    m2 <- mean(exc, na.rm = TRUE); s2 <- sd(exc, na.rm = TRUE)
    abs(m1 - m2) / sqrt((s1^2 + s2^2) / 2)
  }
}

smd_results <- data.frame(
  Variable = compare_vars,
  SMD = round(sapply(compare_vars, function(v) {
    compute_smd(combined_all, v, is_cat = v %in% cat_vars)
  }), 3),
  row.names = NULL
)
print(smd_results)

tab_mat <- print(
  tab_selection, showAllLevels = TRUE, printToggle = FALSE,
  quote = FALSE, noSpaces = TRUE, nonnormal = cont_nonnormal
)

smd_vec <- setNames(smd_results$SMD, smd_results$Variable)
smd_col <- rep(NA_character_, nrow(tab_mat))
names(smd_col) <- rownames(tab_mat)
for (v in compare_vars) {
  match_rows <- grep(paste0("^", v), rownames(tab_mat))
  if (length(match_rows) > 0) smd_col[match_rows[1]] <- as.character(smd_vec[v])
}

write.csv(
  cbind(tab_mat, SMD = smd_col),
  file.path(dir_tables, "Supplementary_Table4_SelectionBias.csv"),
  row.names = TRUE
)

## -- IPW analysis --
## Note: year of birth is excluded from the propensity score model
## (near-complete separation: excluded participants disproportionately
## lack a recorded birth date, reflecting pregnancy loss or delivery
## outside the study hospital).

covar_eligible <- bind_rows(
  covar_included_all %>% mutate(included = 1L),
  covar_excluded_all %>% mutate(included = 0L)
) %>%
  recode_vars() %>%
  mutate(age = as.numeric(maternalage), bmi = as.numeric(prepregnancybmi)) %>%
  filter(complete.cases(dplyr::select(
    ., age, bmi, parity_cat, insurance_cat, maritalstatus, season, srei
  )))

cat("IPW model sample size (complete cases):", nrow(covar_eligible), "\n")

ipw_fit <- weightit(
  included ~ age + bmi + parity_cat + insurance_cat +
    maritalstatus + season + srei,
  data     = covar_eligible,
  method   = "glm",
  link     = "logit",
  estimand = "ATE"
)

print(summary(ipw_fit))
print(quantile(ipw_fit$weights, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm = TRUE))

cat("Weights >10:", sum(ipw_fit$weights > 10, na.rm = TRUE), "\n")
cat("Weights >5: ", sum(ipw_fit$weights > 5,  na.rm = TRUE), "\n")

weights_df <- data.frame(
  subject_id = covar_eligible$subject_id[covar_eligible$included == 1],
  ipw        = ipw_fit$weights[covar_eligible$included == 1]
)

dat_ipw <- dat %>%
  left_join(weights_df, by = "subject_id") %>%
  left_join(preg_avg %>% filter(subject_id %in% keep_ids), by = "subject_id")

n_dropped_ipw <- sum(is.na(dat_ipw$ipw))
cat("Participants excluded from IPW-weighted comparison due to missing covariates:",
    n_dropped_ipw, "\n")

dat_ipw <- dat_ipw %>% filter(!is.na(ipw))

cat("IPW-weighted comparison sample size:", nrow(dat_ipw), "\n")

lm_unweighted <- lm(
  bw_z_intergrowth ~ preg_avg_pm25 + parity + bmi + age +
    insurance + maritalstatus + season + year,
  data = dat_ipw
)
lm_ipw <- lm(
  bw_z_intergrowth ~ preg_avg_pm25 + parity + bmi + age +
    insurance + maritalstatus + season + year,
  data = dat_ipw, weights = ipw
)

format_coef <- function(mod, var) {
  b <- coef(summary(mod))[var, ]
  sprintf(
    "β = %.4f, SE = %.4f, 95%% CI: %.4f, %.4f",
    b["Estimate"], b["Std. Error"],
    b["Estimate"] - 1.96 * b["Std. Error"],
    b["Estimate"] + 1.96 * b["Std. Error"]
  )
}
cat("\n── Unweighted PM2.5 coefficient ──\n")
cat(format_coef(lm_unweighted, "preg_avg_pm25"), "\n")

cat("\n── IPW-weighted PM2.5 coefficient ──\n")
cat(format_coef(lm_ipw, "preg_avg_pm25"), "\n")

########################################
## Sensitivity analysis 5: E-value
########################################

sd_outcome <- sd(dat$bw_z_intergrowth, na.rm = TRUE)
cat("SD of BW-for-GA z-score:", sd_outcome, "\n")

fit_summary <- summary(fit_adj)

est_high <- fit_summary$cumulative[fit_summary$cumulative$group == "≥75th percentile", "mean"]
se_high  <- fit_summary$cumulative[fit_summary$cumulative$group == "≥75th percentile", "sd"]
est_low  <- fit_summary$cumulative[fit_summary$cumulative$group == "<75th percentile", "mean"]
se_low   <- fit_summary$cumulative[fit_summary$cumulative$group == "<75th percentile", "sd"]

cat("\nExtracted estimates:\n")
cat("High disadvantage: mean =", est_high, ", SD =", se_high, "\n")
cat("Low disadvantage:  mean =", est_low,  ", SD =", se_low,  "\n")

eval_high <- evalues.MD(est = est_high, se = se_high, sd = sd_outcome)
cat("\nE-value: high structural disadvantage group\n")
print(eval_high)

eval_low <- evalues.MD(est = est_low, se = se_low, sd = sd_outcome)
cat("\nE-value: lower structural disadvantage group\n")
print(eval_low)

########################################
## Table 1 + Flowchart
########################################

ids_in_overlap <- unique(df_overlap$subject_id)

race_lookup <- df_overlap %>%
  filter(subject_id %in% keep_ids) %>%
  group_by(subject_id) %>%
  summarise(
    raceethnicitycombined = first_nonNA(raceethnicitycombined),
    .groups = "drop"
  )

preg_avg_tbl <- preg_avg %>%
  filter(subject_id %in% keep_ids)

tbl_dat <- dat %>%
  dplyr::select(
    subject_id, srei1, srei,
    bmi, age, gestationalagedays,
    birthweight_grams,
    sex, parity, insurance,
    maritalstatus, season, year
  ) %>%
  left_join(preg_avg_tbl, by = "subject_id") %>%
  left_join(race_lookup, by = "subject_id") %>%
  mutate(
    race_label = case_when(
      raceethnicitycombined == 2 ~ "Asian",
      raceethnicitycombined == 3 ~ "Black",
      raceethnicitycombined == 4 ~ "Hispanic",
      raceethnicitycombined == 6 ~ "White",
      raceethnicitycombined == 1 ~ "Unknown",
      raceethnicitycombined %in% c(5, 7, 8) ~ "Other",
      TRUE ~ "Unknown"
    ) %>%
      factor(levels = c("Asian", "Black", "Hispanic", "White", "Other", "Unknown")),
    srei_group = srei1
  )

N_total <- nrow(tbl_dat)
N_high  <- sum(tbl_dat$srei_group == "≥75th percentile", na.rm = TRUE)
N_low   <- sum(tbl_dat$srei_group == "<75th percentile", na.rm = TRUE)

cat(sprintf("N total: %s | N high-SREI: %s | N low-SREI: %s\n",
            comma(N_total), comma(N_high), comma(N_low)))

hi <- function(x) x[tbl_dat$srei_group == "≥75th percentile"]
lo <- function(x) x[tbl_dat$srei_group == "<75th percentile"]

add_cat <- function(rows, key, label, var, p_val) {
  rows[[paste0(key, "_hdr")]] <- make_row(
    label, "", "", "", fmt_p(p_val), is_header = TRUE
  )
  for (lv in levels(tbl_dat[[var]])) {
    rows[[paste0(key, "_", lv)]] <- make_row(
      paste0("  ", lv),
      n_pct(tbl_dat[[var]] == lv),
      n_pct(hi(tbl_dat[[var]]) == lv),
      n_pct(lo(tbl_dat[[var]]) == lv),
      "",
      indent = TRUE
    )
  }
  rows
}

########################################
## Build Table 1 rows
########################################

rows <- list()

p_pm25 <- wilcox_p("preg_avg_pm25",      "srei_group", tbl_dat)
p_bwt  <- wilcox_p("birthweight_grams",  "srei_group", tbl_dat)
p_age  <- wilcox_p("age",                "srei_group", tbl_dat)
p_bmi  <- wilcox_p("bmi",                "srei_group", tbl_dat)
p_ga   <- wilcox_p("gestationalagedays", "srei_group", tbl_dat)

rows[["pm25"]] <- make_row(
  "Pregnancy-average PM2.5, ug/m3",
  med_iqr(tbl_dat$preg_avg_pm25),
  med_iqr(hi(tbl_dat$preg_avg_pm25)),
  med_iqr(lo(tbl_dat$preg_avg_pm25)),
  fmt_p(p_pm25)
)

p_srei <- wilcox_p("srei", "srei_group", tbl_dat)

rows[["srei_score"]] <- make_row(
  "SREI score",
  med_iqr(tbl_dat$srei),
  med_iqr(hi(tbl_dat$srei)),
  med_iqr(lo(tbl_dat$srei)),
  fmt_p(p_srei)
)

rows[["bwt"]] <- make_row(
  "Birthweight, grams",
  med_iqr(tbl_dat$birthweight_grams),
  med_iqr(hi(tbl_dat$birthweight_grams)),
  med_iqr(lo(tbl_dat$birthweight_grams)),
  fmt_p(p_bwt)
)

rows[["age"]] <- make_row(
  "Maternal age, years",
  med_iqr(tbl_dat$age),
  med_iqr(hi(tbl_dat$age)),
  med_iqr(lo(tbl_dat$age)),
  fmt_p(p_age)
)

rows[["bmi"]] <- make_row(
  "Pre-pregnancy BMI, kg/m2",
  med_iqr(tbl_dat$bmi),
  med_iqr(hi(tbl_dat$bmi)),
  med_iqr(lo(tbl_dat$bmi)),
  fmt_p(p_bmi)
)

rows[["ga"]] <- make_row(
  "Gestational age, days",
  med_iqr(tbl_dat$gestationalagedays),
  med_iqr(hi(tbl_dat$gestationalagedays)),
  med_iqr(lo(tbl_dat$gestationalagedays)),
  fmt_p(p_ga)
)

rows <- add_cat(rows, "sex",    "Sex",
                "sex",           chisq_p("sex",           "srei_group", tbl_dat))
rows <- add_cat(rows, "parity", "Parity",
                "parity",        chisq_p("parity",        "srei_group", tbl_dat))
rows <- add_cat(rows, "race",   "Race/Ethnicity",
                "race_label",    chisq_p("race_label",    "srei_group", tbl_dat))
rows <- add_cat(rows, "ins",    "Insurance",
                "insurance",     chisq_p("insurance",     "srei_group", tbl_dat))
rows <- add_cat(rows, "mar",    "Marital status",
                "maritalstatus", chisq_p("maritalstatus", "srei_group", tbl_dat))
rows <- add_cat(rows, "season", "Season of gestation",
                "season",        chisq_p("season",        "srei_group", tbl_dat))
rows <- add_cat(rows, "year",   "Year of birth",
                "year",          chisq_p("year",          "srei_group", tbl_dat))

tbl1 <- bind_rows(rows) %>%
  dplyr::select(Characteristic, Overall, Low_SREI, High_SREI, p_value, is_header, indent)

is_header_rows <- which(tbl1$is_header == TRUE)
indent_rows    <- which(tbl1$indent == TRUE)

tbl1 <- tbl1 %>% dplyr::select(-is_header, -indent)

colnames(tbl1) <- c(
  "Characteristic",
  sprintf("Overall\n(N = %s)", comma(N_total)),
  sprintf("Lower structural disadvantage\n(SREI <75th percentile)\n(n = %s)", comma(N_low)),
  sprintf("High structural disadvantage\n(SREI ≥75th percentile)\n(n = %s)", comma(N_high)),
  "p-value"
)

########################################
## Format and save Table 1
########################################

gt_tbl <- gt::gt(tbl1) %>%
  gt::tab_style(
    style     = gt::cell_text(weight = "bold"),
    locations = gt::cells_body(rows = is_header_rows)
  ) %>%
  gt::tab_style(
    style     = gt::cell_text(indent = gt::px(16)),
    locations = gt::cells_body(columns = "Characteristic", rows = indent_rows)
  ) %>%
  gt::cols_align(align = "left", columns = "Characteristic") %>%
  gt::cols_align(align = "center", columns = 2:5) %>%
  gt::tab_options(
    table.font.name                   = "Times New Roman",
    table.font.size                   = gt::px(10),
    data_row.padding                  = gt::px(3),
    table.border.top.style            = "solid",
    table.border.top.width            = gt::px(2),
    table.border.top.color            = "black",
    table.border.bottom.style         = "solid",
    table.border.bottom.width         = gt::px(2),
    table.border.bottom.color         = "black",
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = gt::px(1),
    column_labels.border.bottom.color = "black",
    column_labels.font.weight         = "bold"
  ) %>%
  gt::tab_header(
    title = gt::md(sprintf(
      "**Table 1.** Cohort Characteristics by Neighborhood Structural Disadvantage (N = %s)",
      comma(N_total)
    ))
  ) %>%
  gt::tab_source_note(
    source_note = paste0(
      "Values are median (IQR) or n (%). ",
      "p-values from Wilcoxon rank-sum test (continuous) or ",
      "chi-square test (categorical). ",
      "SREI = Structural Racism Effect Index."
    )
  )

print(gt_tbl)

ft_tbl1 <- flextable::flextable(tbl1) %>%
  flextable::bold(i = is_header_rows, bold = TRUE) %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::align(align = "left", part = "all") %>%
  flextable::align(j = 2:5, align = "center", part = "all") %>%
  flextable::padding(i = indent_rows, j = 1, padding.left = 20, part = "body") %>%
  flextable::border_remove() %>%
  flextable::hline_top(
    part   = "header",
    border = officer::fp_border(color = "black", width = 1.5)
  ) %>%
  flextable::hline_bottom(
    part   = "header",
    border = officer::fp_border(color = "black", width = 1)
  ) %>%
  flextable::hline_bottom(
    part   = "body",
    border = officer::fp_border(color = "black", width = 1.5)
  ) %>%
  flextable::border_inner_h(
    part   = "body",
    border = officer::fp_border(color = "#CCCCCC", width = 0.5)
  ) %>%
  flextable::bold(part = "header", bold = TRUE) %>%
  flextable::add_footer_lines(
    "Values are presented as median (IQR) or n (%). p-values from Wilcoxon rank-sum test for continuous variables and chi-square test for categorical variables. SREI = Structural Racism Effect Index."
  ) %>%
  flextable::fontsize(i = 1, size = 8, part = "footer") %>%
  flextable::font(fontname = "Times New Roman", part = "footer") %>%
  flextable::align(align = "left", part = "footer") %>%
  flextable::set_table_properties(width = 1, layout = "autofit") %>%
  flextable::set_caption(
    "Table 1. Cohort Characteristics by Neighborhood Structural Disadvantage"
  )

doc_tbl1 <- officer::read_docx() %>%
  flextable::body_add_flextable(ft_tbl1)

print(
  doc_tbl1,
  target = file.path(dir_tables, "Table1_demographics.docx")
)

########################################
## Exclusion flowchart
########################################

n_enrolled  <- n_distinct(birth_data$subject_id)
n_pm25      <- n_distinct(df_overlap$subject_id)
n_no_pm25   <- n_enrolled - n_pm25
n_final     <- length(keep_ids)
n_excluded  <- n_pm25 - n_final

eligible   <- birth_data %>% filter(subject_id %in% ids_in_overlap)
n_preterm  <- sum(eligible$gestationalagedays < 259, na.rm = TRUE)
n_miss_bwt <- sum(is.na(eligible$birthweight_grams))
n_twins    <- sum(as.character(eligible$twin) == "1", na.rm = TRUE)
n_sparse   <- sum(week_counts$n_weeks_nonNA < 35)

cat("\nFlowchart counts:\n")
cat(sprintf("  Enrolled:            %s\n", comma(n_enrolled)))
cat(sprintf("  No PM2.5 coverage:   %s\n", comma(n_no_pm25)))
cat(sprintf("  Had PM2.5 data:      %s\n", comma(n_pm25)))
cat(sprintf("  Preterm (<37 wks):   %s\n", comma(n_preterm)))
cat(sprintf("  Missing birthweight: %s\n", comma(n_miss_bwt)))
cat(sprintf("  Confirmed twin:      %s\n", comma(n_twins)))
cat(sprintf("  <35 wks PM2.5:       %s\n", comma(n_sparse)))
cat(sprintf("  Total excluded:      %s\n", comma(n_excluded)))
cat(sprintf("  Final sample:        %s\n", comma(n_final)))

mx  <- 0.30
bx  <- 0.18
by  <- 0.07
ex  <- 0.74
ebx <- 0.22
y1  <- 0.88
y2  <- 0.62
y3  <- 0.36
ym1 <- (y1 - by + y2 + by) / 2
ym2 <- (y2 - by + y3 + by) / 2

flowchart <- ggplot() +
  annotate(
    "rect",
    xmin = mx - bx, xmax = mx + bx,
    ymin = y1 - by, ymax = y1 + by,
    fill = "white", color = "black", linewidth = 0.8
  ) +
  annotate(
    "text", x = mx, y = y1,
    label = sprintf("Generation C enrolled\n(N = %s)", comma(n_enrolled)),
    size = 3.2, hjust = 0.5, vjust = 0.5
  ) +
  annotate(
    "segment",
    x = mx, xend = mx, y = y1 - by, yend = y2 + by,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = mx, xend = ex - ebx, y = ym1, yend = ym1,
    linewidth = 0.5
  ) +
  annotate(
    "rect",
    xmin = ex - ebx, xmax = ex + ebx,
    ymin = ym1 - 0.075, ymax = ym1 + 0.075,
    fill = "white", color = "black", linewidth = 0.6
  ) +
  annotate(
    "text", x = ex, y = ym1,
    label = sprintf("No geocoded address\nor PM2.5 coverage \n(n = %s)", comma(n_no_pm25)),
    size = 3.2, hjust = 0.5, vjust = 0.5, lineheight = 1.3
  ) +
  annotate(
    "rect",
    xmin = mx - bx, xmax = mx + bx,
    ymin = y2 - by, ymax = y2 + by,
    fill = "white", color = "black", linewidth = 0.8
  ) +
  annotate(
    "text", x = mx, y = y2,
    label = sprintf("Had PM2.5 data \n(n = %s)", comma(n_pm25)),
    size = 3.2, hjust = 0.5, vjust = 0.5
  ) +
  annotate(
    "segment",
    x = mx, xend = mx, y = y2 - by, yend = y3 + by,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = mx, xend = ex - ebx, y = ym2, yend = ym2,
    linewidth = 0.5
  ) +
  annotate(
    "rect",
    xmin = ex - ebx, xmax = ex + ebx,
    ymin = ym2 - 0.11, ymax = ym2 + 0.11,
    fill = "white", color = "black", linewidth = 0.6
  ) +
  annotate(
    "richtext", x = ex, y = ym2,
    label = sprintf(
      paste0(
        "**Excluded (n = %s total)***:<br>",
        "Preterm (<37 wks): n = %s<br>",
        "Missing birthweight: n = %s<br>",
        "Confirmed twin: n = %s<br>",
        "<35 wks PM2.5: n = %s"
      ),
      comma(n_excluded), comma(n_preterm),
      comma(n_miss_bwt), comma(n_twins), comma(n_sparse)
    ),
    size = 3.2, hjust = 0.5, vjust = 0.5, lineheight = 1.3,
    fill = NA, label.color = NA
  ) +
  annotate(
    "rect",
    xmin = mx - bx, xmax = mx + bx,
    ymin = y3 - by, ymax = y3 + by,
    fill = "white", color = "black", linewidth = 0.8
  ) +
  annotate(
    "text", x = mx, y = y3,
    label = sprintf("Final analytic sample\n(n = %s)", comma(n_final)),
    size = 3.2, hjust = 0.5, vjust = 0.5
  ) +
  annotate(
    "text", x = 0.5, y = 0.12,
    label = "*Exclusion categories are not mutually exclusive",
    size = 2.7, hjust = 0.5, color = "grey30"
  ) +
  xlim(0.05, 1.00) +
  ylim(0.04, 1.00) +
  theme_void() +
  theme(
    plot.title  = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(title = "Figure 1. Participant Flowchart")

print(flowchart)

ggsave(
  filename = file.path(dir_figures, "Flowchart_exclusions.pdf"),
  plot = flowchart,
  width = 7,
  height = 8
)

ggsave(
  filename = file.path(dir_figures, "Flowchart_exclusions.png"),
  plot = flowchart,
  width = 7,
  height = 8,
  dpi = 300
)

########################################
## End
########################################

sessionInfo()
