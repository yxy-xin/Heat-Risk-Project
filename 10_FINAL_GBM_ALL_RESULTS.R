# ============================================================
# 10_FINAL_GBM_ALL_RESULTS.R
# Identical to 09 in all model specs, parameters, SHAP, clustering
# and manual-review export.
#
# Added outputs for heat_injury_analysis_details.docx and poster:
#   [DOC-1]  confirmed_heat_yearly_by_state.xlsx
#   [DOC-2]  all_claims_yearly_by_state.xlsx
#   [DOC-3]  variable_oneway_summary.xlsx  (cat counts + numeric percentiles)
#   [DOC-4]  predictor_list.xlsx
#   [DOC-5]  gbm_config.txt
#   [DOC-6]  feature_importance_all.xlsx
#   [DOC-7]  gains_chart_train_test.png
#   [DOC-8]  avse_pred_overall.png
#   [DOC-9]  avse_top20_<feature>.png  (x20, also bundled in all_shap_plots.pdf)
#   [DOC-10] dist_prediction_table.xlsx  (Y=0 / Y=1 counts by score band)
#   [DOC-11] underreported_by_year.xlsx
#   [DOC-12] underreported_by_state.xlsx
#
# Everything else (model, SHAP, score dist, clustering, manual review)
# is 09 code reproduced verbatim — nothing changed.
# ============================================================


# -------------------------------------------------------
# 0. Packages
# -------------------------------------------------------
required_packages <- c(
  "dplyr", "tidyr", "tidyverse", "purrr", "stringr", "lubridate",
  "slider", "readr", "xgboost", "Matrix", "caret", "janitor",
  "PRROC", "writexl", "naniar", "ggplot2"
)

missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

invisible(lapply(required_packages, library, character.only = TRUE))

options(scipen = 999)


# -------------------------------------------------------
# 1. Paths
# -------------------------------------------------------
root_dir <- "C:/Users/User/OneDrive - UNSW/Honours"
bom_dir  <- file.path(root_dir, "6. Processing BoM Data")
proj_dir <- file.path(root_dir, "3. Project 1")

years <- 1984:2025


# -------------------------------------------------------
# 2. Helper: read one year of BoM files (from 08)
# -------------------------------------------------------
read_one_year_climate <- function(y, bom_dir) {
  
  message("Reading year: ", y)
  
  f_tmax <- file.path(bom_dir, paste0("poa_daily_weighted_max_temp_", y, ".rds"))
  f_tmin <- file.path(bom_dir, paste0("poa_daily_weighted_min_temp_", y, ".rds"))
  f_vp9  <- file.path(bom_dir, paste0("poa_daily_vp9am_", y, ".rds"))
  f_vp3  <- file.path(bom_dir, paste0("poa_daily_vp_3pm_", y, ".rds"))
  
  stopifnot(file.exists(f_tmax), file.exists(f_tmin),
            file.exists(f_vp9),  file.exists(f_vp3))
  
  tmax_raw <- readRDS(f_tmax)
  tmin_raw <- readRDS(f_tmin)
  vp3_raw  <- readRDS(f_vp3)
  vp9_raw  <- readRDS(f_vp9)
  
  tmax_col <- intersect(c("tmax_weighted", "temp_max_weighted"), names(tmax_raw))
  if (length(tmax_col) == 0) stop(paste0("No recognised max temp column in year ", y))
  tmax <- tmax_raw |>
    dplyr::select(ID, Date, POA_CODE21, dplyr::all_of(tmax_col[1])) |>
    dplyr::rename(temp_max_weighted = dplyr::all_of(tmax_col[1]))
  
  tmin_col <- intersect(c("tmin_weighted", "temp_min_weighted"), names(tmin_raw))
  if (length(tmin_col) == 0) stop(paste0("No recognised min temp column in year ", y))
  tmin <- tmin_raw |>
    dplyr::select(ID, Date, POA_CODE21, dplyr::all_of(tmin_col[1])) |>
    dplyr::rename(temp_min_weighted = dplyr::all_of(tmin_col[1]))
  
  vp3_col <- intersect(c("vp15_weighted", "vp_3pm_weighted"), names(vp3_raw))
  if (length(vp3_col) == 0) stop(paste0("No recognised 3pm vapour pressure column in year ", y))
  vp3 <- vp3_raw |>
    dplyr::select(ID, Date, POA_CODE21, dplyr::all_of(vp3_col[1])) |>
    dplyr::rename(vp_3pm_weighted = dplyr::all_of(vp3_col[1]))
  
  vp9_col <- intersect(c("vp09_weighted", "vp_9am_weighted"), names(vp9_raw))
  if (length(vp9_col) == 0) stop(paste0("No recognised 9am vapour pressure column in year ", y))
  if (!"POA_CODE21" %in% names(vp9_raw)) {
    vp9_raw <- vp9_raw |>
      dplyr::left_join(tmax |> dplyr::select(ID, Date, POA_CODE21), by = c("ID", "Date"))
  }
  vp9 <- vp9_raw |>
    dplyr::select(ID, Date, POA_CODE21, dplyr::all_of(vp9_col[1])) |>
    dplyr::rename(vp_9am_weighted = dplyr::all_of(vp9_col[1]))
  
  out <- tmax |>
    dplyr::left_join(tmin, by = c("ID", "Date", "POA_CODE21")) |>
    dplyr::left_join(vp3,  by = c("ID", "Date", "POA_CODE21")) |>
    dplyr::left_join(vp9,  by = c("ID", "Date", "POA_CODE21")) |>
    dplyr::mutate(
      year  = lubridate::year(Date),
      month = lubridate::month(Date),
      doy   = lubridate::yday(Date)
    )
  
  return(out)
}


# -------------------------------------------------------
# 3. Build full 1984-2025 climate history (from 08)
# -------------------------------------------------------
clim_hist <- purrr::map_dfr(years, read_one_year_climate, bom_dir = bom_dir)

print(dim(clim_hist))
print(names(clim_hist))
print(summary(clim_hist[, c("temp_max_weighted", "temp_min_weighted", "vp_3pm_weighted")]))


# -------------------------------------------------------
# 4. Compute simplified WBGT (from 08)
#    WBGT = 0.567 * Ta + 0.393 * e + 3.94
# -------------------------------------------------------
clim_hist <- clim_hist |>
  dplyr::mutate(
    wbgt_clim = 0.567 * temp_max_weighted + 0.393 * vp_3pm_weighted + 3.94
  )
# Named wbgt_clim to avoid conflict with wbgt already in the claims dataset


# -------------------------------------------------------
# 5. Define POA-specific TMRED (from 08)
# -------------------------------------------------------
tmred_poa <- clim_hist |>
  dplyr::group_by(POA_CODE21) |>
  dplyr::summarise(
    tmred_wbgt_mean   = mean(wbgt_clim, na.rm = TRUE),
    tmred_wbgt_median = median(wbgt_clim, na.rm = TRUE),
    .groups = "drop"
  )

get_modal_bin_midpoint <- function(x, binwidth = 0.5) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  bins     <- floor(x / binwidth) * binwidth
  mode_bin <- as.numeric(names(sort(table(bins), decreasing = TRUE)[1]))
  mode_bin + binwidth / 2
}

tmred_modal <- clim_hist |>
  dplyr::group_by(POA_CODE21) |>
  dplyr::summarise(
    tmred_wbgt_modal = get_modal_bin_midpoint(wbgt_clim, binwidth = 0.5),
    .groups = "drop"
  )

tmred_poa <- tmred_poa |>
  dplyr::left_join(tmred_modal, by = "POA_CODE21")

clim_hist <- clim_hist |>
  dplyr::left_join(tmred_poa, by = "POA_CODE21") |>
  dplyr::mutate(tmred_wbgt = tmred_wbgt_modal)


# -------------------------------------------------------
# 6. Excess-above-TMRED and lag features (from 08)
# -------------------------------------------------------
clim_features_min <- clim_hist |>
  dplyr::arrange(POA_CODE21, Date) |>
  dplyr::group_by(POA_CODE21) |>
  dplyr::mutate(
    wbgt_excess_tmred = pmax(wbgt_clim - tmred_wbgt, 0),
    wbgt_lag1         = dplyr::lag(wbgt_clim, 1),
    wbgt_lag2         = dplyr::lag(wbgt_clim, 2)
  ) |>
  dplyr::ungroup() |>
  # Only keep the derived features for joining — wbgt is already in the claims data
  dplyr::select(POA_CODE21, Date, tmred_wbgt, wbgt_excess_tmred, wbgt_lag1, wbgt_lag2)

saveRDS(
  clim_features_min,
  file.path(proj_dir, "poa_daily_tmred_lag12_1984_2025.rds")
)


###############################################################################
###############################################################################
# CLAIMS PROCESSING  (06.2 structure)
###############################################################################
###############################################################################

# -------------------------------------------------------
# 7. Load claims (06.2 style)
# -------------------------------------------------------
claims_file       <- file.path(proj_dir, "swa_bom_WBGT_data.RDS")
swa_data_filtered <- readRDS(claims_file)

cat("Claims loaded:", nrow(swa_data_filtered), "rows\n")


# -------------------------------------------------------
# 8. Heat month flag (06.2)
# -------------------------------------------------------
swa_data_filtered <- swa_data_filtered %>%
  mutate(
    heat_month = if_else(
      month_occurence %in% c("10", "11", "12", "01", "02", "03"),
      1L, 0L
    )
  )


# -------------------------------------------------------
# 9. Recode non-heat mechanisms (06.2 + 08 type fix)
# -------------------------------------------------------
swa_data_filtered <- swa_data_filtered %>%
  mutate(
    heat_mechanism = if_else(
      mechanism_sub_major %in% c(
        "Stepping, kneeling or sitting on objects",
        "Hitting stationary objects",
        "Being hit by falling objects",
        "Being hit by a person accidentally",
        "Being trapped between stationary and moving objects",
        "Repetitive movement, low muscle loading",
        "Contact with hot objects",
        "Contact with electricity",
        "Vehicle incident",
        "Rollover"
      ),
      0L,
      as.integer(heat_mechanism)   # keep as integer to avoid string comparison later
    )
  )


# -------------------------------------------------------
# 10. Confirmed heat label (06.2)
# -------------------------------------------------------
swa_data_filtered <- swa_data_filtered %>%
  mutate(
    confirmed_heat = case_when(
      # Clear non-heat
      agency_minor == "Fire, flame and smoke" ~ 0,
      nature_major == "Burn"                  ~ 0,
      
      # Clear heat
      nature_classification == "Injuries" & nature_minor == "Heat stress/heat stroke"               ~ 1,
      nature_classification == "Injuries" & mechanism_sub_major == "Exposure to environmental heat" ~ 1,
      nature_classification == "Injuries" & agency_sub_minor == "Sun"                               ~ 1,
      nature_classification == "Injuries" & breakdown_agency_sub_minor == "Sun"                     ~ 1,
      
      TRUE ~ 0
    )
  )

cat("Confirmed heat cases:", sum(swa_data_filtered$confirmed_heat), "\n")
cat("Total claims:        ", nrow(swa_data_filtered), "\n")

# -------------------------------------------------------
# [DOC-1] Time-series confirmed heat claims by accident year × state
# [DOC-2] Time-series all claims by accident year × state
# -------------------------------------------------------
# 'accident_year' is not a column — derive it from the date of occurrence
swa_data_filtered <- swa_data_filtered %>%
  mutate(accident_year = lubridate::year(as.Date(`Date of Occurence (D1)`)))

confirmed_yearly_state <- swa_data_filtered %>%
  filter(confirmed_heat == 1) %>%
  count(accident_year, `Jurisdiction (SWA2)`, name = "n") %>%
  tidyr::pivot_wider(names_from = `Jurisdiction (SWA2)`, values_from = n, values_fill = 0) %>%
  arrange(accident_year)

all_claims_yearly_state <- swa_data_filtered %>%
  count(accident_year, `Jurisdiction (SWA2)`, name = "n") %>%
  tidyr::pivot_wider(names_from = `Jurisdiction (SWA2)`, values_from = n, values_fill = 0) %>%
  arrange(accident_year)

writexl::write_xlsx(confirmed_yearly_state, "confirmed_heat_yearly_by_state.xlsx")
writexl::write_xlsx(all_claims_yearly_state, "all_claims_yearly_by_state.xlsx")
cat("[DOC-1] Saved: confirmed_heat_yearly_by_state.xlsx\n")
cat("[DOC-2] Saved: all_claims_yearly_by_state.xlsx\n")


# -------------------------------------------------------
# 11. Not-heat label and filter (06.2)
# -------------------------------------------------------
swa_data_filtered <- swa_data_filtered %>%
  mutate(
    not_heat = case_when(
      confirmed_heat == 1                              ~ 0,
      nature_classification == "Diseases and conditions" ~ 1,
      heat_mechanism == 0                              ~ 1,
      heat_month == 0                                  ~ 1,
      agency_minor == "Fire, flame and smoke"          ~ 1,
      nature_major == "Burn"                           ~ 1,
      wbgt_percentile < 0.5                            ~ 1,
      TRUE                                             ~ 0
    )
  )

cat("Confirmed not-heat:", sum(swa_data_filtered$not_heat), "\n")
stopifnot(nrow(swa_data_filtered[swa_data_filtered$not_heat == 1 &
                                   swa_data_filtered$confirmed_heat == 1, ]) == 0)

semi_supervised <- swa_data_filtered %>%
  filter(not_heat == 0)

cat("Semi-supervised rows:", nrow(semi_supervised), "\n")


# -------------------------------------------------------
# 12. Merge TMRED/lag features (08 pipeline -> 06.2 claims)
# -------------------------------------------------------
semi_supervised <- semi_supervised %>%
  mutate(
    claim_date = as.Date(`Date of Occurence (D1)`),
    claim_poa  = stringr::str_pad(
      as.character(postcode_2021), width = 4, side = "left", pad = "0"
    )
  ) %>%
  dplyr::left_join(
    clim_features_min,
    by = c("claim_date" = "Date", "claim_poa" = "POA_CODE21")
  )

merge_check <- semi_supervised %>%
  summarise(
    n                  = n(),
    missing_excess     = sum(is.na(wbgt_excess_tmred)),
    missing_lag1       = sum(is.na(wbgt_lag1)),
    missing_lag2       = sum(is.na(wbgt_lag2))
  )
print(merge_check)


# -------------------------------------------------------
# 13. Helper columns (06.2 style)
# -------------------------------------------------------
semi_supervised <- semi_supervised %>%
  mutate(
    heat_age = if_else(`Age (SWA3)` >= 16 & `Age (SWA3)` <= 25, 1L, 0L),
    
    heat_industry = case_when(
      division == "Mining"                                     ~ 1L,
      division == "Construction"                               ~ 1L,
      division == "Agriculture, Forestry and Fishing"          ~ 1L,
      division == "Manufacturing"                              ~ 1L,
      division == "Electricity, Gas, Water and Waste Services" ~ 1L,
      division == "Transport, Postal and Warehousing"          ~ 1L,
      TRUE                                                     ~ 0L
    ),
    
    heat_occupations = case_when(
      major == "Labourers"                       ~ 1L,
      major == "Technicians and Trades Workers"  ~ 1L,
      major == "Machinery Operators and Drivers" ~ 1L,
      TRUE                                       ~ 0L
    ),
    
    heat_wbgt = if_else(wbgt_percentile > 0.9, 1L, 0L)
  )


# -------------------------------------------------------
# 14. Column selection
#     06.2 predictor set + wbgt_excess_tmred, wbgt_lag1, wbgt_lag2
# -------------------------------------------------------
columns_to_keep <- c(
  "Sex (C4)",
  "Age (SWA3)",
  "month_occurence",
  "Jurisdiction (SWA2)",
  "mechanism_sub_major",
  "nature_minor",
  "agency_sub_major",
  "subdivision",
  "sub_major",
  "wbgt",
  "wbgt_percentile",
  "wbgt_3day_exposure",
  "heat_wbgt",
  "heat_age",
  "heat_industry",
  "heat_occupations",
  "wbgt_excess_tmred",   # NEW: TMRED-anchored climatological excess
  "wbgt_lag1",           # NEW: prior-day WBGT
  "wbgt_lag2",           # NEW: two-day prior WBGT
  "confirmed_heat"
)

semi_supervised_clean <- semi_supervised %>%
  dplyr::select(all_of(columns_to_keep)) %>%
  janitor::clean_names()


# -------------------------------------------------------
# 15. Missing value treatment (06.2 + 08 for new features)
# -------------------------------------------------------
miss_var_summary(semi_supervised_clean)

age_median <- median(semi_supervised_clean$age_swa3, na.rm = TRUE)
semi_supervised_clean$age_swa3[is.na(semi_supervised_clean$age_swa3)] <- age_median

semi_supervised_clean <- semi_supervised_clean %>%
  mutate(
    heat_age          = if_else(age_swa3 >= 16 & age_swa3 <= 25, 1L, 0L),
    wbgt_excess_tmred = if_else(is.na(wbgt_excess_tmred), 0,               wbgt_excess_tmred),
    wbgt_lag1         = if_else(is.na(wbgt_lag1),         wbgt_percentile, wbgt_lag1),
    wbgt_lag2         = if_else(is.na(wbgt_lag2),         wbgt_lag1,       wbgt_lag2)
  )

# Drop rows missing core WBGT columns (claims pre-1984)
semi_supervised_clean <- semi_supervised_clean %>%
  drop_na(wbgt, wbgt_percentile, heat_wbgt, wbgt_3day_exposure,
          wbgt_excess_tmred, wbgt_lag1, wbgt_lag2)

# Remove label-adjacent nature_minor category (from 08)
semi_supervised_clean <- semi_supervised_clean %>%
  filter(!grepl("Effects of weather, exposure, air pressure", nature_minor))

miss_var_summary(semi_supervised_clean)
cat("Final modelling rows:", nrow(semi_supervised_clean), "\n")

# -------------------------------------------------------
# [DOC-3] One-way summary of all modelling variables
# -------------------------------------------------------
pred_cols <- setdiff(names(semi_supervised_clean), "confirmed_heat")

cat_cols <- pred_cols[sapply(semi_supervised_clean[pred_cols], function(x)
  is.character(x) | is.factor(x) | (is.integer(x) & length(unique(x[!is.na(x)])) <= 20)
)]
num_cols <- setdiff(pred_cols, cat_cols)

cat_summary_df <- bind_rows(lapply(cat_cols, function(col) {
  semi_supervised_clean %>%
    count(value = as.character(.data[[col]]), name = "count") %>%
    mutate(variable = col, pct = round(count / sum(count) * 100, 2)) %>%
    select(variable, value, count, pct)
}))

num_summary_df <- semi_supervised_clean %>%
  select(all_of(num_cols)) %>%
  summarise(across(everything(), list(
    mean   = ~round(mean(.x,           na.rm = TRUE), 4),
    median = ~round(median(.x,         na.rm = TRUE), 4),
    p05    = ~round(quantile(.x, 0.05, na.rm = TRUE), 4),
    q25    = ~round(quantile(.x, 0.25, na.rm = TRUE), 4),
    q75    = ~round(quantile(.x, 0.75, na.rm = TRUE), 4),
    p95    = ~round(quantile(.x, 0.95, na.rm = TRUE), 4),
    min    = ~round(min(.x,            na.rm = TRUE), 4),
    max    = ~round(max(.x,            na.rm = TRUE), 4)
  ), .names = "{.col}__{.fn}")) %>%
  tidyr::pivot_longer(everything(), names_to = c("variable", "stat"), names_sep = "__") %>%
  tidyr::pivot_wider(names_from = stat, values_from = value)

writexl::write_xlsx(
  list(Categorical = cat_summary_df, Numeric = num_summary_df),
  "variable_oneway_summary.xlsx"
)
cat("[DOC-3] Saved: variable_oneway_summary.xlsx\n")

# -------------------------------------------------------
# [DOC-4] Predictor list with type and group
# -------------------------------------------------------
predictor_list_df <- data.frame(
  variable = pred_cols,
  type     = ifelse(pred_cols %in% cat_cols, "Categorical", "Numeric"),
  group    = case_when(
    grepl("wbgt|tmred|lag",                      pred_cols) ~ "WBGT / Climate",
    grepl("heat_month|month_",                   pred_cols) ~ "Temporal",
    grepl("heat_industry|subdivision|sub_major",  pred_cols) ~ "Industry / Occupation",
    grepl("heat_occup",                          pred_cols) ~ "Industry / Occupation",
    grepl("heat_age|age|sex",                    pred_cols) ~ "Demographic",
    grepl("mechanism|nature|agency",             pred_cols) ~ "Injury coding",
    TRUE                                                    ~ "Other"
  ),
  stringsAsFactors = FALSE
)
writexl::write_xlsx(predictor_list_df, "predictor_list.xlsx")
cat("[DOC-4] Saved: predictor_list.xlsx\n")


###############################################################################
###############################################################################
# MODEL FITTING
###############################################################################
###############################################################################

required_gbm_packages <- c("xgboost", "Matrix", "caret", "janitor", "PRROC", "writexl", "naniar")
missing_gbm <- required_gbm_packages[
  !required_gbm_packages %in% installed.packages()[, "Package"]
]
if (length(missing_gbm) > 0) install.packages(missing_gbm, dependencies = TRUE)

library(xgboost); library(Matrix); library(caret)
library(PRROC);   library(writexl); library(dplyr); library(ggplot2)


# -------------------------------------------------------
# 16. Train / test split (06.2: prop = 0.4 downsampling)
# -------------------------------------------------------
set.seed(123)
train_indices <- caret::createDataPartition(
  semi_supervised_clean$confirmed_heat, p = 0.8, list = FALSE
)

train_data <- semi_supervised_clean[train_indices, ]
test_data  <- semi_supervised_clean[-train_indices, ]

# Keep all positives; downsample negatives to 40% (06.2 original proportion)
positive_train <- train_data %>% filter(confirmed_heat == 1)

set.seed(123)
negative_train <- train_data %>%
  filter(confirmed_heat == 0) %>%
  slice_sample(prop = 0.4)   # KEY: 0.4 from 06.2, not 0.02 from 08

train_downsampled <- bind_rows(positive_train, negative_train)
cat("Downsampled training rows:", nrow(train_downsampled), "\n")


# -------------------------------------------------------
# 17. Sparse matrices
# -------------------------------------------------------
train_label <- train_downsampled$confirmed_heat
test_label  <- test_data$confirmed_heat

x_train <- train_downsampled[, setdiff(names(train_downsampled), "confirmed_heat")]
x_test  <- test_data[,        setdiff(names(test_data),         "confirmed_heat")]

train_matrix <- Matrix::sparse.model.matrix(~ . - 1, data = x_train)
test_matrix  <- Matrix::sparse.model.matrix(~ . - 1, data = x_test)


# -------------------------------------------------------
# 18. Remove label-leakage columns (08 pattern approach)
# -------------------------------------------------------
remove_patterns <- c(
  "Exposure to environmental heat",
  "Heat stress/heat stroke"
)

remove_cols_train <- unique(unlist(lapply(remove_patterns, function(pat) {
  grep(pat, colnames(train_matrix), value = TRUE, fixed = TRUE)
})))

remove_cols_test <- unique(unlist(lapply(remove_patterns, function(pat) {
  grep(pat, colnames(test_matrix), value = TRUE, fixed = TRUE)
})))

if (length(remove_cols_train) > 0)
  train_matrix <- train_matrix[, !colnames(train_matrix) %in% remove_cols_train, drop = FALSE]
if (length(remove_cols_test) > 0)
  test_matrix  <- test_matrix[,  !colnames(test_matrix)  %in% remove_cols_test,  drop = FALSE]


# -------------------------------------------------------
# 19. Align test columns to train
# -------------------------------------------------------
train_cols <- colnames(train_matrix)

missing_cols <- setdiff(train_cols, colnames(test_matrix))
if (length(missing_cols) > 0) {
  zeros <- Matrix::Matrix(
    0, nrow = nrow(test_matrix), ncol = length(missing_cols),
    sparse = TRUE, dimnames = list(NULL, missing_cols)
  )
  test_matrix <- cbind(test_matrix, zeros)
}
test_matrix <- test_matrix[, train_cols, drop = FALSE]

cat("Column alignment check -- identical order?:",
    identical(colnames(train_matrix), colnames(test_matrix)), "\n")

dtrain <- xgb.DMatrix(data = train_matrix, label = train_label)
dtest  <- xgb.DMatrix(data = test_matrix,  label = test_label)


# -------------------------------------------------------
# 20. Model parameters (06.2)
# -------------------------------------------------------
scale_pos_weight <- sum(train_label == 0) / sum(train_label == 1)

params <- list(
  objective        = "binary:logistic",
  eval_metric      = "aucpr",
  eta              = 0.03,
  max_depth        = 1,
  min_child_weight = 1.0,
  colsample_bytree = 0.8,
  subsample        = 0.8,
  nthread          = 1,
  scale_pos_weight = scale_pos_weight
)

# -------------------------------------------------------
# [DOC-5] GBM configuration summary
# -------------------------------------------------------
writeLines(c(
  "=== GBM Configuration (XGBoost) ===",
  "Objective:           binary:logistic",
  "Eval metric:         AUC-PR",
  "Learning rate (eta): 0.03",
  "Max depth:           1",
  "Min child weight:    1.0",
  "Colsample bytree:    0.8",
  "Subsample:           0.8",
  paste0("scale_pos_weight:    ", round(scale_pos_weight, 4),
         "  (auto: n_neg / n_pos in downsampled train)"),
  "nrounds (final):     3000  (same as 09; CV run for diagnostic only)",
  "",
  "=== Downsampling ===",
  "Positives:  all kept",
  "Negatives:  prop = 0.4  (06.2 baseline)",
  "",
  "=== Train / Test split ===",
  paste0("Train rows (pre-downsample):  ", nrow(train_data)),
  paste0("Train rows (post-downsample): ", nrow(train_downsampled)),
  paste0("  - Positives: ", sum(train_label == 1)),
  paste0("  - Negatives: ", sum(train_label == 0)),
  paste0("Test rows:                    ", nrow(test_data)),
  paste0("  - Positives: ", sum(test_label == 1)),
  paste0("  - Negatives: ", sum(test_label == 0))
), "gbm_config.txt")
cat("[DOC-5] Saved: gbm_config.txt\n")


# -------------------------------------------------------
# 21. Cross-validation (retained from 08 for diagnostics)
# -------------------------------------------------------
set.seed(123)
cv <- xgb.cv(
  params                = params,
  data                  = dtrain,
  nrounds               = 3000,
  nfold                 = 5,
  early_stopping_rounds = 50,
  maximize              = TRUE,
  verbose               = 1
)

cv_log <- cv$evaluation_log

best_round_by_test   <- cv_log$iter[which.max(cv_log$test_aucpr_mean)]
best_round_by_smooth <- {
  smoothed <- stats::loess(test_aucpr_mean ~ iter, data = cv_log, span = 0.1)
  cv_log$iter[which.max(predict(smoothed))]
}

cat("Best round (raw test AUC-PR peak):     ", best_round_by_test, "\n")
cat("Best round (smoothed test AUC-PR peak):", best_round_by_smooth, "\n")

# CV learning curve
p_cv <- ggplot(cv_log, aes(x = iter)) +
  geom_ribbon(aes(ymin = test_aucpr_mean - test_aucpr_std,
                  ymax = test_aucpr_mean + test_aucpr_std),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = train_aucpr_mean, color = "Train"),   linewidth = 0.8) +
  geom_line(aes(y = test_aucpr_mean,  color = "CV Test"), linewidth = 0.8) +
  geom_vline(xintercept = 402, linetype = "dashed", color = "darkred", linewidth = 0.7) +
  annotate("text", x = 420, y = min(cv_log$test_aucpr_mean),
           label = "nrounds = 402\n(06.2 baseline)",
           color = "darkred", size = 3.5, hjust = 0) +
  scale_color_manual(values = c("Train" = "firebrick", "CV Test" = "steelblue")) +
  labs(
    title    = "XGBoost CV learning curve -- AUC-PR",
    subtitle = "Shaded band = +/-1 SD across 5 folds | dashed = 06.2 nrounds = 402",
    x        = "Boosting round",
    y        = "AUC-PR",
    color    = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

print(p_cv)
ggsave("cv_learning_curve.png", p_cv, width = 10, height = 5.5, dpi = 300)

# Overfitting diagnostic
train_at_402 <- cv_log$train_aucpr_mean[cv_log$iter == 402]
test_at_402  <- cv_log$test_aucpr_mean[cv_log$iter  == 402]

cat("\n=== Overfitting diagnostic at nrounds = 402 ===\n")
cat("Train AUC-PR:", round(train_at_402, 4), "\n")
cat("Test  AUC-PR:", round(test_at_402,  4), "\n")
cat("Gap:         ", round(train_at_402 - test_at_402, 4), "\n")


# -------------------------------------------------------
# 22. Final model fit (06.2: nrounds = 402)
# -------------------------------------------------------
set.seed(123)
model <- xgb.train(
  params    = params,
  data      = dtrain,
  nrounds   = 3000,
  watchlist = list(train = dtrain),
  verbose   = 1
)

# -------------------------------------------------------
# 23. Test set evaluation
# -------------------------------------------------------
pred_test <- predict(model, dtest)

pr_test <- PRROC::pr.curve(
  scores.class0 = pred_test[test_label == 1],
  scores.class1 = pred_test[test_label == 0],
  curve = TRUE
)
cat("Test AUC-PR:", pr_test$auc.integral, "\n")

# Score distributions by label
plot_data <- data.frame(
  pred    = c(predict(model, dtrain), pred_test),
  dataset = c(rep("Train", nrow(train_downsampled)), rep("Test", nrow(test_data))),
  label   = c(train_label, test_label)
)
plot_data$label <- factor(plot_data$label, labels = c("Unlabelled", "Confirmed"))

ggplot(plot_data, aes(x = pred)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "steelblue", colour = "white") +
  facet_grid(label ~ dataset, scales = "free_y") +
  labs(x = "Predicted Score", y = "Count",
       title = "Distribution of Predicted Heat Scores") +
  theme_minimal()

ggsave("dist_predicted_scores.png")


# -------------------------------------------------------
# 24. Feature importance
# -------------------------------------------------------
importance_model <- xgb.importance(
  feature_names = colnames(train_matrix),
  model         = model
)

print(head(importance_model, 20))

importance_model %>%
  mutate(rank = row_number()) %>%
  filter(grepl("wbgt|agency_sub_major|heat_industry|heat_age|heat_occup", Feature)) %>%
  select(rank, Feature, Gain, Cover, Frequency) %>%
  print()

# [DOC-6] Full importance table
writexl::write_xlsx(
  importance_model %>% as.data.frame() %>% mutate(rank = row_number()),
  "feature_importance_all.xlsx"
)
cat("[DOC-6] Saved: feature_importance_all.xlsx\n")

top20_features <- importance_model$Feature[1:min(20, nrow(importance_model))]


# -------------------------------------------------------
# 25. Score all unlabelled claims (from 08)
# -------------------------------------------------------
unlabelled_claims <- semi_supervised_clean %>%
  filter(confirmed_heat == 0)

cat("Total unlabelled claims:", nrow(unlabelled_claims), "\n")

x_unlabelled      <- unlabelled_claims %>% select(-confirmed_heat)
unlabelled_matrix <- Matrix::sparse.model.matrix(~ . - 1, data = x_unlabelled)

remove_cols_unlabelled <- unique(unlist(lapply(remove_patterns, function(pat) {
  grep(pat, colnames(unlabelled_matrix), value = TRUE, fixed = TRUE)
})))
if (length(remove_cols_unlabelled) > 0) {
  unlabelled_matrix <- unlabelled_matrix[,
                                         !colnames(unlabelled_matrix) %in% remove_cols_unlabelled, drop = FALSE]
}

missing_cols <- setdiff(colnames(train_matrix), colnames(unlabelled_matrix))
if (length(missing_cols) > 0) {
  zeros <- Matrix::Matrix(0, nrow = nrow(unlabelled_matrix), ncol = length(missing_cols),
                          sparse = TRUE, dimnames = list(NULL, missing_cols))
  unlabelled_matrix <- cbind(unlabelled_matrix, zeros)
}
extra_cols <- setdiff(colnames(unlabelled_matrix), colnames(train_matrix))
if (length(extra_cols) > 0) {
  unlabelled_matrix <- unlabelled_matrix[,
                                         !colnames(unlabelled_matrix) %in% extra_cols, drop = FALSE]
}
unlabelled_matrix <- unlabelled_matrix[, colnames(train_matrix), drop = FALSE]

cat("Unlabelled matrix aligned?:",
    identical(colnames(unlabelled_matrix), colnames(train_matrix)), "\n")

dunlabelled <- xgboost::xgb.DMatrix(data = unlabelled_matrix)
unlabelled_claims$heat_prob <- predict(model, dunlabelled)

unlabelled_claims <- unlabelled_claims %>%
  mutate(
    prob_bin = cut(
      heat_prob,
      breaks         = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
      include.lowest = TRUE,
      right          = TRUE
    )
  )

bin_summary <- unlabelled_claims %>%
  count(prob_bin, name = "unlabelled_claims") %>%
  arrange(prob_bin)
print(bin_summary)

thresholds <- c(0.001, 0.005, 0.01, 0.014, 0.02, 0.05, 0.1)
cum_summary <- data.frame(
  threshold       = thresholds,
  n_at_or_below   = sapply(thresholds, function(x) sum(unlabelled_claims$heat_prob <= x)),
  pct_at_or_below = sapply(thresholds, function(x) mean(unlabelled_claims$heat_prob <= x) * 100)
)
print(cum_summary)

# Export all unlabelled claims
write_xlsx(
  unlabelled_claims %>% filter(heat_prob > 0.5) %>% arrange(desc(heat_prob)),
  "high_score_unlabelled_claims.xlsx"
)

# -------------------------------------------------------
# [DOC-7] Gains chart (train vs test)
# -------------------------------------------------------
pred_train <- predict(model, dtrain)

make_gains <- function(pred, label, n_bins = 10) {
  data.frame(pred = pred, label = label) %>%
    arrange(desc(pred)) %>%
    mutate(cum_pos = cumsum(label), tot_pos = sum(label), cum_pct = row_number() / n()) %>%
    mutate(bin = ntile(row_number(), n_bins)) %>%
    group_by(bin) %>%
    summarise(pct_population = max(cum_pct), gain = max(cum_pos) / max(tot_pos), .groups = "drop")
}

p_gains <- bind_rows(
  make_gains(pred_train, train_label) %>% mutate(split = "Train"),
  make_gains(pred_test,  test_label)  %>% mutate(split = "Test")
) %>%
  ggplot(aes(x = pct_population, y = gain, color = split, linetype = split)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c("Train" = "firebrick", "Test" = "steelblue")) +
  labs(
    title    = "Gains Chart -- Train vs Test",
    subtitle = "% of confirmed heat cases captured vs. % of claims reviewed (sorted by descending score)",
    x = "% of Claims Reviewed", y = "% of Confirmed Heat Cases Captured",
    color = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("gains_chart_train_test.png", p_gains, width = 9, height = 6, dpi = 300)
cat("[DOC-7] Saved: gains_chart_train_test.png\n")

# -------------------------------------------------------
# [DOC-8] AvsE overall -- prediction percentile bins
# -------------------------------------------------------
test_pred_df        <- as.data.frame(as.matrix(test_matrix))
test_pred_df$pred   <- pred_test
test_pred_df$actual <- test_label

test_pred_df_corrected <- test_pred_df
test_pred_df_corrected_wbgt <- test_pred_df %>%
  filter(wbgt_percentile >= 0.5)

avse_overall_df <- data.frame(pred = pred_test, actual = test_label) %>%
  mutate(bin = ntile(pred, 100)) %>%
  group_by(bin) %>%
  summarise(avg_pred = mean(pred), avg_actual = mean(actual), n = n(), .groups = "drop")

p_avse_overall <- ggplot(avse_overall_df, aes(x = bin)) +
  geom_line(aes(y = avg_pred,   color = "Average Prediction"), linewidth = 1) +
  geom_line(aes(y = avg_actual, color = "Average Actual"),     linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Average Prediction" = "steelblue", "Average Actual" = "firebrick")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = "AvsE -- Overall: Prediction vs Actual by Prediction Percentile",
    subtitle = "X: prediction percentile bin (1-100) | Y: proportion confirmed heat",
    x = "Prediction Percentile Bin", y = "Proportion of Confirmed Cases", color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("avse_pred_overall.png", p_avse_overall, width = 9, height = 6, dpi = 300)
cat("[DOC-8] Saved: avse_pred_overall.png\n")

# -------------------------------------------------------
# [DOC-9] AvsE for top-20 predictors (REVISED)
# heat_month removed from model — no longer in top20_features
# wbgt / wbgt_percentile: restrict to valid region (>= 0.5 percentile)
#   because Y=0 pool excludes wbgt_percentile < 0.5 by design,
#   so the biased region is uninformative and misleading.
# -------------------------------------------------------
# Features affected by the wbgt_percentile < 0.5 exclusion filter
# — restrict AvsE to the unbiased region for these
wbgt_biased_feats <- c("wbgt_percentile", "wbgt", "wbgt_excess_tmred",
                       "wbgt_lag1", "wbgt_lag2", "wbgt_3day_exposure")

avse_plot_list <- list()
winter_month_feats <- paste0("month_occurence", c("04","05","06","07","08","09"))
for (feat in top20_features) {
  if (!feat %in% names(test_pred_df)) next
  
  if (feat %in% winter_month_feats) {
    message("Skipping AvsE for ", feat, " — labelling artefact (excluded from Y=0 pool)")
    next
  }
  
  # For WBGT features: restrict to wbgt_percentile >= 0.5 (valid region)
  df_feat <- if (feat %in% wbgt_biased_feats) {
    test_pred_df %>% filter(wbgt_percentile >= 0.5)
  } else {
    test_pred_df
  }
  
  vals      <- df_feat[[feat]][!is.na(df_feat[[feat]])]
  is_binary <- all(unique(vals) %in% c(0, 1))
  
  if (is_binary) {
    avse_data <- df_feat %>%
      select(bin = !!sym(feat), pred, actual) %>%
      mutate(bin = factor(as.character(bin))) %>%
      group_by(bin) %>%
      summarise(avg_pred = mean(pred), avg_actual = mean(actual), n = n(), .groups = "drop")
    
    p_avse <- ggplot(avse_data, aes(x = bin)) +
      geom_col(aes(y = avg_actual),            fill = "firebrick", alpha = 0.45, width = 0.5) +
      geom_point(aes(y = avg_pred),            color = "steelblue", size = 4) +
      geom_line(aes(y = avg_pred, group = 1),  color = "steelblue", linewidth = 0.9)
  } else {
    avse_data <- df_feat %>%
      select(feat_val = !!sym(feat), pred, actual) %>%
      filter(!is.na(feat_val)) %>%
      mutate(bin = cut(feat_val, breaks = 10, include.lowest = TRUE)) %>%
      group_by(bin) %>%
      summarise(avg_pred = mean(pred), avg_actual = mean(actual), n = n(), .groups = "drop")
    
    p_avse <- ggplot(avse_data, aes(x = bin)) +
      geom_col(aes(y = avg_actual),            fill = "firebrick", alpha = 0.45) +
      geom_point(aes(y = avg_pred),            color = "steelblue", size = 3) +
      geom_line(aes(y = avg_pred, group = 1),  color = "steelblue", linewidth = 0.9)
  }
  
  # Add note to subtitle for restricted WBGT plots
  wbgt_note <- if (feat %in% wbgt_biased_feats)
    " | Restricted to wbgt_percentile >= 0.5 (valid region — see methodology note)"
  else ""
  
  p_avse <- p_avse +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title    = paste0("AvsE: ", feat),
      subtitle = paste0("Bars = actual proportion confirmed heat | Points/line = avg prediction",
                        wbgt_note),
      x = feat, y = "Proportion Confirmed Heat"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x   = element_text(angle = 30, hjust = 1, size = 9),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8.5)
    )
  
  safe_name <- gsub("[^A-Za-z0-9_]", "_", feat)
  ggsave(paste0("avse_top20_", safe_name, ".png"), p_avse, width = 9, height = 5.5, dpi = 300)
  avse_plot_list[[feat]] <- p_avse
}
cat("[DOC-9] Saved avse_top20_<feature>.png for", length(avse_plot_list), "features\n")
# -------------------------------------------------------
# [DOC-10] Distribution of prediction table (Y=0 / Y=1 by score band)
# Built from ALL scored claims, not just the downsampled train+test split
# -------------------------------------------------------

# Confirmed positives from the test set (unaffected by downsampling)
confirmed_pred <- data.frame(
  pred   = pred_test[test_label == 1],
  actual = 1L
)

# All unlabelled claims scored by the model (the full ~514k pool)
unlabelled_pred <- data.frame(
  pred   = unlabelled_claims$heat_prob,
  actual = 0L
)

dist_pred_table <- bind_rows(confirmed_pred, unlabelled_pred) %>%
  mutate(prob_bin = cut(pred,
                        breaks         = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                        include.lowest = TRUE, right = TRUE)) %>%
  group_by(prob_bin) %>%
  summarise(
    `Obs with Y = 0 (unlabelled)` = sum(actual == 0),
    `Obs with Y = 1 (confirmed)`  = sum(actual == 1),
    .groups = "drop"
  ) %>%
  rename(`Prediction Band` = prob_bin)

writexl::write_xlsx(dist_pred_table, "dist_prediction_table.xlsx")
cat("[DOC-10] Saved: dist_prediction_table.xlsx\n")

###############################################################################
###############################################################################
# SHAP ANALYSIS  (from 08)
###############################################################################
###############################################################################

required_shap_packages <- c("shapviz", "ggplot2")
missing_shap <- required_shap_packages[
  !required_shap_packages %in% installed.packages()[, "Package"]
]
if (length(missing_shap) > 0) install.packages(missing_shap, dependencies = TRUE)

library(shapviz)
library(ggplot2)

set.seed(123)
id           <- sample(seq_len(nrow(test_matrix)), size = min(1000, nrow(test_matrix)))
X_pred_small <- as.matrix(test_matrix[id, , drop = FALSE])

sv <- shapviz(
  model,
  X_pred = X_pred_small,
  X      = as.data.frame(X_pred_small)
)

# Beeswarm -- top 15 features
sv_importance(sv, kind = "beeswarm", max_display = 15) +
  labs(title = "SHAP Beeswarm -- WBGT TMRED Base Model")
ggsave("shap_beeswarm.png", width = 10, height = 7, dpi = 300)

# Bar -- mean |SHAP|
sv_importance(sv, kind = "bar", max_display = 15) +
  labs(title = "Mean Absolute SHAP Importance -- WBGT TMRED Base Model")
ggsave("shap_bar.png", width = 9, height = 6, dpi = 300)

# Dependence -- WBGT features
sv_dependence(sv, v = "wbgt_percentile")
ggsave("shap_dep_wbgt_percentile.png", width = 7, height = 5, dpi = 300)

sv_dependence(sv, v = "wbgt_excess_tmred")
ggsave("shap_dep_wbgt_excess_tmred.png", width = 7, height = 5, dpi = 300)

sv_dependence(sv, v = "wbgt_lag1")
ggsave("shap_dep_wbgt_lag1.png", width = 7, height = 5, dpi = 300)

sv_dependence(sv, v = "wbgt_lag2")
ggsave("shap_dep_wbgt_lag2.png", width = 7, height = 5, dpi = 300)

sv_dependence(sv, v = "wbgt")
ggsave("shap_dep_wbgt_raw.png", width = 7, height = 5, dpi = 300)

# Outdoor agency dependence (if present)
agency_col <- grep("agency_sub_majorOutdoor", colnames(X_pred_small), value = TRUE)
if (length(agency_col) > 0) {
  sv_dependence(sv, v = agency_col[1])
  ggsave("shap_dep_agency_outdoor.png", width = 7, height = 5, dpi = 300)
}

# --- PDF compile -- all SHAP plots + AvsE top-20 appended at end ---
pdf("all_shap_plots.pdf", width = 12, height = 8)

grid::grid.newpage()
grid::grid.text(
  "WBGT TMRED Base Model -- SHAP Diagnostic Report",
  x = 0.5, y = 0.5,
  gp = grid::gpar(fontsize = 18, fontface = "bold")
)

print(sv_importance(sv, kind = "beeswarm", max_display = 15) +
        labs(title = "SHAP Beeswarm -- Top 15 Features"))
print(sv_importance(sv, kind = "bar", max_display = 20) +
        labs(title = "Mean Absolute SHAP Importance -- Top 20"))

dep_vars <- c("wbgt", "wbgt_percentile", "wbgt_excess_tmred",
              "wbgt_lag1", "wbgt_lag2", "wbgt_3day_exposure")
if (length(agency_col) > 0) dep_vars <- c(dep_vars, agency_col[1])
dep_vars <- dep_vars[dep_vars %in% colnames(X_pred_small)]

for (feat in dep_vars) {
  tryCatch({
    print(sv_dependence(sv, v = feat) +
            labs(title = paste0("SHAP Dependence: ", feat)))
  }, error = function(e) {
    grid::grid.newpage()
    grid::grid.text(paste0("Could not plot: ", feat, "\n", conditionMessage(e)),
                    x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 12, col = "red"))
  })
}

# Waterfall -- one confirmed heat case
confirmed_idx <- which(test_label[id] == 1)
if (length(confirmed_idx) > 0) {
  grid::grid.newpage()
  grid::grid.text("Waterfall -- Confirmed Heat Claim",
                  x = 0.5, y = 0.5,
                  gp = grid::gpar(fontsize = 16, fontface = "bold"))
  print(sv_waterfall(sv, row_id = confirmed_idx[1]) +
          labs(title = "SHAP Waterfall -- Confirmed Heat Claim"))
  print(sv_force(sv, row_id = confirmed_idx[1]) +
          labs(title = "SHAP Force -- Confirmed Heat Claim"))
}

# [DOC-9 in PDF] AvsE overall + top-20
grid::grid.newpage()
grid::grid.text("AvsE Charts -- Overall and Top-20 Predictors",
                x = 0.5, y = 0.5,
                gp = grid::gpar(fontsize = 16, fontface = "bold"))

print(p_avse_overall)

for (feat in top20_features) {
  if (feat %in% names(avse_plot_list)) {
    tryCatch(
      print(avse_plot_list[[feat]]),
      error = function(e) message("Could not print AvsE for: ", feat)
    )
  }
}

dev.off()
cat("SHAP PDF saved: all_shap_plots.pdf\n")


# --- WBGT-related feature SHAP ranks ---
shap_matrix   <- sv$S
mean_abs_shap <- colMeans(abs(shap_matrix))

imp_full <- data.frame(
  feature    = names(mean_abs_shap),
  importance = as.numeric(mean_abs_shap),
  row.names  = NULL,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(importance)) %>%
  mutate(global_rank = row_number())

cat("Total features:", nrow(imp_full), "\n")
print(head(imp_full, 5))

wbgt_ranks <- imp_full %>%
  filter(grepl("wbgt|heat_age|heat_industry|heat_occup|age_swa3",
               feature, ignore.case = TRUE)) %>%
  select(global_rank, feature, mean_abs_shap = importance) %>%
  arrange(global_rank) %>%
  mutate(
    pct_outranked = round((1 - global_rank / nrow(imp_full)) * 100, 1),
    in_top_25pct  = global_rank <= nrow(imp_full) * 0.25
  )

cat("\n=== WBGT-related feature SHAP ranks ===\n")
cat("Total features in model:", nrow(imp_full), "\n")
print(wbgt_ranks)


# -------------------------------------------------------
# [DOC-10] Prediction distribution table — by SCORE PERCENTILE bins
# and rescore all 1,495 confirmed heat cases
# -------------------------------------------------------

# --- Step 1: Rescore all confirmed positives ---
confirmed_all <- semi_supervised_clean %>% filter(confirmed_heat == 1)
x_confirmed   <- confirmed_all %>% select(-confirmed_heat)
conf_matrix   <- Matrix::sparse.model.matrix(~ . - 1, data = x_confirmed)

remove_conf <- unique(unlist(lapply(remove_patterns, function(pat)
  grep(pat, colnames(conf_matrix), value = TRUE, fixed = TRUE))))
if (length(remove_conf) > 0)
  conf_matrix <- conf_matrix[, !colnames(conf_matrix) %in% remove_conf, drop = FALSE]

miss_conf <- setdiff(colnames(train_matrix), colnames(conf_matrix))
if (length(miss_conf) > 0) {
  zeros <- Matrix::Matrix(0, nrow = nrow(conf_matrix), ncol = length(miss_conf),
                          sparse = TRUE, dimnames = list(NULL, miss_conf))
  conf_matrix <- cbind(conf_matrix, zeros)
}
conf_matrix <- conf_matrix[, colnames(train_matrix), drop = FALSE]
confirmed_scores <- predict(model, xgb.DMatrix(conf_matrix))

# --- Step 2: Compute percentile BREAKS from unlabelled scores only ---
# Percentiles are defined on the unlabelled pool (the 514k claims being ranked)
# Confirmed positives are then placed into those same bins for reference
pct_breaks <- quantile(
  unlabelled_claims$heat_prob,
  probs = seq(0, 1, by = 0.1),
  names = TRUE
)
# Force endpoints to avoid edge NAs
pct_breaks[1]   <- 0
pct_breaks[11]  <- 1

cat("\nPercentile break points (from unlabelled scores):\n")
print(round(pct_breaks, 4))

# --- Step 3: Assign percentile bins ---
unlabelled_pctbin <- cut(
  unlabelled_claims$heat_prob,
  breaks         = pct_breaks,
  include.lowest = TRUE,
  right          = TRUE,
  labels         = paste0("P", seq(10, 100, 10))   # P10 = bottom 10%, P100 = top 10%
)

confirmed_pctbin <- cut(
  confirmed_scores,
  breaks         = pct_breaks,
  include.lowest = TRUE,
  right          = TRUE,
  labels         = paste0("P", seq(10, 100, 10))
)

# --- Step 4: Build the distribution table ---
dist_pred_table <- bind_rows(
  data.frame(
    score       = unlabelled_claims$heat_prob,
    pct_bin     = unlabelled_pctbin,
    actual      = 0L
  ),
  data.frame(
    score       = confirmed_scores,
    pct_bin     = confirmed_pctbin,
    actual      = 1L
  )
) %>%
  group_by(pct_bin) %>%
  summarise(
    score_min                     = round(min(score[actual == 0], na.rm = TRUE), 4),
    score_max                     = round(max(score[actual == 0], na.rm = TRUE), 4),
    `Obs with Y = 0 (unlabelled)` = sum(actual == 0),
    `Obs with Y = 1 (confirmed)`  = sum(actual == 1),
    .groups = "drop"
  ) %>%
  rename(
    `Percentile Bin`   = pct_bin,
    `Score Range (min)` = score_min,
    `Score Range (max)` = score_max
  )

cat("\nDOC-10 distribution table:\n")
print(dist_pred_table)
cat("Y=0 total:", sum(dist_pred_table[[4]]), "\n")
cat("Y=1 total:", sum(dist_pred_table[[5]]), "  (should be 1495)\n")

writexl::write_xlsx(dist_pred_table, "dist_prediction_table.xlsx")
cat("[DOC-10] Saved: dist_prediction_table.xlsx\n")

# -------------------------------------------------------
# Fine-grained manual review: P85 to P100 in 1% increments
# -------------------------------------------------------

# --- Step 1: Compute 1% break points from P85 to P100 ---
fine_probs  <- seq(0.85, 1.00, by = 0.01)
fine_breaks <- quantile(unlabelled_claims$heat_prob, probs = fine_probs)
fine_breaks[length(fine_breaks)] <- 1   # force ceiling only

cat("\nFine percentile break points (P85-P100):\n")
print(round(fine_breaks, 6))

# Labels: P86 = 85th-86th percentile slice, ..., P100 = 99th-100th
fine_labels <- paste0("P", seq(86, 100, 1))

# --- Step 2: Assign fine bins to unlabelled claims ---
unlabelled_claims <- unlabelled_claims %>%
  mutate(
    pct_bin_fine = cut(
      heat_prob,
      breaks         = fine_breaks,
      include.lowest = TRUE,
      right          = TRUE,
      labels         = fine_labels
    )
  )

unlabelled_fine <- unlabelled_claims %>%
  filter(!is.na(pct_bin_fine))

cat("Claims in P85-P100:", nrow(unlabelled_fine), "\n")
cat("Per-bin counts:\n")
print(table(unlabelled_fine$pct_bin_fine))

# --- Step 3: Attach raw lookup columns ---
raw_unlabelled_source <- semi_supervised %>%
  filter(not_heat == 0) %>%
  drop_na(wbgt_percentile) %>%
  filter(!grepl("Effects of weather, exposure, air pressure", nature_minor)) %>%
  filter(confirmed_heat == 0)

stopifnot(nrow(raw_unlabelled_source) == nrow(unlabelled_claims))

raw_id_cols_present <- intersect(
  c("Claim Number", "Date of Occurence (D1)", "postcode_2021",
    "Jurisdiction (SWA2)", "Occupation Title", "Injury Description",
    "division", "major", "nature_minor", "mechanism_sub_major",
    "agency_sub_major", "agency_sub_minor",
    "wbgt", "wbgt_percentile", "wbgt_excess_tmred",
    "wbgt_lag1", "wbgt_lag2"),
  names(raw_unlabelled_source)
)

raw_lookup_fine <- bind_cols(
  raw_unlabelled_source %>% select(all_of(raw_id_cols_present)),
  data.frame(
    heat_prob    = unlabelled_claims$heat_prob,
    pct_bin_fine = unlabelled_claims$pct_bin_fine
  )
) %>%
  filter(!is.na(pct_bin_fine))

cat("raw_lookup_fine rows:", nrow(raw_lookup_fine), "\n")

# --- Step 4: Score reference sheet ---
# Also include confirmed heat cases per bin for reference
confirmed_pctbin_fine <- cut(
  confirmed_scores,
  breaks         = fine_breaks,
  include.lowest = TRUE,
  right          = TRUE,
  labels         = fine_labels
)

confirmed_fine_counts <- data.frame(
  pct_bin_fine  = confirmed_pctbin_fine,
  confirmed     = 1L
) %>%
  filter(!is.na(pct_bin_fine)) %>%
  count(pct_bin_fine, name = "n_confirmed_heat")

score_ref_fine <- unlabelled_fine %>%
  group_by(pct_bin_fine) %>%
  summarise(
    score_min = round(min(heat_prob), 6),
    score_max = round(max(heat_prob), 6),
    n_claims  = n(),
    .groups   = "drop"
  ) %>%
  left_join(confirmed_fine_counts, by = "pct_bin_fine") %>%
  mutate(n_confirmed_heat = if_else(is.na(n_confirmed_heat), 0L, n_confirmed_heat)) %>%
  rename(
    `Percentile Bin`      = pct_bin_fine,
    `Score Min`           = score_min,
    `Score Max`           = score_max,
    `N Unlabelled Claims` = n_claims,
    `N Confirmed Heat`    = n_confirmed_heat
  )

print(score_ref_fine)

# --- Step 5: Build per-bin review sheets ---
review_cols_fine <- c(
  "Claim Number", "Date of Occurence (D1)", "postcode_2021",
  "Jurisdiction (SWA2)", "Occupation Title", "major", "division",
  "Injury Description", "mechanism_sub_major", "nature_minor",
  "agency_sub_major", "agency_sub_minor",
  "wbgt", "wbgt_percentile", "wbgt_excess_tmred", "wbgt_lag1", "wbgt_lag2",
  "heat_prob", "pct_bin_fine",
  "Is_Heat_Related", "Reviewer_Notes"
)

prep_fine_sheet <- function(df, bin_label, n_sample = 50, seed = 42) {
  set.seed(seed)
  pool    <- df %>% filter(pct_bin_fine == bin_label)
  n_avail <- nrow(pool)
  
  # Review all claims if <=100, otherwise sample 50
  sampled <- if (n_avail <= 100) {
    pool %>% arrange(desc(heat_prob))
  } else {
    pool %>%
      slice_sample(n = n_sample, replace = FALSE) %>%
      arrange(desc(heat_prob))
  }
  
  sampled <- sampled %>%
    mutate(Is_Heat_Related = NA_character_, Reviewer_Notes = NA_character_)
  
  cat(sprintf("Bin %-4s: %5d available | %d in sheet%s\n",
              bin_label, n_avail, nrow(sampled),
              if (n_avail <= 100) " [ALL CLAIMS — small bin]" else ""))
  
  sampled %>% select(all_of(intersect(review_cols_fine, names(sampled))))
}

fine_bin_sheets <- lapply(fine_labels, function(b)
  prep_fine_sheet(raw_lookup_fine, b, n_sample = 50))
names(fine_bin_sheets) <- paste0("Bin_", fine_labels)

# --- Step 6: Confirmed heat anchors ---
set.seed(42)
confirmed_anchors <- semi_supervised %>%
  filter(confirmed_heat == 1) %>%
  slice_sample(n = 20) %>%
  mutate(
    heat_prob       = 1.0,
    pct_bin_fine    = "Confirmed (anchor)",
    Is_Heat_Related = "YES — confirmed by TOOCS coding",
    Reviewer_Notes  = "Use these as calibration anchors"
  ) %>%
  select(all_of(intersect(review_cols_fine, names(.))))

# --- Step 7: Instructions ---
instructions_fine <- data.frame(
  Section = c(
    "PURPOSE",
    "BIN MEANING",
    "SAMPLE SIZE NOTE",
    "HOW TO CALCULATE pi",
    "HOW TO CALCULATE FINAL %",
    "HOW TO CALCULATE FINAL %",
    "Is_Heat_Related", "Is_Heat_Related", "Is_Heat_Related",
    "KEY FIELDS", "KEY FIELDS", "KEY FIELDS", "KEY FIELDS", "KEY FIELDS"
  ),
  Detail = c(
    paste0("This workbook covers the top 15% of model scores (P85-P100) in 1% slices. ",
           "These are the claims the model considers most likely to be hidden heat injuries. ",
           "The Score_Reference sheet shows how many confirmed heat cases also fall in each bin."),
    paste0("P86 = claims scoring between the 85th and 86th percentile of the ~514k ",
           "unlabelled score distribution. P100 = 99th-100th percentile (highest scores). ",
           "See Score_Reference for exact score thresholds and confirmed heat counts per bin."),
    paste0("Bins with <=100 unlabelled claims show ALL claims (flagged in console output). ",
           "For larger bins, 50 are randomly sampled. Use N Unlabelled Claims from ",
           "Score_Reference as N_i in your burden calculation — not the sample size."),
    paste0("For each bin: pi_i = n_YES / (n_YES + n_NO). ",
           "Treat UNSURE consistently — either exclude from denominator or count as 0.5."),
    paste0("Estimated hidden heat = SUM over all reviewed bins of (N_i * pi_i), ",
           "where N_i = N Unlabelled Claims in Score_Reference for that bin."),
    paste0("Final % = (1495 + estimated hidden heat) / (1495 + 514822) * 100. ",
           "1495 = confirmed heat cases. 514822 = total unlabelled pool. ",
           "Report P85-P100 contribution separately from the full-range estimate."),
    "YES    — claim is clearly heat-related",
    "NO     — claim is clearly not heat-related",
    "UNSURE — cannot determine; exclude or treat as 0.5",
    "Injury Description  — heat, sun, sweating, dizziness, dehydration",
    "mechanism_sub_major — 'Exposure to environmental heat' is direct evidence",
    "nature_minor        — 'Heat stress/heat stroke' is direct evidence",
    "wbgt                — >28 is hot, >32 is very hot for outdoor workers",
    "wbgt_excess_tmred   — degrees above local long-run normal; >5 is anomalous"
  )
)

# --- Step 8: Write workbook ---
writexl::write_xlsx(
  c(
    list(
      Instructions      = instructions_fine,
      Score_Reference   = as.data.frame(score_ref_fine),
      Confirmed_Anchors = as.data.frame(confirmed_anchors)
    ),
    lapply(fine_bin_sheets, as.data.frame)
  ),
  path = "manual_review_fine_P85_P100.xlsx"
)

cat("\n=== Saved: manual_review_fine_P85_P100.xlsx ===\n")
purrr::walk2(names(fine_bin_sheets), fine_bin_sheets, function(nm, df) {
  cat(sprintf("  %-10s — %d claims in sheet\n", nm, nrow(df)))
})


###############################################################################
# VISUALISATION SECTION
# Requires from memory:
#   - pi_table               (quantile π values, P86-P100)
#   - unlabelled_claims      (with heat_prob, pct_bin_fine, jurisdiction_swa2)
#   - raw_unlabelled_source  (original rows with accident_year, Jurisdiction SWA2, postcode_2021)
#   - confirmed_heat_meta    (semi_supervised filtered to confirmed_heat == 1)
#
# Outputs:
#   [VIZ-1]  plot_underreported_by_year.png  + underreported_by_year_q.xlsx
#   [VIZ-2]  plot_underreported_by_state.png + underreported_by_state_q.xlsx
#   [VIZ-3]  leaflet postcode map            + underreported_by_poa_q.xlsx
###############################################################################

pi_table <- data.frame(
  pct_bin_fine = paste0("P", seq(86, 100)),
  pi = c(
    0.40,   # P86
    0.36,   # P87
    0.40,   # P88
    0.42,   # P89
    0.42,   # P90
    0.50,   # P91
    0.42,   # P92
    0.48,   # P93
    0.60,   # P94
    0.54,   # P95
    0.58,   # P96
    0.52,   # P97
    0.60,   # P98
    0.68,   # P99
    0.36    # P100
  ),
  stringsAsFactors = FALSE
)


library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(scales)

# -------------------------------------------------------
# STEP 1: Attach quantile π to every unlabelled claim
#         using row-position alignment (verified by stopifnot below)
# -------------------------------------------------------
stopifnot(nrow(raw_unlabelled_source) == nrow(unlabelled_claims))

raw_unlabelled_with_pi <- bind_cols(
  raw_unlabelled_source %>%
    mutate(
      accident_year = lubridate::year(as.Date(`Date of Occurence (D1)`))
    ),
  unlabelled_claims %>% select(heat_prob, pct_bin_fine)
) %>%
  left_join(pi_table %>% select(pct_bin_fine, pi), by = "pct_bin_fine") %>%
  mutate(pi = if_else(is.na(pi), 0, pi))   # bins outside P86-P100 → π = 0

confirmed_heat_meta <- semi_supervised %>%
  filter(confirmed_heat == 1) %>%
  mutate(accident_year = lubridate::year(as.Date(`Date of Occurence (D1)`)))

cat("Total estimated hidden (quantile method):",
    round(sum(raw_unlabelled_with_pi$pi, na.rm = TRUE)), "\n")


###############################################################################
# [VIZ-1] UNDERRECOGNISED HEAT INJURIES BY ACCIDENT YEAR
# FIXED: now uses raw_unlabelled_with_pi (quantile method) not old underreported_year
###############################################################################

# --- Calculate by year ---
potential_by_year_q <- raw_unlabelled_with_pi %>%
  group_by(accident_year) %>%
  summarise(potential_heat = sum(pi, na.rm = TRUE), .groups = "drop")

confirmed_by_year_q <- confirmed_heat_meta %>%
  count(accident_year, name = "confirmed_heat_n")

underreported_year_q <- confirmed_by_year_q %>%
  full_join(potential_by_year_q, by = "accident_year") %>%
  mutate(
    confirmed_heat_n = if_else(is.na(confirmed_heat_n), 0L, confirmed_heat_n),
    potential_heat   = if_else(is.na(potential_heat),   0,  potential_heat),
    total_estimated  = confirmed_heat_n + potential_heat
  ) %>%
  filter(accident_year >= 2008, accident_year <= 2023) %>%
  arrange(accident_year)

writexl::write_xlsx(underreported_year_q, "underreported_by_year_q.xlsx")
cat("[VIZ-1] Saved: underreported_by_year_q.xlsx\n")

# --- Plot ---
plot_year_long <- underreported_year_q %>%
  select(accident_year, confirmed_heat_n, potential_heat) %>%
  pivot_longer(
    cols      = c(confirmed_heat_n, potential_heat),
    names_to  = "type",
    values_to = "n"
  ) %>%
  mutate(type = factor(type,
                       levels = c("potential_heat", "confirmed_heat_n"),
                       labels = c("Estimated hidden heat injuries",
                                  "Confirmed heat injuries")))

p_year <- ggplot(plot_year_long,
                 aes(x = factor(accident_year), y = n, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(
    data = underreported_year_q,
    aes(x = factor(accident_year), y = total_estimated,
        label = comma(round(total_estimated))),
    inherit.aes = FALSE,
    vjust = -0.4, size = 3, colour = "grey30"
  ) +
  scale_fill_manual(
    values = c(
      "Confirmed heat injuries"        = "#2E75B6",
      "Estimated hidden heat injuries" = "#F4B942"
    )
  ) +
  scale_y_continuous(
    labels = comma_format(),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title    = "Estimated heat-related workers' compensation claims by accident year",
    subtitle = paste0("Confirmed cases (TOOCS-coded) and estimated hidden cases ",
                      "(manual review precision x model-scored unlabelled claims)"),
    x        = "Accident year",
    y        = "Number of claims",
    fill     = NULL,
    caption  = paste0(
      "NDS coverage: 2008-09 to 2023-24 financial years. ",
      "Hidden heat estimates derived from manual review of top-15% ",
      "model-scored claims (P86-P100). pi weighted by bin size (~5,129 per bin)."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "top",
    axis.text.x        = element_text(angle = 45, hjust = 1),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(size = 10, colour = "grey40"),
    plot.caption       = element_text(size = 8,  colour = "grey50"),
    panel.grid.major.x = element_blank()
  )

ggsave("plot_underreported_by_year.png", p_year, width = 12, height = 6, dpi = 300)
cat("[VIZ-1] Saved: plot_underreported_by_year.png\n")


###############################################################################
# [VIZ-2] UNDERRECOGNISED HEAT INJURIES BY STATE
# Already correct in original — uses raw_unlabelled_with_pi (quantile method)
###############################################################################

potential_by_state_q <- raw_unlabelled_with_pi %>%
  group_by(`Jurisdiction (SWA2)`) %>%
  summarise(potential_heat = round(sum(pi, na.rm = TRUE)), .groups = "drop")

confirmed_by_state_q <- confirmed_heat_meta %>%
  count(`Jurisdiction (SWA2)`, name = "confirmed_heat_n")

underreported_state_q <- confirmed_by_state_q %>%
  full_join(potential_by_state_q, by = "Jurisdiction (SWA2)") %>%
  mutate(
    confirmed_heat_n     = if_else(is.na(confirmed_heat_n), 0L, confirmed_heat_n),
    potential_heat       = if_else(is.na(potential_heat),   0,  potential_heat),
    total_estimated      = confirmed_heat_n + potential_heat,
    hidden_per_confirmed = round(potential_heat / pmax(confirmed_heat_n, 1), 1)
  ) %>%
  filter(!is.na(`Jurisdiction (SWA2)`)) %>%
  arrange(desc(total_estimated))

cat("\n=== Underrecognised by state (quantile method) ===\n")
print(underreported_state_q)
cat("Sum of all state totals:", sum(underreported_state_q$total_estimated), "\n")

writexl::write_xlsx(underreported_state_q, "underreported_by_state_q.xlsx")
cat("[VIZ-2] Saved: underreported_by_state_q.xlsx\n")

plot_state_long <- underreported_state_q %>%
  select(`Jurisdiction (SWA2)`, confirmed_heat_n, potential_heat) %>%
  pivot_longer(
    cols      = c(confirmed_heat_n, potential_heat),
    names_to  = "type",
    values_to = "n"
  ) %>%
  mutate(
    type  = factor(type,
                   levels = c("potential_heat", "confirmed_heat_n"),
                   labels = c("Estimated hidden heat injuries",
                              "Confirmed heat injuries")),
    state = factor(`Jurisdiction (SWA2)`,
                   levels = underreported_state_q$`Jurisdiction (SWA2)`)
  )

p_state <- ggplot(plot_state_long,
                  aes(x = state, y = n, fill = type)) +
  geom_col(width = 0.65) +
  geom_text(
    data = underreported_state_q,
    aes(x = `Jurisdiction (SWA2)`, y = total_estimated,
        label = comma(total_estimated)),
    inherit.aes = FALSE,
    hjust = -0.15, size = 3.2, colour = "grey30"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Confirmed heat injuries"        = "#2E75B6",
      "Estimated hidden heat injuries" = "#F4B942"
    )
  ) +
  scale_y_continuous(
    labels = comma_format(),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title    = "Estimated heat-related workers' compensation claims by state",
    subtitle = "Confirmed cases and estimated hidden cases, ordered by total estimated burden",
    x        = NULL,
    y        = "Number of claims",
    fill     = NULL,
    caption  = paste0(
      "Hidden heat estimates derived from manual review of top-15% ",
      "model-scored claims (P86-P100)."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "top",
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(size = 10, colour = "grey40"),
    plot.caption       = element_text(size = 8,  colour = "grey50"),
    panel.grid.major.y = element_blank()
  )

ggsave("plot_underreported_by_state.png", p_state, width = 10, height = 6, dpi = 300)
cat("[VIZ-2] Saved: plot_underreported_by_state.png\n")

###############################################################################
# VIZ3 PREP — Run this BEFORE VIZ3_build_html.R
#
# Requires objects already in memory from main script:
#   raw_unlabelled_with_pi   (has postcode_2021, accident_year, pi)
#   confirmed_heat_meta      (has postcode_2021, accident_year)
#   poa_shp_path             (path to POA shapefile)
#   postcode_mappings_path   (optional, for SA2 names)
#
# Produces objects:
#   poa_geojson_str   — GeoJSON string with all year data embedded as properties
#   total_breaks_r    — quantile colour breaks for total mode
#   ratio_breaks_r    — quantile colour breaks for ratio mode
###############################################################################
sa4_corr_path <- "C:/Users/User/OneDrive - UNSW/Honours/3. Project 1/CG_POA_2021_SA4_2021.xlsx"
sa4_shp_path <- "C:/Users/User/OneDrive - UNSW/Honours/3. Project 1/SA4_2021_AUST_GDA2020.shp"
poa_shp_path <- "C:/Users/User/OneDrive - UNSW/Honours/3. Project 1/Postcode Mappings/POA_2021_AUST_GDA2020_SHP/POA_2021_AUST_GDA2020.shp"
library(dplyr); library(tidyr); library(sf)
library(readxl); library(rmapshaper); library(stringr); library(jsonlite)

# ── 1. Load POA → SA4 correspondence ─────────────────────────────────────────
poa_sa4 <- read_excel(sa4_corr_path) %>%
  transmute(
    claim_poa = str_pad(as.character(POA_CODE_2021), 4, "left", "0"),
    sa4_code  = str_pad(as.character(SA4_CODE_2021), 3, "left", "0"),
    sa4_name  = SA4_NAME_2021,
    ratio     = RATIO_FROM_TO
  ) %>%
  filter(!is.na(ratio), ratio > 0)

cat("Correspondence: ", nrow(poa_sa4), "rows |",
    n_distinct(poa_sa4$claim_poa), "POAs →",
    n_distinct(poa_sa4$sa4_code), "SA4s\n")

# ── 2. Aggregate potential heat by POA × year, then disaggregate to SA4 ──────
potential_poa_year <- raw_unlabelled_with_pi %>%
  mutate(claim_poa = str_pad(as.character(postcode_2021), 4, "left", "0")) %>%
  group_by(claim_poa, accident_year) %>%
  summarise(p = round(sum(pi, na.rm = TRUE)), .groups = "drop") %>%
  filter(accident_year >= 2008, accident_year <= 2023)

confirmed_poa_year <- confirmed_heat_meta %>%
  mutate(claim_poa = str_pad(as.character(postcode_2021), 4, "left", "0")) %>%
  count(claim_poa, accident_year, name = "c") %>%
  filter(accident_year >= 2008, accident_year <= 2023)

poa_year_full <- potential_poa_year %>%
  full_join(confirmed_poa_year, by = c("claim_poa", "accident_year")) %>%
  mutate(
    p = if_else(is.na(p), 0, as.numeric(p)),
    c = if_else(is.na(c), 0L, c)
  )

# Disaggregate POA × year → SA4 × year using RATIO_FROM_TO
sa4_year <- poa_year_full %>%
  left_join(poa_sa4 %>% select(claim_poa, sa4_code, sa4_name, ratio),
            by = "claim_poa", relationship = "many-to-many") %>%
  filter(!is.na(sa4_code)) %>%
  group_by(sa4_code, sa4_name, accident_year) %>%
  summarise(
    p = sum(p * ratio, na.rm = TRUE),
    c = sum(c * ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p = round(p), c = round(c))

# ── 3. All-years totals per SA4 ───────────────────────────────────────────────
sa4_all <- sa4_year %>%
  group_by(sa4_code, sa4_name) %>%
  summarise(
    all_p = sum(p), all_c = sum(c), .groups = "drop"
  ) %>%
  mutate(
    all_t = all_p + all_c,
    all_r = round(all_p / pmax(all_c, 1), 1),
    # NEW metric requested by supervisor: per-1000-claims rate
    # Need total NDS claims per SA4 for denominator — use confirmed + unlabelled
    # as proxy (we don't have total claims per SA4 directly)
  )

cat("SA4 totals sum:", sum(sa4_all$all_t), "| SA4 areas:", nrow(sa4_all), "\n")

# ── 4. Pivot year data wide (same structure as POA version) ───────────────────
sa4_wide <- sa4_year %>%
  pivot_wider(
    id_cols    = sa4_code,
    names_from = accident_year,
    values_from = c(p, c),
    names_glue  = "y{accident_year}_{.value}",
    values_fill = 0
  )

# ── 5. Load + simplify SA4 shapefile ─────────────────────────────────────────
sa4_sf <- sf::st_read(sa4_shp_path, quiet = TRUE) %>%
  # SA4 shapefiles use SA4_CODE21 (without underscore before 21)
  select(SA4_CODE21, geometry) %>%
  mutate(SA4_CODE21 = str_pad(as.character(SA4_CODE21), 3, "left", "0")) %>%
  ms_simplify(keep = 0.05, keep_shapes = TRUE, snap = TRUE) %>%
  sf::st_transform(crs = 4326)

cat("SA4 shapefile:", nrow(sa4_sf), "areas\n")

# ── 6. Join all data to SA4 shapefile ────────────────────────────────────────
sa4_map <- sa4_sf %>%
  left_join(sa4_all,  by = c("SA4_CODE21" = "sa4_code")) %>%
  left_join(sa4_wide, by = c("SA4_CODE21" = "sa4_code")) %>%
  mutate(
    sa4_name = if_else(is.na(sa4_name), "", as.character(sa4_name)),
    across(where(is.numeric), ~if_else(is.na(.), 0, .))
  ) %>%
  # Rename to POA_CODE21 so VIZ3_build_html.R works unchanged
  rename(POA_CODE21 = SA4_CODE21)

# ── 7. Export GeoJSON string (feeds into VIZ3_build_html.R unchanged) ─────────
tmp <- tempfile(fileext = ".geojson")
sf::st_write(sa4_map, tmp, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
poa_geojson_str <- paste(readLines(tmp, warn = FALSE), collapse = "")
file.remove(tmp)
cat("SA4 GeoJSON:", round(nchar(poa_geojson_str) / 1e6, 1), "MB\n")

# ── 8. Colour breaks for SA4 data ────────────────────────────────────────────
nz_t <- sa4_all$all_t[sa4_all$all_t > 0]
nz_r <- sa4_all$all_r[sa4_all$all_r > 0]
total_breaks_r <- as.numeric(
  c(0, quantile(nz_t, c(.2,.4,.6,.8,.95), na.rm=TRUE), max(nz_t, na.rm=TRUE)))
ratio_breaks_r <- as.numeric(
  c(0, quantile(nz_r, c(.2,.4,.6,.8,.95), na.rm=TRUE), max(nz_r, na.rm=TRUE)))

cat("total_breaks_r:", round(total_breaks_r), "\n")
cat("ratio_breaks_r:", round(ratio_breaks_r, 1), "\n")
cat("Data prep complete — now run VIZ3_build_html.R\n")

# ── 9. Export SA4 summary table ───────────────────────────────────────────────
writexl::write_xlsx(
  sa4_all %>%
    arrange(desc(all_t)) %>%
    transmute(
      `SA4 Code`                              = sa4_code,
      `SA4 Name`                              = sa4_name,
      `Confirmed heat injuries`               = all_c,
      `Estimated hidden heat (Sum n_i * p_i)` = all_p,
      `Total estimated`                       = all_t,
      `Ratio (hidden / confirmed)`            = all_r
    ),
  "underreported_by_SA4.xlsx"
)
cat("Saved: underreported_by_SA4.xlsx\n")

###############################################################################
# VIZ3 BUILD HTML — Run AFTER VIZ3_prep.R
# Fix: controls rendered as plain HTML divs (not JS-built innerHTML)
#      so option value="..." double quotes never sit inside a JS string
###############################################################################

# Computes total NDS claims per SA4 as denominator for rate per 1,000

total_nds_sa4 <- swa_data_filtered %>%
  mutate(claim_poa = str_pad(as.character(postcode_2021), 4, "left", "0")) %>%
  count(claim_poa, name = "n_total") %>%
  left_join(poa_sa4 %>% select(claim_poa, sa4_code, ratio),
            by = "claim_poa", relationship = "many-to-many") %>%
  filter(!is.na(sa4_code)) %>%
  group_by(sa4_code) %>%
  summarise(total_nds = round(sum(n_total * ratio, na.rm = TRUE)), .groups = "drop")

sa4_map <- sa4_map %>%
  left_join(total_nds_sa4 %>% rename(POA_CODE21 = sa4_code), by = "POA_CODE21") %>%
  mutate(
    total_nds = if_else(is.na(total_nds), 1, as.numeric(total_nds)),
    all_rate  = round(all_t / pmax(total_nds, 1) * 1000, 2)
  )

nz_rate <- sa4_map$all_rate[sa4_map$all_rate > 0]
rate_breaks_r <- as.numeric(
  c(0, quantile(nz_rate, c(.2,.4,.6,.8,.95), na.rm=TRUE), max(nz_rate, na.rm=TRUE)))

tmp <- tempfile(fileext = ".geojson")
sf::st_write(sa4_map, tmp, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
poa_geojson_str <- paste(readLines(tmp, warn = FALSE), collapse = "")
file.remove(tmp)

# ── BUILD SCRIPT ─────────────────────────────────────────────────────────────
library(jsonlite)

total_breaks_js <- jsonlite::toJSON(total_breaks_r, auto_unbox = TRUE)
ratio_breaks_js <- jsonlite::toJSON(ratio_breaks_r, auto_unbox = TRUE)
rate_breaks_js  <- jsonlite::toJSON(rate_breaks_r,  auto_unbox = TRUE)

make_opts <- function(selected) {
  paste0(sapply(2008:2023, function(y)
    sprintf("<option value='%d'%s>%d</option>",
            y, if (y == selected) " selected" else "", y)), collapse = "")
}
opts_from <- make_opts(2008)
opts_to   <- make_opts(2023)

html_template <- r"[<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Heat Injury Map - Australia (SA4)</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<style>
  *,*::before,*::after { box-sizing:border-box; margin:0; padding:0; }
  html,body { width:100%; height:100%; }
  #map { width:100%; height:100vh; }
 
  .map-panel {
    position:absolute; z-index:1000; background:white;
    border-radius:6px; font-family:Arial,sans-serif; font-size:12px;
    line-height:1.6; box-shadow:0 1px 8px rgba(0,0,0,.35);
  }
  #header-panel  { top:10px; left:54px; padding:6px 12px; max-width:420px; }
  #control-panel { top:10px; right:10px; padding:10px 14px 12px; min-width:250px; }
  #legend-panel  { bottom:24px; right:10px; padding:8px 12px; min-width:150px; }
  #tooltip-panel {
    position:absolute; z-index:9999; background:white;
    border-radius:6px; box-shadow:0 2px 10px rgba(0,0,0,.25);
    padding:8px 12px; font-family:Arial,sans-serif; font-size:12px;
    pointer-events:none; max-width:280px; line-height:1.5; display:none;
  }
 
  #control-panel h3 { font-size:13px; color:#222; margin-bottom:8px; }
  .sec-label { font-weight:bold; color:#444; margin-bottom:3px; display:block; }
  #control-panel label { display:block; cursor:pointer; margin-bottom:2px; }
  .year-row { display:flex; align-items:center; gap:6px; margin-bottom:5px; }
  .year-row select {
    flex:1; padding:3px 5px; border:1px solid #ccc;
    border-radius:3px; font-size:12px; cursor:pointer;
  }
  #range-display { color:#1a73e8; font-weight:bold; font-size:12px; }
  #mode-desc {
    margin-top:6px; font-size:10px; color:#888;
    border-top:1px solid #eee; padding-top:5px;
  }
 
  .legend-title { font-weight:bold; margin-bottom:4px; font-size:12px; }
  .legend-item  { display:flex; align-items:center; gap:6px; margin-bottom:1px; }
  .legend-swatch { width:14px; height:14px; border-radius:2px; flex-shrink:0; }
 
  .tt-label { font-size:13px; font-weight:bold; color:#00274D; }
  .tt-sa4   { font-size:11px; color:#555; margin-bottom:3px; }
  .tt-hr    { border:none; border-top:1px solid #e0e0e0; margin:4px 0; }
  .tt-metric{ font-size:12px; margin-bottom:1px; }
  .tt-sub   { font-size:10px; color:#888; margin-top:3px; font-style:italic; }
</style>
</head>
<body>
 
<div id="map"></div>
 
<!-- Header -->
<div class="map-panel" id="header-panel">
  <b>Heat-related workers&#8217; compensation claims by SA4 region</b><br>
  <span style="font-size:10px;color:#666;">
    Australia 2008&#8211;2023 &nbsp;|&nbsp; P86&#8211;P100 quantile method &nbsp;|&nbsp; SA4: Statistical Area Level 4 (ABS 2021)
  </span>
</div>
 
<!-- Controls -->
<div class="map-panel" id="control-panel">
  <h3>Map Controls</h3>
 
  <span class="sec-label">Display Mode</span>
  <label>
    <input type="radio" name="dmode" value="total" checked onchange="setMode(this.value)">
    &nbsp;Total estimated heat injuries
  </label>
  <label>
    <input type="radio" name="dmode" value="ratio" onchange="setMode(this.value)">
    &nbsp;Hidden &divide; confirmed ratio
  </label>
  <label style="margin-bottom:10px;">
    <input type="radio" name="dmode" value="rate" onchange="setMode(this.value)">
    &nbsp;Total estimated per 1,000 NDS claims
  </label>
 
  <span class="sec-label" style="margin-top:6px;">Year Range</span>
  <div class="year-row">
    <select id="year-from" onchange="setRange()"><<OPTS_FROM>></select>
    <span style="color:#666;">to</span>
    <select id="year-to"   onchange="setRange()"><<OPTS_TO>></select>
  </div>
  <div id="range-display">All years (2008&#8211;2023)</div>
  <div id="mode-desc">Showing: total estimated heat injuries (all years)</div>
</div>
 
<!-- Legend -->
<div class="map-panel" id="legend-panel">
  <div id="legend-content"></div>
</div>
 
<!-- Tooltip -->
<div id="tooltip-panel"></div>
 
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script>
 
// ── Embedded data ─────────────────────────────────────────────────────────────
var GEOJSON      = <<GEOJSON>>;
var T_BREAKS     = <<T_BREAKS>>;
var R_BREAKS     = <<R_BREAKS>>;
var RATE_BREAKS  = <<RATE_BREAKS>>;
 
var T_COLORS    = ["#FFFFCC","#FED976","#FEB24C","#FD8D3C","#F03B20","#BD0026"];
var R_COLORS    = ["#EFF3FF","#C6DBEF","#9ECAE1","#6BAED6","#2171B5","#084594"];
var RATE_COLORS = ["#F2F0F7","#DADAEB","#BCBDDC","#9E9AC8","#756BB1","#54278F"];
 
// ── State ─────────────────────────────────────────────────────────────────────
var currentFrom = 2008;
var currentTo   = 2023;
var currentMode = "total";
 
// ── Colour helper ─────────────────────────────────────────────────────────────
function getColor(val, breaks, colors) {
  if (!val || val <= 0) return "#d3d3d3";
  for (var i = breaks.length - 2; i >= 0; i--) {
    if (val > breaks[i]) return colors[Math.min(i + 1, colors.length - 1)];
  }
  return colors[0];
}
 
// ── Aggregate from GeoJSON properties across year range ───────────────────────
function getAgg(props, from, to) {
  var p = 0, c = 0;
  if (from === 2008 && to === 2023) {
    p = props.all_p || 0;
    c = props.all_c || 0;
  } else {
    for (var yr = from; yr <= to; yr++) {
      p += (props["y" + yr + "_p"] || 0);
      c += (props["y" + yr + "_c"] || 0);
    }
  }
  var total_nds = props.total_nds || 1;
  return {
    potential : Math.round(p),
    confirmed : Math.round(c),
    total     : Math.round(p + c),
    ratio     : Math.round(c) > 0 ? Math.round(p / c * 10) / 10 : 0,
    rate      : Math.round((p + c) / total_nds * 1000 * 10) / 10
  };
}
 
// ── Map setup ─────────────────────────────────────────────────────────────────
var map = L.map("map", { zoomControl:true }).setView([-27, 134], 4);
 
L.tileLayer("https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png", {
  attribution: "&copy; OpenStreetMap &copy; CARTO", maxZoom:18
}).addTo(map);
 
// ── Feature style ─────────────────────────────────────────────────────────────
function featureStyle(feature) {
  var agg = getAgg(feature.properties, currentFrom, currentTo);
  var val, breaks, colors;
  if (currentMode === "total") {
    val = agg.total; breaks = T_BREAKS; colors = T_COLORS;
  } else if (currentMode === "ratio") {
    val = agg.ratio; breaks = R_BREAKS; colors = R_COLORS;
  } else {
    val = agg.rate;  breaks = RATE_BREAKS; colors = RATE_COLORS;
  }
  return {
    fillColor   : getColor(val, breaks, colors),
    fillOpacity : 0.75,
    color       : "white",
    weight      : 0.6,
    smoothFactor: 0.5
  };
}
 
// ── Tooltip ───────────────────────────────────────────────────────────────────
var ttDiv = document.getElementById("tooltip-panel");
 
function showTooltip(e, feature) {
  var p   = feature.properties;
  var agg = getAgg(p, currentFrom, currentTo);
  var rng = currentFrom === currentTo
            ? String(currentFrom)
            : currentFrom + "&#8211;" + currentTo;
 
  var metricLine;
  if (currentMode === "total") {
    metricLine = "<span class='tt-metric'><b>" + agg.total.toLocaleString() +
                 "</b> total estimated heat injuries</span>";
  } else if (currentMode === "ratio") {
    metricLine = "<span class='tt-metric'><b>" + agg.ratio +
                 "&times;</b> hidden &divide; confirmed</span>";
  } else {
    metricLine = "<span class='tt-metric'><b>" + agg.rate +
                 "</b> estimated heat injuries per 1,000 NDS claims</span>";
  }
 
  ttDiv.innerHTML =
    "<div class='tt-label'>" + (p.sa4_name || p.POA_CODE21) + "</div>" +
    "<div class='tt-sa4'>SA4 code: " + p.POA_CODE21 + "</div>" +
    "<hr class='tt-hr'>" +
    metricLine +
    "<span class='tt-metric'>Confirmed: <b>" + agg.confirmed.toLocaleString() + "</b>" +
    " &nbsp;|&nbsp; Hidden: <b>" + agg.potential.toLocaleString() + "</b></span>" +
    "<span class='tt-metric'>Ratio: <b>" + agg.ratio + "&times;</b>" +
    " &nbsp;|&nbsp; Rate: <b>" + agg.rate + "</b> per 1,000</span>" +
    "<div class='tt-sub'>Year range: " + rng + "</div>";
 
  ttDiv.style.display = "block";
  moveTooltip(e);
}
 
function moveTooltip(e) {
  ttDiv.style.left = (e.originalEvent.pageX + 14) + "px";
  ttDiv.style.top  = (e.originalEvent.pageY - 28) + "px";
}
function hideTooltip() { ttDiv.style.display = "none"; }
 
// ── GeoJSON layer ─────────────────────────────────────────────────────────────
var geojsonLayer = L.geoJSON(GEOJSON, {
  style: featureStyle,
  onEachFeature: function(feature, layer) {
    layer.on({
      mouseover: function(e) {
        e.target.setStyle({ weight:2, color:"#333", fillOpacity:0.9 });
        showTooltip(e, feature);
      },
      mousemove: function(e) { moveTooltip(e); },
      mouseout : function(e) {
        geojsonLayer.resetStyle(e.target);
        hideTooltip();
      }
    });
  }
}).addTo(map);
 
// ── Update all ────────────────────────────────────────────────────────────────
function updateMap() {
  geojsonLayer.setStyle(featureStyle);
  renderLegend();
  renderDesc();
}
 
// ── Legend ────────────────────────────────────────────────────────────────────
function renderLegend() {
  var breaks, colors, title;
  if (currentMode === "total") {
    breaks = T_BREAKS; colors = T_COLORS;
    title  = "Total estimated<br>heat injuries";
  } else if (currentMode === "ratio") {
    breaks = R_BREAKS; colors = R_COLORS;
    title  = "Hidden &divide; confirmed<br>ratio";
  } else {
    breaks = RATE_BREAKS; colors = RATE_COLORS;
    title  = "Estimated heat injuries<br>per 1,000 NDS claims";
  }
  var h = "<div class='legend-title'>" + title + "</div>";
  h += "<div class='legend-item'><div class='legend-swatch' style='background:#d3d3d3'></div><span>0</span></div>";
  for (var i = 0; i < colors.length; i++) {
    var lo  = i === 0 ? 0 : breaks[i];
    var hi  = breaks[i + 1];
    var lbl = hi !== undefined
              ? (currentMode === "rate" ? lo.toFixed(1) + "&#8211;" + hi.toFixed(1)
                                       : Math.round(lo) + "&#8211;" + Math.round(hi))
              : "&gt;" + (currentMode === "rate" ? lo.toFixed(1) : Math.round(lo));
    h += "<div class='legend-item'><div class='legend-swatch' style='background:" +
         colors[i] + "'></div><span>" + lbl + "</span></div>";
  }
  document.getElementById("legend-content").innerHTML = h;
}
 
// ── Mode description ──────────────────────────────────────────────────────────
function renderDesc() {
  var rng = currentFrom === 2008 && currentTo === 2023
            ? "all years"
            : currentFrom + "&#8211;" + currentTo;
  var modeStr = currentMode === "total" ? "total estimated heat injuries"
              : currentMode === "ratio" ? "hidden divided by confirmed"
              : "estimated heat injuries per 1,000 NDS claims";
  document.getElementById("mode-desc").innerHTML = "Showing: " + modeStr + " (" + rng + ")";
  document.getElementById("range-display").innerHTML =
    currentFrom === 2008 && currentTo === 2023
    ? "All years (2008&#8211;2023)"
    : currentFrom + "&#8211;" + currentTo;
}
 
// ── Control handlers ──────────────────────────────────────────────────────────
function setMode(val) { currentMode = val; updateMap(); }
 
function setRange() {
  var f = parseInt(document.getElementById("year-from").value);
  var t = parseInt(document.getElementById("year-to").value);
  if (f > t) { document.getElementById("year-to").value = f; t = f; }
  currentFrom = f; currentTo = t;
  updateMap();
}
 
// Prevent map drag over panels
["control-panel","legend-panel","header-panel"].forEach(function(id) {
  var el = document.getElementById(id);
  L.DomEvent.disableScrollPropagation(el);
  L.DomEvent.disableClickPropagation(el);
});
 
// Initial render
renderLegend();
renderDesc();
 
</script>
</body>
</html>]"

# ── Substitute placeholders ───────────────────────────────────────────────────
html <- html_template
html <- gsub("<<GEOJSON>>",    poa_geojson_str, html, fixed = TRUE)
html <- gsub("<<T_BREAKS>>",   total_breaks_js, html, fixed = TRUE)
html <- gsub("<<R_BREAKS>>",   ratio_breaks_js, html, fixed = TRUE)
html <- gsub("<<RATE_BREAKS>>",rate_breaks_js,  html, fixed = TRUE)
html <- gsub("<<OPTS_FROM>>",  opts_from,        html, fixed = TRUE)
html <- gsub("<<OPTS_TO>>",    opts_to,           html, fixed = TRUE)

writeLines(html, "index.html", useBytes = TRUE)
cat(sprintf("[VIZ-3 SA4] Saved: index.html (%.1f MB)\n",
            file.info("index.html")$size / 1e6))
cat("Deploy: drag index.html to app.netlify.com/drop\n")

###############################################################################
### SENSITIVITY CHECKS
# Independent of the model: apply simple criteria to ALL claims,
# count how many plausibly could be heat-related, compare to model's 1.12%
#
# Requires:
#   swa_data_filtered  — full claims dataset after basic filtering (pre-modelling)
#                        has: wbgt, wbgt_percentile, wbgt_3day_exposure,
#                             heat_industry, heat_month, division, major,
#                             confirmed_heat, wbgt_excess_tmred
#   TOTAL_NDS_CLAIMS   <- 3326133

TOTAL_NDS_CLAIMS  <- 3326133L
MODEL_HIDDEN      <- 37339
MODEL_PREVALENCE  <- MODEL_HIDDEN / TOTAL_NDS_CLAIMS * 100   # 1.12%
CONFIRMED_N       <- nrow(confirmed_heat_meta)
PI_MEAN           <- mean(pi_table$pi)   # mean manual-review precision ≈ 0.49

sense_row <- function(label, filter_expr, data = swa_data_filtered) {
  # eval_tidy uses data masking AND keeps the calling environment → finds thr
  mask <- rlang::eval_tidy(rlang::enquo(filter_expr), data = data)
  n    <- sum(mask, na.rm = TRUE)
  pa   <- round(n * PI_MEAN)
  tibble(
    scenario           = label,
    n_raw              = n,
    prevalence_raw_pct = round(n  / TOTAL_NDS_CLAIMS * 100, 3),
    n_pi_adjusted      = pa,
    prevalence_adj_pct = round(pa / TOTAL_NDS_CLAIMS * 100, 3),
    vs_model_pct       = round(pa / TOTAL_NDS_CLAIMS * 100 - MODEL_PREVALENCE, 3)
  )
}

# ── Print WBGT distribution to guide threshold choice ─────────────────────────
cat("=== WBGT distribution across all unlabelled claims ===\n")
print(round(quantile(swa_data_filtered$wbgt,
                     c(.50,.75,.85,.90,.95,.99), na.rm = TRUE), 1))

cat("\n=== wbgt_3day_exposure distribution ===\n")
print(round(quantile(swa_data_filtered$wbgt_3day_exposure,
                     c(.50,.75,.85,.90,.95,.99), na.rm = TRUE), 1))

cat("\n=== hours_worked_per_week distribution ===\n")
print(table(cut(swa_data_filtered$"Hours Worked per Week (C8)",
                breaks = c(0, 20, 35, 40, 48, Inf),
                labels = c("<20", "20-35", "35-40", "40-48", "48+"))))

cat("\n=== severity_indicator values ===\n")
print(table(swa_data_filtered$`Serious Claims (SWA6)`, useNA = "ifany"))

cat("\n=== bodily_location_of_injury (top 10) ===\n")
print(sort(table(swa_data_filtered$`Bodily Location of Injury (D5)`),
           decreasing = TRUE)[1:10])

# ── SENSE CHECK SCENARIOS ─────────────────────────────────────────────────────
# PART A — WBGT THRESHOLD RANGE SCAN
# Fixed criteria: heat_industry == 1 AND heat_occupations == 1
# Vary: wbgt threshold from 25°C to 32°C
# Literature anchors:
#   25°C  — ISO 7243 action limit for HEAVY unacclimatised work
#   27.5°C — ACGIH TLV for HEAVY acclimatised work (continuous)
#   28°C  — ISO 7243 / ACGIH action limit for MODERATE outdoor work (most cited)
#   32°C  — "very hot for outdoor workers" (your script note)
wbgt_thresholds <- c(25, 26, 27, 27.5, 28, 29, 30, 32)

results_wbgt <- bind_rows(lapply(wbgt_thresholds, function(thr) {
  n <- with(swa_data_filtered,
            sum(wbgt > thr & heat_industry == 1 & heat_occupation == 1,
                na.rm = TRUE))
  tibble(
    scenario       = sprintf("WBGT > %.1f°C + outdoor industry + occupation", thr),
    n_passing      = n,
    prevalence_pct = round(n / TOTAL_NDS_CLAIMS * 100, 3)
  )
}))

cat(sprintf("Model estimate for comparison: %.3f%%\n\n", MODEL_PREVALENCE))
print(results_wbgt)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(results_wbgt, aes(x = wbgt_thresholds, y = prevalence_pct)) +
  geom_line(colour = "#00274D", linewidth = 1) +
  geom_point(size = 3.5, colour = "#00274D") +
  geom_hline(yintercept = MODEL_PREVALENCE, linetype = "dashed",
             colour = "#C0392B", linewidth = 0.8) +
  annotate("text", x = 25.1, y = MODEL_PREVALENCE + 0.05,
           label = sprintf("Model: %.2f%%", MODEL_PREVALENCE),
           colour = "#C0392B", size = 3.5, hjust = 0) +
  geom_vline(xintercept = 28, linetype = "dotted", colour = "grey50") +
  annotate("text", x = 28, y = max(results_wbgt$prevalence_pct) * 0.97,
           label = "28°C\n(ISO moderate)", size = 2.8, colour = "grey40",
           hjust = -0.1) +
  scale_x_continuous(breaks = wbgt_thresholds) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Sensitivity Check: WBGT Threshold vs. Implied Prevalence",
    subtitle = "Criteria: WBGT > threshold + outdoor industry + outdoor occupation\nRaw count ÷ 3,326,133 total NDS claims — no model scores used",
    x        = "WBGT threshold (°C)",
    y        = "Implied prevalence (% of all NDS claims)",
    caption  = "Red dashed = model estimate (1.12%). Criteria-passing claims treated as potential heat injuries."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 10))

ggsave("sensitivity_A_wbgt_range.png", width = 10, height = 5.5, dpi = 300)

# PART B: 3-day cumulative heatwave scan
# ── Print distribution to report actual percentile values in output ───────────
cat("=== wbgt_3day_exposure distribution across all claims ===\n")
q_vals <- quantile(swa_data_filtered$wbgt_3day_exposure,
                   probs = c(.50, .75, .85, .90, .95, .99), na.rm = TRUE)
print(round(q_vals, 1))

p85_3d <- as.numeric(q_vals["85%"])
p90_3d <- as.numeric(q_vals["90%"])
p95_3d <- as.numeric(q_vals["95%"])

# ── Define thresholds with literature source labels ───────────────────────────
threshold_def <- tribble(
  ~threshold, ~label,
  75,         "75 (= 25°C × 3 days | ISO lower action limit)",
  81,         "81 (= 27°C × 3 days | ACGIH heavy work limit)",
  84,         "84 (= 28°C × 3 days | ISO moderate work standard)",
  90,         "90 (= 30°C × 3 days | sustained hot)",
  96,         "96 (= 32°C × 3 days | sustained very hot)"
)

# Add percentile-based thresholds (Varghese/BoM relative approach)
pct_defs <- tribble(
  ~threshold,        ~label,
  round(p85_3d, 0), sprintf("%.0f (= P85 of 3-day WBGT | BoM relative threshold)", p85_3d),
  round(p90_3d, 0), sprintf("%.0f (= P90 of 3-day WBGT | BoM relative threshold)", p90_3d),
  round(p95_3d, 0), sprintf("%.0f (= P95 of 3-day WBGT | moderate/high heatwave, Varghese 2019)", p95_3d)
)

# Combine, deduplicate by threshold value
all_thresholds <- bind_rows(threshold_def, pct_defs) %>%
  arrange(threshold) %>%
  distinct(threshold, .keep_all = TRUE)   # remove exact duplicates if P85/90/95 land on 75/81/84/90/96

cat("\n=== Thresholds to be tested ===\n")
print(all_thresholds)

results_3day <- bind_rows(lapply(seq_len(nrow(all_thresholds)), function(i) {
  thr <- all_thresholds$threshold[i]
  lbl <- all_thresholds$label[i]
  
  n <- with(swa_data_filtered,
            sum(wbgt_3day_exposure > thr &
                  heat_industry    == 1   &
                  heat_occupation == 1,
                na.rm = TRUE))
  
  tibble(
    threshold      = thr,
    scenario       = lbl,
    n_passing      = n,
    prevalence_pct = round(n / TOTAL_NDS_CLAIMS * 100, 3)
  )
}))

cat(sprintf("\nModel estimate: %.3f%%\n\n", MODEL_PREVALENCE))
print(results_3day %>% select(scenario, n_passing, prevalence_pct), n = Inf)

# ── Save ──────────────────────────────────────────────────────────────────────
writexl::write_xlsx(results_3day, "sensitivity_B_3day_scan.xlsx")
cat("\nSaved: sensitivity_B_3day_scan.xlsx\n")

# ── Plot ──────────────────────────────────────────────────────────────────────
# Distinguish fixed (literature) vs percentile-based thresholds
results_3day <- results_3day %>%
  mutate(
    type = case_when(
      threshold %in% c(round(p85_3d), round(p90_3d), round(p95_3d)) ~
        "Relative (percentile-based, BoM/Varghese approach)",
      TRUE ~ "Fixed (ISO 7243 / ACGIH standards × 3 days)"
    )
  )

p_3day <- ggplot(results_3day,
                 aes(x = threshold, y = prevalence_pct,
                     colour = type, shape = type)) +
  geom_line(aes(group = 1), colour = "grey70", linewidth = 0.7) +
  geom_point(size = 4) +
  geom_hline(yintercept = MODEL_PREVALENCE, linetype = "dashed",
             colour = "#C0392B", linewidth = 0.8) +
  annotate("text", x = min(results_3day$threshold) + 0.5,
           y = MODEL_PREVALENCE + 0.04,
           label = sprintf("Model estimate: %.2f%%", MODEL_PREVALENCE),
           colour = "#C0392B", size = 3.5, hjust = 0) +
  
  # ISO 28°C × 3 reference line
  geom_vline(xintercept = 84, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 84, y = max(results_3day$prevalence_pct) * 0.98,
           label = "84\n(28°C × 3\nISO moderate)",
           size = 2.8, colour = "grey40", hjust = -0.1, vjust = 1) +
  
  scale_colour_manual(values = c(
    "Fixed (ISO 7243 / ACGIH standards × 3 days)"                    = "#00274D",
    "Relative (percentile-based, BoM/Varghese approach)"              = "#E67E22"
  )) +
  scale_shape_manual(values = c(
    "Fixed (ISO 7243 / ACGIH standards × 3 days)"                    = 16,
    "Relative (percentile-based, BoM/Varghese approach)"              = 17
  )) +
  scale_x_continuous(breaks = results_3day$threshold) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  
  labs(
    title    = "Part B: Implied Prevalence vs. 3-Day Cumulative WBGT Threshold",
    subtitle = paste0(
      "Fixed criteria: outdoor industry (heat_industry=1) + outdoor occupation (heat_occupations=1)\n",
      "Blue = fixed ISO/ACGIH standard thresholds × 3 days  |  ",
      "Orange = relative percentile thresholds (BoM/Varghese 2019 approach)\n",
      sprintf("P85=%.0f, P90=%.0f, P95=%.0f",
              p85_3d, p90_3d, p95_3d)
    ),
    x       = "wbgt_3day_exposure threshold (sum of 3 daily WBGT values, °C)",
    y       = "Implied prevalence (% of all NDS claims)",
    colour  = "Threshold type",
    shape   = "Threshold type",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    plot.caption     = element_text(colour = "grey50", size = 8),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

ggsave("sensitivity_B_3day_scan.png", p_3day,
       width = 11, height = 6.5, dpi = 300)
cat("Saved: sensitivity_B_3day_scan.png\n")

###############################################################################
## SEIFA Scores Analysis
library(dplyr); library(tidyr); library(readxl)
library(ggplot2); library(scales); library(writexl); library(stringr)

seifa_path <- "C:/Users/User/OneDrive - UNSW/Honours/3. Project 1/Postal Area, Indexes, SEIFA 2021.xlsx"

all_raw <- read_excel(seifa_path, sheet = "Table 1", skip = 5,
                      col_names = TRUE, col_types = "text")

seifa <- all_raw %>%
  rename(
    claim_poa        = `2021 Postal Area (POA) Code`,
    irsd_score       = `Score...2`,
    irsd_decile      = `Decile...3`,
    ier_score        = `Score...4`,
    ier_decile       = `Decile...5`,
    ieo_score        = `Score...6`,
    ieo_decile       = `Decile...7`,
    irsad_score      = `Score...8`,
    irsad_decile     = `Decile...9`,
    population       = `Usual Resident Population`,
    caution_flag     = `Data should be used with caution - area not well represented by SA1s`,
    cross_state_flag = `POA crosses state or territory boundaries`
  ) %>%
  mutate(
    claim_poa     = str_pad(as.character(claim_poa), 4, "left", "0"),
    across(c(irsd_score, irsd_decile, irsad_score, irsad_decile,
             ier_score,  ier_decile,  ieo_score,   ieo_decile,
             population),
           ~suppressWarnings(as.numeric(.))),
    irsd_decile   = as.integer(irsd_decile),
    irsad_decile  = as.integer(irsad_decile),
    irsd_quintile = ceiling(irsd_decile  / 2),
    irsad_quintile= ceiling(irsad_decile / 2),
    caution       = !is.na(caution_flag)     & caution_flag     == "Y",
    cross_state   = !is.na(cross_state_flag) & cross_state_flag == "Y"
  ) %>%
  filter(!is.na(irsd_decile)) %>%
  select(-caution_flag, -cross_state_flag)

cat("SEIFA loaded:", nrow(seifa), "postcodes\n")
cat("IRSD deciles:", min(seifa$irsd_decile), "-", max(seifa$irsd_decile), "\n")
cat("Caution:", sum(seifa$caution), "| Cross-state:", sum(seifa$cross_state), "\n")
print(head(seifa, 3))

# Attach to claims (2014-2023)
swa_seifa <- swa_data_filtered %>%
  filter(accident_year %in% 2014:2023) %>%
  mutate(claim_poa = str_pad(as.character(postcode_2021), 4, "left", "0")) %>%
  left_join(seifa, by = "claim_poa")

cat(sprintf("\nCoverage in 2014-2023 NDS claims: %.1f%%\n",
            mean(!is.na(swa_seifa$irsd_decile)) * 100))

# ── 4. Confirmed heat injuries by SEIFA decile ────────────────────────────────
confirmed_by_decile <- swa_seifa %>%
  filter(confirmed_heat == 1, !is.na(irsd_decile)) %>%
  group_by(irsd_decile) %>%
  summarise(confirmed = n(), .groups = "drop")

# Total NDS claims by SEIFA decile (denominator for rate)
total_by_decile <- swa_seifa %>%
  filter(!is.na(irsd_decile)) %>%
  group_by(irsd_decile) %>%
  summarise(total_nds = n(), .groups = "drop")

STUDY_YEARS <- 2014:2023
# ── 5. Estimated hidden heat injuries by SEIFA decile ────────────────────────
hidden_by_decile <- raw_unlabelled_with_pi %>%
  filter(accident_year %in% STUDY_YEARS) %>%
  mutate(claim_poa = str_pad(as.character(postcode_2021), 4, "left", "0")) %>%
  left_join(seifa %>% select(claim_poa, irsd_decile), by = "claim_poa") %>%
  filter(!is.na(irsd_decile)) %>%
  group_by(irsd_decile) %>%
  summarise(hidden_est = round(sum(pi, na.rm = TRUE)), .groups = "drop")

# ── 6. Combine into summary table by decile ───────────────────────────────────
seifa_summary <- total_by_decile %>%
  left_join(confirmed_by_decile, by = "irsd_decile") %>%
  left_join(hidden_by_decile,    by = "irsd_decile") %>%
  mutate(
    confirmed  = if_else(is.na(confirmed),  0L, confirmed),
    hidden_est = if_else(is.na(hidden_est), 0,  hidden_est),
    total_est  = confirmed + hidden_est,
    # Rate per 1,000 NDS claims
    rate_per_1000 = round(total_est / total_nds * 1000, 2),
    # Ratio hidden to confirmed
    ratio      = round(hidden_est / pmax(confirmed, 1), 1),
    decile_label = paste0("Decile ", irsd_decile,
                          if_else(irsd_decile == 1, "\n(most disadvantaged)",
                                  if_else(irsd_decile == 10, "\n(most advantaged)", "")))
  ) %>%
  arrange(irsd_decile)

cat("\n=== Summary by IRSD Decile (2014-2023) ===\n")
print(seifa_summary %>%
        select(irsd_decile, total_nds, confirmed, hidden_est, total_est,
               rate_per_1000, ratio))

writexl::write_xlsx(seifa_summary, "seifa_decile_summary_2014_2023.xlsx")
cat("Saved: seifa_decile_summary_2014_2023.xlsx\n")

# ── 7. PLOT A: Total estimated heat injuries by SEIFA decile (stacked) ─────
plot_stack <- seifa_summary %>%
  select(irsd_decile, decile_label, confirmed, hidden_est) %>%
  pivot_longer(c(confirmed, hidden_est),
               names_to = "type", values_to = "n") %>%
  mutate(type = factor(type,
                       levels = c("hidden_est", "confirmed"),
                       labels = c("Estimated hidden", "Confirmed")))

p_stack <- ggplot(plot_stack,
                  aes(x = factor(irsd_decile), y = n, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(data = seifa_summary,
            aes(x = factor(irsd_decile), y = total_est,
                label = comma(total_est)),
            inherit.aes = FALSE, vjust = -0.4, size = 3, colour = "grey30") +
  scale_fill_manual(values = c("Confirmed" = "#2E75B6",
                               "Estimated hidden" = "#F4B942")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  scale_x_discrete(labels = c("1\n(most\ndisadvantaged",
                              "2","3","4","5","6","7","8","9","10\n(most\nadvantaged)")) +
  labs(
    title    = "Estimated heat-related workers' compensation claims by SEIFA decile",
    subtitle = "IRSD Decile 1 = most disadvantaged | 2014–2023 | Confirmed (blue) + estimated hidden (gold)",
    x        = "IRSD Decile (socio-economic disadvantage)",
    y        = "Number of claims",
    fill     = NULL,
    caption  = "SEIFA 2021 (ABS). Hidden estimates: quantile method (P86-P100), pi from manual review."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

ggsave("seifa_plot_A_stacked.png", p_stack, width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_A_stacked.png\n")

# ── 8. PLOT B: Rate per 1,000 NDS claims by SEIFA decile ─────────────────────
p_rate <- ggplot(seifa_summary,
                 aes(x = factor(irsd_decile), y = rate_per_1000)) +
  geom_col(fill = "#8B0000", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rate_per_1000)),
            vjust = -0.4, size = 3.2, colour = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_x_discrete(labels = c("1\n(most\ndisadvantaged",
                              "2","3","4","5","6","7","8","9","10\n(most\nadvantaged)")) +
  labs(
    title    = "Estimated heat injury rate per 1,000 NDS claims by SEIFA decile",
    subtitle = "Total estimated (confirmed + hidden) ÷ total NDS claims in that decile × 1,000 | 2014–2023",
    x        = "IRSD Decile (socio-economic disadvantage)",
    y        = "Estimated heat injuries per 1,000 NDS claims",
    caption  = "SEIFA 2021 (ABS). Rate standardises for different workforce sizes across deciles."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

ggsave("seifa_plot_B_rate.png", p_rate, width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_B_rate.png\n")

# ── 9. PLOT C: Time series by IRSD quintile (2014-2023) ──────────────────────
# Supervisor asked for trend over 2014-2023 like the state/year analysis

# Confirmed by quintile × year
confirmed_qy <- swa_seifa %>%
  filter(confirmed_heat == 1, !is.na(irsd_quintile)) %>%
  count(accident_year, irsd_quintile, name = "confirmed")

# Hidden by quintile × year
hidden_qy <- raw_unlabelled_with_pi %>%
  filter(accident_year %in% STUDY_YEARS) %>%
  mutate(claim_poa = str_pad(as.character(postcode_2021), 4, "left", "0")) %>%
  left_join(seifa %>% select(claim_poa, irsd_quintile), by = "claim_poa") %>%
  filter(!is.na(irsd_quintile)) %>%
  group_by(accident_year, irsd_quintile) %>%
  summarise(hidden_est = round(sum(pi, na.rm = TRUE)), .groups = "drop")

# Total NDS by quintile × year (denominator)
total_qy <- swa_seifa %>%
  filter(!is.na(irsd_quintile)) %>%
  count(accident_year, irsd_quintile, name = "total_nds")

seifa_yr_quintile <- total_qy %>%
  left_join(confirmed_qy, by = c("accident_year","irsd_quintile")) %>%
  left_join(hidden_qy,    by = c("accident_year","irsd_quintile")) %>%
  mutate(
    confirmed  = if_else(is.na(confirmed),  0L, confirmed),
    hidden_est = if_else(is.na(hidden_est), 0,  hidden_est),
    total_est  = confirmed + hidden_est,
    rate_per_1000 = round(total_est / total_nds * 1000, 2),
    quintile_label = factor(irsd_quintile,
                            labels = c("Q1 (most disadvantaged)","Q2","Q3","Q4","Q5 (most advantaged)"))
  )

p_trend <- ggplot(seifa_yr_quintile,
                  aes(x = accident_year, y = rate_per_1000,
                      colour = quintile_label, group = quintile_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(
    "Q1 (most disadvantaged)" = "#B2182B",
    "Q2"                      = "#EF8A62",
    "Q3"                      = "#F7F7F7",
    "Q4"                      = "#67A9CF",
    "Q5 (most advantaged)"    = "#2166AC"
  )) +
  scale_x_continuous(breaks = 2014:2023) +
  labs(
    title    = "Estimated heat injury rate by SEIFA quintile over time",
    subtitle = "Estimated heat injuries per 1,000 NDS claims | IRSD quintiles | 2014–2023",
    x        = "Accident year",
    y        = "Rate per 1,000 NDS claims",
    colour   = "IRSD Quintile",
    caption  = "SEIFA 2021 (ABS). Rate = (confirmed + estimated hidden) ÷ total NDS claims × 1,000."
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        axis.text.x   = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave("seifa_plot_C_trend.png", p_trend, width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_C_trend.png\n")

# ── 10. PLOT D: WBGT distribution by SEIFA decile ────────────────────────────
wbgt_seifa <- swa_seifa %>%
  filter(!is.na(irsd_decile), !is.na(wbgt)) %>%
  mutate(decile_f = factor(irsd_decile))

p_wbgt <- ggplot(wbgt_seifa, aes(x = decile_f, y = wbgt, fill = decile_f)) +
  geom_violin(trim = TRUE, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8,
               show.legend = FALSE) +
  scale_fill_manual(values = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(10),
                    guide = "none") +
  scale_x_discrete(labels = c("1\n(most\ndisadvantaged",
                              "2","3","4","5","6","7","8","9","10\n(most\nadvantaged)")) +
  labs(
    title    = "WBGT distribution by SEIFA decile",
    subtitle = "Do more disadvantaged areas have systematically higher heat exposure? | 2014–2023",
    x        = "IRSD Decile",
    y        = "WBGT on day of claim (°C)",
    caption  = "Violin = full distribution; box = IQR + median. All NDS claims 2014-2023."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

ggsave("seifa_plot_D_wbgt.png", p_wbgt, width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_D_wbgt.png\n")

# ── 11. Summary statistics table ─────────────────────────────────────────────
seifa_stats <- swa_seifa %>%
  filter(!is.na(irsd_decile)) %>%
  group_by(irsd_decile) %>%
  summarise(
    n_claims       = n(),
    mean_wbgt      = round(mean(wbgt, na.rm = TRUE), 1),
    median_wbgt    = round(median(wbgt, na.rm = TRUE), 1),
    pct_outdoor    = round(mean(heat_industry == 1, na.rm = TRUE) * 100, 1),
    pct_confirmed  = round(mean(confirmed_heat == 1) * 100, 3),
    .groups = "drop"
  ) %>%
  left_join(seifa_summary %>% select(irsd_decile, rate_per_1000, ratio,
                                     total_est),
            by = "irsd_decile")

cat("\n=== SEIFA EDA Summary Statistics ===\n")
print(seifa_stats)

writexl::write_xlsx(
  list("By decile (2014-2023)"  = seifa_summary,
       "By quintile x year"     = seifa_yr_quintile,
       "WBGT stats by decile"   = seifa_stats),
  "seifa_eda_all.xlsx"
)
cat("Saved: seifa_eda_all.xlsx\n")

###############################################################################
# SEIFA ANALYSIS — Part 2
# A: Underrecognition rate by SEIFA decile
# B: Model score distribution by SEIFA decile
#
# Requires in memory:
#   swa_seifa            — claims 2014-2023 with irsd_decile attached
#   unlabelled_claims    — has heat_prob, pct_bin_fine, postcode_2021
#   confirmed_heat_meta  — confirmed Y=1 claims with postcode_2021
#   seifa                — seifa lookup (claim_poa, irsd_decile, irsd_quintile)
###############################################################################

library(dplyr); library(ggplot2); library(scales); library(tidyr)

STUDY_YEARS <- 2014:2023

# ── A. UNDERRECOGNITION RATE BY SEIFA DECILE ─────────────────────────────────
# Underrecognition rate = hidden_est / total_est
# i.e. what proportion of total estimated heat injuries went unrecorded

# Already have seifa_summary from previous code — just add the rate
underrecog_by_decile <- seifa_summary %>%
  mutate(
    underrecog_rate    = round(hidden_est / total_est * 100, 1),
    confirmed_rate     = round(confirmed  / total_est * 100, 1)
  ) %>%
  select(irsd_decile, confirmed, hidden_est, total_est,
         underrecog_rate, confirmed_rate, ratio)

cat("=== Underrecognition rate by IRSD decile ===\n")
print(underrecog_by_decile)

# Plot A1: underrecognition rate as a bar chart
p_underrecog <- ggplot(underrecog_by_decile,
                       aes(x = factor(irsd_decile), y = underrecog_rate)) +
  geom_col(fill = "#2E75B6", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = paste0(underrecog_rate, "%")),
            vjust = -0.4, size = 3.2, colour = "grey30") +
  geom_hline(yintercept = mean(underrecog_by_decile$underrecog_rate),
             linetype = "dashed", colour = "#C0392B", linewidth = 0.8) +
  annotate("text", x = 0.7,
           y = mean(underrecog_by_decile$underrecog_rate) + 0.3,
           label = sprintf("Mean: %.1f%%",
                           mean(underrecog_by_decile$underrecog_rate)),
           colour = "#C0392B", size = 3.2, hjust = 0) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  scale_x_discrete(labels = c("1\n(most\ndisadvantaged",
                              "2","3","4","5","6","7","8","9","10\n(most\nadvantaged)")) +
  labs(
    title    = "Underrecognition rate by IRSD decile",
    subtitle = "Underrecognition rate = estimated hidden / total estimated heat injuries × 100\n2014-2023 | P86-P100 quantile method",
    x        = "IRSD Decile",
    y        = "Underrecognition rate (%)",
    caption  = "Red dashed = mean rate across all deciles."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

ggsave("seifa_plot_underrecog_rate.png", p_underrecog,
       width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_underrecog_rate.png\n")

# Plot A2: hidden:confirmed ratio by decile (companion plot)
p_ratio <- ggplot(underrecog_by_decile,
                  aes(x = factor(irsd_decile), y = ratio)) +
  geom_col(fill = "#8B0000", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(ratio, "x")),
            vjust = -0.4, size = 3.2, colour = "grey30") +
  scale_x_discrete(labels = c("1\n(most\ndisadvantaged",
                              "2","3","4","5","6","7","8","9","10\n(most\nadvantaged)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Hidden-to-confirmed ratio by IRSD decile",
    subtitle = "Number of estimated hidden heat injuries per each confirmed case | 2014-2023",
    x        = "IRSD Decile",
    y        = "Hidden : confirmed ratio",
    caption  = "Higher ratio = greater relative underrecognition compared to confirmed coding."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

ggsave("seifa_plot_hidden_ratio.png", p_ratio,
       width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_hidden_ratio.png\n")

# ── B. MODEL SCORE DISTRIBUTION BY SEIFA DECILE ──────────────────────────────
# Join heat_prob from unlabelled_claims to SEIFA via postcode

score_seifa <- raw_unlabelled_source %>%
  mutate(claim_poa = str_pad(as.character(postcode_2021), 4, "left", "0")) %>%
  bind_cols(unlabelled_claims %>% select(heat_prob)) %>%
  left_join(seifa %>% select(claim_poa, irsd_decile, irsd_quintile),
            by = "claim_poa") %>%
  filter(!is.na(irsd_decile), !is.na(heat_prob))

cat("score_seifa rows:", nrow(score_seifa), "\n")
cat("IRSD deciles present:", sort(unique(score_seifa$irsd_decile)), "\n")

score_summary <- score_seifa %>%
  group_by(irsd_decile) %>%
  summarise(
    n_claims     = n(),
    mean_score   = round(mean(heat_prob), 4),
    median_score = round(median(heat_prob), 4),
    pct_p95plus  = round(mean(heat_prob > quantile(unlabelled_claims$heat_prob,
                                                   0.95), na.rm = TRUE) * 100, 2),
    .groups = "drop"
  )
print(score_summary)

# Plot B1: mean model score by SEIFA decile
p_score_mean <- ggplot(score_summary,
                       aes(x = factor(irsd_decile), y = mean_score)) +
  geom_col(fill = "#F4B942", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.4f", mean_score)),
            vjust = -0.4, size = 3, colour = "grey30") +
  scale_x_discrete(labels = c("1\n(most\ndisadvantaged",
                              "2","3","4","5","6","7","8","9","10\n(most\nadvantaged)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Mean model heat probability score by IRSD decile",
    subtitle = "XGBoost heat probability score (unlabelled claims only) | 2014-2023",
    x        = "IRSD Decile",
    y        = "Mean heat probability score",
    caption  = "Higher score = claim more likely to be an unrecognised heat injury per the model."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

ggsave("seifa_plot_score_mean.png", p_score_mean,
       width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_score_mean.png\n")

# Plot B2: score distribution violin by SEIFA decile
p_score_violin <- ggplot(
  score_seifa %>% filter(heat_prob > 0.05),  # trim near-zero for clarity
  aes(x = factor(irsd_decile), y = heat_prob, fill = factor(irsd_decile))
) +
  geom_violin(trim = TRUE, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               show.legend = FALSE) +
  scale_fill_manual(
    values = colorRampPalette(c("#B2182B","#F7F7F7","#2166AC"))(10)
  ) +
  scale_x_discrete(labels = c("1\n(most\ndisadvantaged",
                              "2","3","4","5","6","7","8","9","10\n(most\nadvantaged)")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Distribution of model heat probability scores by IRSD decile",
    subtitle = "Unlabelled claims scoring > 0.05 shown | 2014-2023",
    x        = "IRSD Decile",
    y        = "Heat probability score",
    caption  = "Violin = full score distribution; box = IQR + median."
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        panel.grid.major.x = element_blank())

ggsave("seifa_plot_score_violin.png", p_score_violin,
       width = 12, height = 6, dpi = 300)
cat("Saved: seifa_plot_score_violin.png\n")


# ── C. Export combined summary ────────────────────────────────────────────────
full_seifa_summary <- underrecog_by_decile %>%
  left_join(score_summary, by = "irsd_decile")

writexl::write_xlsx(full_seifa_summary, "seifa_underrecog_score_summary.xlsx")
cat("Saved: seifa_underrecog_score_summary.xlsx\n")

###############################################################################
# INDUSTRY × YEAR SUMMARY
# Total estimated heat claims (confirmed + hidden) by industry and year
# Includes hidden:confirmed ratio per industry × year
# Industry groupings: ANZSIC Division and Subdivision
###############################################################################

library(dplyr); library(tidyr); library(writexl)

make_industry_summary <- function(industry_var) {
  
  label <- industry_var   # used in sheet names
  
  # ── Hidden cases ───────────────────────────────────────────────────────────
  hidden_ind_yr <- raw_unlabelled_with_pi %>%
    filter(accident_year >= 2008, accident_year <= 2023,
           !is.na(.data[[industry_var]])) %>%
    group_by(industry = .data[[industry_var]], accident_year) %>%
    summarise(hidden_est = round(sum(pi, na.rm = TRUE)), .groups = "drop")
  
  # ── Confirmed cases ────────────────────────────────────────────────────────
  confirmed_ind_yr <- confirmed_heat_meta %>%
    filter(accident_year >= 2008, accident_year <= 2023,
           !is.na(.data[[industry_var]])) %>%
    count(industry = .data[[industry_var]], accident_year, name = "confirmed")
  
  # ── Combine ────────────────────────────────────────────────────────────────
  combined <- confirmed_ind_yr %>%
    full_join(hidden_ind_yr, by = c("industry", "accident_year")) %>%
    mutate(
      confirmed  = if_else(is.na(confirmed),  0L, confirmed),
      hidden_est = if_else(is.na(hidden_est), 0,  hidden_est),
      total_est  = confirmed + hidden_est,
      # ratio: hidden per each confirmed case; NA if confirmed = 0
      ratio      = if_else(confirmed > 0,
                           round(hidden_est / confirmed, 1),
                           NA_real_)
    )
  
  # ── Industry grand totals for row ordering ─────────────────────────────────
  grand_totals <- combined %>%
    group_by(industry) %>%
    summarise(grand_total = round(sum(total_est)), .groups = "drop") %>%
    arrange(desc(grand_total))
  
  row_order <- grand_totals$industry
  
  # ── Pivot: total_est ───────────────────────────────────────────────────────
  pivot_total <- combined %>%
    select(industry, accident_year, total_est) %>%
    pivot_wider(names_from = accident_year, values_from = total_est,
                values_fill = 0, names_sort = TRUE) %>%
    mutate(Grand_Total = rowSums(across(where(is.numeric)))) %>%
    mutate(industry = factor(industry, levels = row_order)) %>%
    arrange(industry) %>%
    mutate(industry = as.character(industry)) %>%
    rename(Industry = industry)
  
  # ── Pivot: hidden:confirmed ratio ─────────────────────────────────────────
  pivot_ratio <- combined %>%
    select(industry, accident_year, ratio) %>%
    pivot_wider(names_from = accident_year, values_from = ratio,
                names_sort = TRUE) %>%
    # Grand average ratio across all years
    left_join(
      combined %>%
        group_by(industry) %>%
        summarise(
          Grand_Avg_Ratio = round(
            sum(hidden_est) / pmax(sum(confirmed), 1), 1
          ),
          .groups = "drop"
        ),
      by = "industry"
    ) %>%
    mutate(industry = factor(industry, levels = row_order)) %>%
    arrange(industry) %>%
    mutate(industry = as.character(industry)) %>%
    rename(Industry = industry)
  
  list(total = pivot_total, ratio = pivot_ratio)
}

# ── Run for Division and Subdivision ─────────────────────────────────────────
out_div  <- make_industry_summary("division")
out_sub  <- make_industry_summary("subdivision")

# ── Export: one Excel file, four sheets ──────────────────────────────────────
writexl::write_xlsx(
  list(
    "Division_Total est."          = out_div$total,
    "Division_Hidden confirmed ratio"    = out_div$ratio,
    "Subdivision_Total est."       = out_sub$total,
    "Subdivision_Hidden confirmed ratio" = out_sub$ratio
  ),
  "industry_year_heat_summary.xlsx"
)
cat("Saved: industry_year_heat_summary.xlsx\n")

# Quick print of division totals
cat("\n=== Division totals (2008-2023) ===\n")
print(out_div$total %>% select(Industry, Grand_Total), n = 25)
cat("\n=== Division hidden:confirmed ratio (all-years avg) ===\n")
print(out_div$ratio %>% select(Industry, Grand_Avg_Ratio), n = 25)