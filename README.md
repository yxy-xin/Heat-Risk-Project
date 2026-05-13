# Multi-day Heat Patterns and Injury Mechanisms in Australian Workers' Compensation Claims

**Author:** Xinyi Yan  
**Supervisors:** A/Prof Fei Huang, Xi Lin (Taylor Fry)  
**Institution:** UNSW Sydney  

---

## Overview

This project develops a Gradient Boosted Machine (GBM) classifier to identify probable heat-related workers' compensation claims that are systematically missed by existing administrative coding conventions. Using a Positive-Unlabelled (PU) learning framework, the model recovers a substantial population of hidden heat claims from the Safe Work Australia (SWA) National Dataset for Compensation-based Statistics (NDS).

Key outputs include:
- A trained XGBoost classifier with SHAP-based feature importance
- Cluster analysis (PCA + K-means) revealing four distinct heat injury profiles
- Estimated underreporting burden of heat-related claims by year and state
- Manual review workbooks for model validation

---

## Data Sources

> ⚠️ **Data is confidential and not included in this repository.** The NDS data is provided by Safe Work Australia under a Data Sharing Agreement (DSA) with UNSW. Do not commit any raw or derived data files.

| Source | Description |
|--------|-------------|
| Safe Work Australia (SWA) NDS | ~3.3M accepted workers' compensation claims, 2008–09 to 2023–24 |
| Bureau of Meteorology (BoM) | Daily gridded climate data (Tmax, Tmin, vapour pressure), 1984–2025 |

Claims are linked to climate data by 2021 ASGS postcode boundaries and date of occurrence. The key exposure metric is **Wet Bulb Globe Temperature (WBGT)**, estimated as:

```
WBGT = 0.567 × Temperature + 0.393 × Vapour Pressure + 3.94
```

---

## Repository Structure

```
├── scripts/
│   ├── 06_2_FINAL_WBGT_GBM.R          # GBM model training (base)
│   ├── 07_WBGT_TMRED_LAG_GBM.R        # GBM with TMRED lag features
│   └── 10_FINAL_GBM_ALL_RESULTS.R     # Final results, scoring, outputs
├── README.md
└── .gitignore
```

---

## Script Overview

### `10_FINAL_GBM_ALL_RESULTS.R`
Consolidates all downstream results:
- Scores all ~514k unlabelled claims
- Generates SHAP feature importance and dependence plots
- Builds percentile-binned distribution tables (DOC-10)
- Produces fine-grained manual review workbooks (P85–P100)
- Estimates underreported heat injuries by year (DOC-11) and state (DOC-12)

---

## Model Specification

| Parameter | Value |
|-----------|-------|
| Algorithm | XGBoost (`binary:logistic`) |
| Evaluation metric | AUC-PR (precision-recall) |
| Learning rate (eta) | 0.03 |
| Max tree depth | 1 |
| nrounds | 3000 (with 5-fold CV, early stopping at 50 rounds) |
| Train/test split | 80/20 stratified |
| Downsampling | All confirmed positives retained; negatives at 40% |

---

## Key Predictors

- `wbgt` — raw daily WBGT at workplace postcode
- `wbgt_percentile` — rolling 30-day WBGT percentile rank
- `wbgt_excess_tmred` — WBGT above postcode-level TMRED baseline (floored at 0)
- `wbgt_lag1`, `wbgt_lag2` — prior-day and two-day-prior WBGT
- `heat_month` — 1 if month is October–March
- `agency_sub_major`, `mechanism_sub_major`, `nature_minor` — TOOCS injury coding variables
- Occupation and industry indicators (ANZSCO/ANZSIC)

---

## Requirements

R version 4.x or later. Install required packages with:

```r
install.packages(c(
  "tidyverse", "xgboost", "Matrix", "writexl",
  "shapviz", "ggplot2", "irlba", "caret"
))
```

> It is recommended to use `renv` to capture exact package versions:
> ```r
> renv::restore()
> ```

---

## Data Acknowledgement

This research uses data from the National Dataset for Compensation-based Statistics (NDS), 2008–09 to 2023–24 financial years, provided by Safe Work Australia under a Data Sharing Agreement. All outputs are aggregated; no individual-level information is disclosed.
