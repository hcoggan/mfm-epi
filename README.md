# Penn Maternal-Fetal Monitoring: Epidemiological Analysis

## Overview

This R script performs an epidemiological analysis of racial disparities in delivery methods and neonatal outcomes (acidemia, Apgar scores) using birth records from the Hospital of the University of Pennsylvania (HUP). The primary focus is on comparing outcomes (delivery method and fetal acidemia) between Black and white patients, controlling for a range of clinical and demographic confounders.

## Research Questions

1. Do delivery method rates differ by race, after adjusting for clinical confounders?
2. Do neonatal acidemia rates and depressed Apgar scores differ by race?
3. Do these racial disparities persist within subgroups (e.g., by delivery method, gestational age, parity)?
4. Are there racial differences in continuous cord pH values or rupture-to-delivery time?

## Dependencies

The following R packages are required:

```r
# Core data manipulation
dplyr, data.table, tibble, tidyr, janitor, readr, stringr, forcats, lubridate

# Statistical modelling
survival, glmnet, xgboost, car, survey, weights, pROC

# Matching
MatchIt, cobalt

# Epidemiology / meta-analysis
metafor, tidycmprsk, comorbidity, exact2x2

# Visualization
ggplot2, ggsurvfit, ggridges, ggraph, ggdendro, ggsignif, cowplot, patchwork

# Reporting
gtsummary, knitr, yardstick, broom

# Utilities
fuzzyjoin, zipcodeR, igraph, tidygraph, httpgd, testit
```

## Configuration

At the top of the script, two key thresholds define the **primary analysis** vs. sensitivity analyses:

| Parameter | Primary Analysis | Sensitivity Analysis |
|---|---|---|
| `pH_threshold` | `7.10` | varies |
| `base_excess_threshold` | `Inf` (pH only) | finite cutoff |

The working directory points to the project root: `/Volumes/chip-lacava/Groups/mfm/penn-maternal-fetal-monitoring/epi-work/mar-res-main/`.

## Data Sources

| File | Description |
|---|---|
| `combined_covariates-2025-11-03.csv` | Raw birth records with demographics, delivery details, and clinical variables |
| `CTGData_labs.csv` | Cord blood gas lab results (pH and base excess) |
| `id_map/...` | ID mapping linking CTG trace filenames to maternal MRNs and timestamps |
| `id_map/CTGData_offsets/0[1-9]_offsetmap.csv` | Timestamp offsets for CTG trace files |

## Pipeline

### Stage 1: Covariate Preprocessing
**Output:** `preprocessed-births-without-outcomes.csv`

- Deduplicates records by baby MRN
- Derives: gestational age (weeks), pregnancy term category, race, delivery method, maternal age, parity, season/time of delivery, stage 1/2 labour durations
- Classifies C-sections as: *Nonreassuring EFM*, *Labor Arrest*, or *Planned*
- Filters to HUP deliveries from April 2017 onward
- Excludes multiple births (same mother, different baby, within 30 days)
  

### Stage 2: Linking Acidemia Outcomes
**Output:** `preprocessed-births-with-outcomes.csv`

- Joins lab results (cord pH, cord base excess) to birth records via maternal MRN and timestamp proximity (within 30 days)
- Identifies prior C-sections both directly (chart flags) and indirectly (delivery history)
- Defines **acidemia** as cord pH ≤ 7.10 (primary) or pH ≤ threshold AND base excess < threshold (sensitivity)
- Removes births with missing gestational age, pH, or base excess
- - Retains only the first birth per mother in the dataset

### Stage 3: Cohort Description
**Output:** `concise_summary_table.csv`, `all_complications.csv`, `unadjusted_apgar.pdf`

- Summarises cohort by delivery method, maternal age, gestational age, race, parity, prior C-section, and outcomes
- Reports median (IQR) for continuous variables
- Plots 5-minute Apgar score distribution

### Stage 4: Unadjusted Comparisons
**Outputs:** `unadjusted_outcome_rates_by_method.pdf`, `unadjusted_method_by_race.pdf`, `unadjusted_outcome_by_race.pdf`, `unadjusted_outcome_by_race_and_delivery_method.pdf`

- **Part A:** Unadjusted outcome rates (acidemia, depressed Apgar) by delivery method vs. vaginal reference
- **Part B:** Unadjusted delivery method rates by race vs. white reference
- **Part C:** Unadjusted acidemia rates by race vs. white reference
- **Part D:** Acidemia rates by race, stratified by delivery method (Black/white only)

Statistical tests: Fisher's exact (small cells) or chi-squared.

### Stage 5: Logistic Regression
**Outputs:** `delivery-method-lr.csv`, `acidemia-outcome-lr.csv`, `lr-delivery-method.pdf`, `lr-acidemia-outcome.pdf`

- Binarises categorical predictors against reference categories (white race, 2017, spring, afternoon, age 20-29, full term)
- **Part A:** Multivariable logistic regression for each delivery method (vs. vaginal)
- **Part B:** Multivariable logistic regression for acidemia and depressed Apgar, adding delivery method as a predictor
- VIF check (threshold < 5) before model fitting; FDR correction applied across predictors

Reference categories: white race, year 2017, spring, afternoon, age 20–29, full term gestation, vaginal delivery.

### Stage 6: Exact Matching Analysis
**Outputs:** `exact-matching-*.csv`, `forest-plots/*.pdf`

Matches Black patients to white patients (nearest-neighbour, with replacement) on: maternal age, pregnancy term, prior C-section, prior births. Additional matching on delivery method where outcome is the endpoint.

- **Part A:** Unstratified racial differences in delivery method after matching
- **Part B:** Stratified racial differences in delivery method (by age, term, parity, prior C-section)
- **Part C:** Racial differences in acidemia outcomes after matching (including on delivery method)
- **Part D:** Acidemia racial differences stratified by clinical subgroup
- **Part E:** Delivery method racial differences stratified by acidemia status

Results visualised as forest plots using `metafor::forest()`.

### Stage 7: Continuous pH Analysis
**Outputs:** `pH_diffs_overall.pdf`, `pH_diffs_overall_by_dm.pdf`, `dm_diffs_overall_by_pH.pdf`, `time_diffs_overall_by_dm.pdf`

- Matches Black to white patients on clinical confounders + delivery method
- Tests for racial differences in mean cord pH (weighted linear regression)
- Plots acidemia category rates by race, overall and stratified by delivery method
- Plots delivery method rates by race, stratified by acidemia category
- Tests for racial differences in rupture-to-delivery time (log-linear model) by delivery method

## Key Variables

| Variable | Description |
|---|---|
| `acidemia` | Cord pH ≤ 7.10 (primary outcome) |
| `depressed_apgar` | 5-minute Apgar score < 7 |
| `delivery_method` | Vaginal / C-Section (Nonreassuring EFM) / C-Section (Labor Arrest) / C-Section (Planned) |
| `race_black` | Binary indicator: Black (1) vs. white (0) |
| `previous_c_section_indicated` | Prior C-section (direct chart flag or inferred from delivery history) |
| `pregnancy_term` | Very preterm (<32w) / Preterm / Early term / Full term / Late term |

## Reproducibility

Set the random seed at the top of the script:
```r
set.seed(210226)
```

All intermediate outputs are written to CSV so individual stages can be re-run independently.
