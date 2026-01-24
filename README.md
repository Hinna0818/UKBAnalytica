# UKBAnalytica

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**UKBAnalytica** is a high-performance R package for processing UK Biobank (UKB) Research Analysis Platform (RAP) data exports. It provides efficient extraction of diagnosis records from multiple sources and generates Cox regression-ready survival datasets.

## Features

- 🚀 **High Performance**: Built on `data.table` for efficient processing of large-scale biobank data
- 📊 **Multiple Data Sources**: Supports ICD-10, ICD-9, self-reported illness, and death registry data
- 🔄 **Flexible Case Definitions**: Extract cases from specific sources for main and sensitivity analyses
- 📈 **Survival Analysis Ready**: Generates datasets with proper prevalent/incident case classification
- 📦 **Predefined Diseases**: Includes validated definitions for common cardiovascular and metabolic diseases

## Installation

```r
# Install from GitHub
devtools::install_github("Hinna0818/UKBAnalytica")
```

## Quick Start

```r
library(UKBAnalytica)
library(data.table)

# Load UKB data
ukb_data <- fread("ukb_data.csv")

# Define diseases of interest
diseases <- get_predefined_diseases()[c("AA", "Hypertension", "Diabetes")]

# Generate analysis-ready dataset (wide output by default)
# Hypertension as primary outcome, adjusting for diabetes history
analysis_dt <- build_survival_dataset(
  dt = ukb_data,
  disease_definitions = diseases,
  sources = "ICD10",          # Main analysis with hospital diagnoses only
  primary_disease = "Hypertension"
)

# Run Cox regression
library(survival)
cox_model <- coxph(
  Surv(outcome_surv_time, outcome_status) ~ 
    Diabetes_history,
  data = analysis_dt
)
summary(cox_model)
```

## Core Functions

### Data Extraction

| Function | Description |
|:---------|:------------|
| `parse_icd10_diagnoses()` | Extract ICD-10 hospital diagnosis records |
| `parse_icd9_diagnoses()` | Extract ICD-9 hospital diagnosis records |
| `parse_self_reported_illnesses()` | Extract self-reported illness records |
| `parse_death_records()` | Extract death registry records |

### Survival Analysis

| Function | Description |
|:---------|:------------|
| `build_survival_dataset()` | Default wide table with full cohort and primary outcome follow-up |
| `extract_cases_by_source()` | Extract cases from specified sources |
| `generate_wide_format()` | Create wide-format disease status table |
| `prepare_analysis_dataset()` | Generate analysis-ready dataset with outcome |

### Disease Definitions

| Function | Description |
|:---------|:------------|
| `create_disease_definition()` | Create custom disease definition |
| `get_predefined_diseases()` | Get predefined disease library |
| `combine_disease_definitions()` | Create composite endpoints (e.g., MACE) |

## Output Format

### Wide Format Table

Each disease generates two columns for flexible use, plus outcome columns:

| Column | Description | Use Case |
|:-------|:------------|:---------|
| `{Disease}_history` | 1 if prevalent case (diagnosed before baseline) | **Covariate adjustment** |
| `{Disease}_incident` | 1 if incident case (diagnosed after baseline) | **Outcome variable** |
| `outcome_status` | 1 if primary disease incident case | **Primary outcome** |
| `outcome_surv_time` | Follow-up time (years) for primary disease | **Primary outcome** |

### Analysis-Ready Dataset

```r
# Example output structure
#    eid  Hypertension_history  Diabetes_history  outcome_status  outcome_surv_time
# 1: 1001                    1                 0               0              14.46
# 2: 1002                    0                 1               1               8.35
```

## Sensitivity Analysis

```r
# Main analysis: ICD-10 only
main_result <- extract_cases_by_source(ukb_data, diseases, sources = "ICD10")

# Sensitivity 1: Hospital records (ICD-10 + ICD-9)
sens1_result <- extract_cases_by_source(ukb_data, diseases, sources = c("ICD10", "ICD9"))

# Sensitivity 2: All sources
sens2_result <- extract_cases_by_source(
  ukb_data, diseases, 
  sources = c("ICD10", "ICD9", "Self-report")
)

# Compare case counts across sources
comparison <- compare_data_sources(ukb_data, diseases)
print(comparison)
```

## UKB Data Fields

| Data Type | Code Field | Date Field |
|:----------|:-----------|:-----------|
| ICD-10 | p41270 | p41280_a* |
| ICD-9 | p41271 | p41281_a* |
| Self-report | p20002_i*_a* | p20008_i*_a* |
| Death | p40001, p40002 | p40000 |

## License

MIT License © 2024 Nan He
