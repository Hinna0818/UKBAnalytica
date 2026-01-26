# UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for UK Biobank RAP Data

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

![package-overview](man/figures/package-overview.png)

**UKBAnalytica** is a high-performance R package for processing UK Biobank (UKB) Research Analysis Platform (RAP) data exports. It provides efficient extraction of diagnosis records from multiple sources and generates Cox regression-ready survival datasets.

## Author

Nan He  
Department of Bioinformatics, School of Basic Medical Sciences, Southern Medical University  
Email: nanh302311@gmail.com

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

### 1. Data Download (RAP Platform)

First, download data from UK Biobank RAP using the included Python scripts:

```bash
# Download demographic/phenotype data
python ukb_data_loader.py demographic --id-file field_ids_demographic.txt -o population.csv

# Download metabolites (non-ratio subset recommended)
python ukb_data_loader.py metabolites --non-ratio -o metabolites.csv

# Download proteomics data
python protein_loader.py -o ./protein_data
```

### 2. Analysis in R

Load your downloaded data and run survival analysis:

```r
library(UKBAnalytica)
library(data.table)

# Load downloaded UKB data
ukb_data <- fread("population.csv")

# Define diseases of interest
diseases <- get_predefined_diseases()[c("AA", "Hypertension", "Diabetes")]

# Generate survival dataset
# - Use self-report for prevalent exclusion (baseline disease history)
# - Use only hospital records for outcome ascertainment (more reliable)
analysis_dt <- build_survival_dataset(
  dt = ukb_data,
  disease_definitions = diseases,
  prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),  # Include self-report
  outcome_sources = c("ICD10", "ICD9", "Death"),                   # Exclude self-report
  primary_disease = "AA"
)

# Cox regression (excluding prevalent cases)
library(survival)
cox_model <- coxph(
  Surv(outcome_surv_time, outcome_status) ~ Hypertension_history + Diabetes_history,
  data = analysis_dt[AA_history == 0]  # Exclude baseline prevalent cases
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

### Variable Preprocessing

| Function | Description |
|:---------|:------------|
| `preprocess_baseline()` | Preprocess baseline characteristics (demographics, lifestyle, etc.) |
| `calculate_blood_pressure()` | Calculate BP from multiple readings |
| `extract_medications()` | Extract medication use (cholesterol, BP, insulin) |
| `calculate_air_pollution()` | Calculate averaged air pollution exposures |
| `get_variable_info()` | View available variables for preprocessing |

## Variable Preprocessing

The package provides standardized preprocessing for common UKB baseline variables:

```r
library(UKBAnalytica)
library(data.table)

# Load UKB data
ukb_data <- fread("population.csv")

# View available variables
get_variable_info("demographics")
get_variable_info("lifestyle")

# Preprocess baseline characteristics
processed <- preprocess_baseline(
  ukb_data,
  variables = c("sex", "age", "ethnicity", "bmi", "smoking", "education"),
  missing_action = "keep"  # or "drop" to remove rows with NA
)

# Calculate blood pressure (combines manual/auto readings)
processed <- calculate_blood_pressure(processed, type = "sbp")
processed <- calculate_blood_pressure(processed, type = "dbp")

# Extract medication use
processed <- extract_medications(processed, c("cholesterol", "blood_pressure"))

# Calculate air pollution exposure
processed <- calculate_air_pollution(processed, c("NO2", "PM2.5"))

# Use custom variable mapping for unlisted fields
custom_mapping <- list(
  my_biomarker = list(ukb_col = "p30000_i0", description = "Custom biomarker")
)
processed <- preprocess_baseline(
  ukb_data,
  variables = c("sex", "age", "my_biomarker"),
  custom_mapping = custom_mapping
)
```

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

The package supports separate source definitions for prevalent case exclusion and outcome ascertainment:

```r
# Main analysis: 
# - Self-report for baseline exclusion (captures pre-existing conditions)
# - Hospital records only for outcome (more reliable dates)
main_dt <- build_survival_dataset(
  ukb_data, 
  diseases,
  prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources = c("ICD10", "Death")
)

# Sensitivity 1: Hospital records only (strictest definition)
strict_dt <- build_survival_dataset(
  ukb_data, 
  diseases,
  prevalent_sources = c("ICD10", "ICD9"),
  outcome_sources = c("ICD10", "ICD9")
)

# Sensitivity 2: All sources including self-report for outcome
# (Note: self-report dates are year-only, less precise)
broad_dt <- build_survival_dataset(
  ukb_data, 
  diseases,
  prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources = c("ICD10", "ICD9", "Self-report", "Death")
)

# Compare case counts across source definitions
comparison <- compare_data_sources(ukb_data, diseases)
```

### Why separate sources for prevalent vs outcome?

| Aspect | Prevalent (History) | Outcome (Incident) |
|:-------|:-------------------|:-------------------|
| **Purpose** | Exclude baseline cases | Define endpoint |
| **Self-report** | ✅ Include (captures pre-existing) | ❌ Exclude (imprecise dates) |
| **Date precision** | Less critical | Critical for survival time |
| **Validation** | Lower bar acceptable | Higher bar required |

## UKB Data Fields

| Data Type | Code Field | Date Field |
|:----------|:-----------|:-----------|
| ICD-10 | p41270 | p41280_a* |
| ICD-9 | p41271 | p41281_a* |
| Self-report | p20002_i*_a* | p20008_i*_a* |
| Death | p40001, p40002 | p40000 |

## License

MIT License © 2024 Nan He

---

## Python Data Loaders (RAP Platform)

This package includes Python scripts for downloading data from UK Biobank Research Analysis Platform (RAP).

### Installation

The Python scripts are located in `inst/python/` after installing the R package. On RAP JupyterLab:

```python
import sys
sys.path.append('/path/to/UKBAnalytica/inst/python')
```

## UK Biobank Data Download Usage Summary

### File Structure
```
inst/
├── python/
│   ├── ukb_data_loader.py          # Demographic & Metabolites downloader (Spark)
│   ├── protein_loader.py           # Proteomics downloader (dx commands)
│   └── field_ids_demographic.txt   # Example demographic field IDs
└── extdata/
    └── metabolites_non_ratio.txt   # Non-ratio metabolites reference (170 fields)
```

### 1. Demographic Data

**Purpose**: Download any UKB fields by specifying field IDs  
**Method**: Spark-based (dxdata)  
**Flexibility**: High - can download any combination of fields

```bash
# Method 1: Specify IDs directly
python ukb_data_loader.py demographic --ids 31,53,21022,21001 -o population.csv

# Method 2: Read IDs from file (recommended for many fields)
python ukb_data_loader.py demographic --id-file field_ids_demographic.txt -o population.csv

# Method 3: Custom ID file
python ukb_data_loader.py demographic --id-file my_custom_ids.txt -o my_data.csv
```

**ID file format** (supports multiple formats):
```txt
# Comments start with #
31      # Sex
53      # Date of assessment
21022   # Age at recruitment

# Comma separated
20116, 20117, 1289

# Space separated  
1299 1309 1319
```

### 2. Metabolites Data (NMR)

**Purpose**: Download NMR metabolomics data  
**Method**: Spark-based (dxdata)  
**Options**: All fields or non-ratio subset

```bash
# All metabolites (251 fields: 20280-20281, 23400-23648)
python ukb_data_loader.py metabolites -o metabolites_all.csv

# Non-ratio metabolites only (170 curated fields)
python ukb_data_loader.py metabolites --non-ratio -o metabolites_non_ratio.csv
```

### 3. Proteomics Data (Olink)

**Purpose**: Download Olink protein expression data  
**Method**: dx extract_dataset commands  
**Features**: Automatic batching, merging, progress tracking

```bash
# Download all proteins (default)
python protein_loader.py -o ./protein_data

# Custom batch size and delay
python protein_loader.py -o ./protein_data --batch-size 100 --delay 3

# Skip merging batch files
python protein_loader.py -o ./protein_data --no-merge

# Different output directory
python protein_loader.py -o /path/to/output
```

## Complete Usage Examples

### Example 1: Basic Demographics + Metabolites + Proteins
```bash
# Download key demographics
python ukb_data_loader.py demographic --ids 31,53,21022,21001 -o demographics.csv

# Download non-ratio metabolites
python ukb_data_loader.py metabolites --non-ratio -o metabolites.csv

# Download proteins
python protein_loader.py -o ./proteins
```

### Example 2: Cardiovascular Study
```bash
# Download CVD-related fields
echo "31,21022,21001,93,94,30870,30780,30760,41270,41280" > cvd_fields.txt
python ukb_data_loader.py demographic --ids "$(cat cvd_fields.txt)" -o cvd_data.csv

# Download metabolites
python ukb_data_loader.py metabolites --non-ratio -o cvd_metabolites.csv
```

## Output Files Summary

| Data Type | Typical Output | Size | Rows | Columns |
|:----------|:---------------|:-----|:-----|:---------|
| Demographics | population.csv | 10-500 MB | ~500K | Variable |
| Metabolites (all) | metabolites.csv | ~300 MB | ~120K | 251 |
| Metabolites (non-ratio) | metabolites_non_ratio.csv | ~200 MB | ~120K | 170 |
| Proteomics | protein_all_merged.csv | ~800 MB | ~50K | ~3000 |

### Python Scripts Summary

| Script | Data Source | Download Method | Best For |
|:-------|:------------|:----------------|:---------|
| `ukb_data_loader.py demographic` | Any UKB field | Spark (fast) | Phenotype data, covariates |
| `ukb_data_loader.py metabolites` | NMR metabolomics | Spark (fast) | Metabolomics studies |
| `protein_loader.py` | Olink proteomics | dx commands | Proteomics studies |

### Common UKB Field IDs Reference

| Category | Field IDs | Description |
|:---------|:----------|:------------|
| **Basic Demographics** | 31, 53, 21022, 21001 | Sex, assessment date, age, BMI |
| **Lifestyle** | 20116, 20117, 1160 | Smoking, alcohol, sleep |
| **Blood Pressure** | 93, 94, 4079, 4080 | Systolic/diastolic BP |
| **Biomarkers** | 30870, 30780, 30760, 30750 | Triglycerides, LDL, HDL, glucose |
| **Hospital Records** | 41270, 41280, 41271, 41281 | ICD-10/9 diagnoses and dates |
| **Death Registry** | 40000, 40001, 40002 | Death date, primary/secondary causes |
| **Brain Imaging** | 24016-24019, 24003-24008 | Brain volumes |

## Troubleshooting

### Common Issues
1. **Spark initialization fails**: Ensure running in RAP JupyterLab
2. **Field not found**: Check field ID format (some need `_i0`, `_a0` suffixes)  
3. **Memory issues**: Use smaller batches for protein download
4. **Timeout**: Increase `--delay` parameter

### Performance Tips
- Use `--non-ratio` for metabolites if you don't need all fields
- Reduce `--batch-size` if protein download fails
- Run downloads during off-peak hours for better performance

