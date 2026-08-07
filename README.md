# UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for UK Biobank RAP Data Analysis

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/Hinna0818/UKBAnalytica?style=flat)](https://github.com/Hinna0818/UKBAnalytica/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/Hinna0818/UKBAnalytica)](https://github.com/Hinna0818/UKBAnalytica/commits/main)
[![Visits](https://hits.sh/github.com/Hinna0818/UKBAnalytica.svg)](https://hits.sh/github.com/Hinna0818/UKBAnalytica/)
<!-- badges: end -->

<img src="man/figures/logo.png" align="right" height="139" alt="UKBAnalytica logo" />

**UKBAnalytica** is a high-performance R package for working with UK Biobank
Research Analysis Platform (RAP) data inside approved RAP projects. It focuses on standardized
phenotyping, survival-ready datasets, scalable preprocessing, and downstream analysis.

**Attention**
The package **does not ship UK Biobank participant-level source records**; examples
use field IDs, runtime-generated fully synthetic toy data, or user-provided
tables that remain within RAP-controlled storage.

**For details, please visit**: [Full documentation for UKBAnalytica](https://hinna0818.github.io/UKBAnalytica/)

![](docs/assets/workflow.png)


## Installation
`UKBAnalytica` is now available on CRAN, which is stable:
```r
install.packages("UKBAnalytica")
```

You can also install the development version of `UKBAnalytica` from GitHub, which is **more recommended**:

```r
# install.packages("devtools")
devtools::install_github("Hinna0818/UKBAnalytica")
```

## Citation

If you use UKBAnalytica in your study, please cite:

> N He, K Mo, G Yu, F He. UKBAnalytica: an integrated R package for scalable phenotyping and reproducible epidemiological analysis within the UK Biobank Research Analysis Platform. medRxiv. 2026.06.19.26356057. doi: https://doi.org/10.64898/2026.06.19.26356057


## Quick start

```r
library(UKBAnalytica)

## Runtime-generated fully synthetic toy data for examples
## In RAP, replace this with your approved project table
ukb_raw <- ukb_demo(n = 200, seed = 42)

## Retrieve a predefined disease definition.
copd_def <- get_predefined_diseases()["COPD"]

## Build a prospective, survival-ready cohort.
surv_dt <- build_survival_dataset(
  dt = ukb_raw,
  disease_definitions = copd_def,
  prevalent_sources = "ICD10",
  outcome_sources = "ICD10",
  primary_disease = "COPD",
  baseline_col = "p53_i0",
  censor_date = as.Date("2023-10-31"), # replace your censor date
  show_flow = TRUE
)

table(surv_dt$outcome_status, useNA = "ifany")
attr(surv_dt, "participant_flow")

## Build several independent survival endpoints in one pass.
multi_defs <- get_predefined_diseases()[
  c("Hypertension", "T2DM", "Dementia")
]
multi_surv <- build_survival_dataset(
  dt = ukb_raw,
  disease_definitions = multi_defs,
  primary_disease = names(multi_defs),
  show_flow = FALSE
)

## Multi-endpoint output uses disease-specific columns such as
## Hypertension_status / Hypertension_surv_time.
```

## Disease Definition Sources

UKBAnalytica builds disease phenotypes by taking the earliest valid evidence
from multiple UK Biobank sources. Predefined disease definitions are available
through `get_predefined_diseases()`, and custom definitions can be created with
`create_disease_definition()`.

Supported sources include:

- `ICD10`: hospital inpatient diagnosis codes.
- `ICD9`: historical hospital inpatient diagnosis codes.
- `Self-report`: touchscreen/verbal interview illness codes.
- `Death`: primary or contributory death-cause ICD-10 codes.
- `OPCS4`: hospital operative procedure codes.
- `CancerRegistry`: cancer registry diagnosis dates and ICD-10 morphology information.
- `FirstOccurrence`: UKB first-occurrence date fields derived from linked health records.
- `Algorithm`: UKB algorithmically-defined outcomes, such as myocardial infarction, stroke, dementia, asthma, COPD, and Parkinson's disease.

If a selected source is not defined for a disease, it is ignored automatically.
For example, cancer registry evidence is used only when the disease definition
has a `cancer_icd10_pattern`; procedure evidence is used only when
`opcs4_pattern` is defined.

## Minimal RAP Extraction Workflow

Run the following inside a UK Biobank RAP R session. Participant-level data
should remain inside approved RAP projects and RAP-controlled storage.

```r
library(UKBAnalytica)

dataset <- rap_find_dataset()
fields <- rap_list_fields()

meta <- ukb_metadata_setup(fields_df = fields)

ids <- get_variable_set("clinical_core", output = "field_id")

dt <- ukb_extract_fields(
  field_id = ids,
  metadata = meta,
  mode = "sync",
  strip_entity_prefix = FALSE
)

dt <- ukb_decode(dt, metadata = meta)
dt <- ukb_clean_missing(dt, action = "na")

ukb_snapshot(dt, "Clinical core extracted and cleaned")
```

## Primary-care GP query workflow

GP linked records are queried separately from `build_survival_dataset()`.
The workflow filters `gp_clinical` in RAP Spark SQL, reads the much smaller
`gp_registrations` table, standardizes both tables, and returns one row per
participant and disease with a coverage-aware `gp_case` value.

The example below reproduces the migraine code used in the official UK
Biobank Primary Care documentation (`F26..` in Read v2 or CTV3). UKB states
that this code is an illustration of the table structure, not a complete
validated migraine phenotype.

```r
library(UKBAnalytica)
library(sparklyr)

master <- paste0("spark://master:", Sys.getenv("SPARK_MASTER_PORT"))
sc <- spark_connect(master)

# Replace these with the current database and approved cohort in your project.
database_name <- "app12345_20260101000000"
cohort_table <- "analysis_cohort" # Spark table/view with an eid column
cohort <- data.frame(
  eid = c(1000001, 1000002, 1000003),
  study_start = as.Date("2010-01-01"),
  baseline_date = as.Date("2011-01-01"),
  study_end = as.Date("2020-12-31")
)
gp_observation_end <- as.Date("2020-12-31") # example; use your release cut-off

migraine_query <- data.frame(
  disease = c("Migraine", "Migraine"),
  code_system = c("Read2", "CTV3"),
  code = c("F26..", "F26..")
) |>
  parse_gp_query()

# Inspect both the gp_clinical and gp_registrations SQL before execution.
migraine_query |>
  run_gp_workflow(
    database = database_name,
    cohort_table = cohort_table,
    collect = "summary",
    dry_run = TRUE
  )

gp_migraine <- migraine_query |>
  run_gp_workflow(
    connection = sc,
    database = database_name,
    cohort_table = cohort_table,
    collect = "summary",
    participants = cohort,
    observation_end = gp_observation_end,
    window_start = "study_start",
    window_end = "study_end",
    index_date = "baseline_date",
    min_lookback_days = 365,
    min_followup_days = 365,
    min_coverage_fraction = 0.8,
    clinical_output = "gp_migraine_summary.csv",
    registration_output = "gp_registration_records.csv"
  )

gp_migraine[, .(
  eid,
  disease,
  gp_case,
  first_gp_date,
  coverage_status,
  coverage_fraction,
  control_eligible,
  gp_case_reason
)]
```

`cohort_table` causes both GP SQL queries to use a Spark left-semi join before
collection. With `collect = "summary"`, matching and earliest-date/count
aggregation are also performed in Spark; this is the preferred mode for an
end-to-end phenotype workflow. Use `collect = "records"` only when individual
GP rows are required for record-level review.

`gp_case` is intentionally three-state:

- `1`: at least one exact Read v2 / CTV3 code matched;
- `0`: no code matched and `control_eligible = TRUE` for the requested window;
- `NA`: observation-window coverage was insufficient or unknown.

The strict control check merges overlapping registration periods before
calculating observed days. `min_lookback_days` and `min_followup_days` require
continuous coverage immediately around `index_date`; `min_coverage_fraction`
checks the complete requested window. If no window is supplied, the workflow
retains the earlier, less strict registration-based behavior for compatibility.

The lower-level registration workflow is also pipe-friendly:

```r
gp_coverage <- gp_registrations_raw |>
  parse_gp_registrations() |>
  summarise_gp_coverage(
    clinical_records = gp_filtered_records,
    participants = cohort_eids,
    observation_end = gp_observation_end
  )

gp_observability <- gp_registrations_raw |>
  parse_gp_registrations() |>
  assess_gp_observability(
    participants = cohort,
    window_start = "study_start",
    window_end = "study_end",
    index_date = "baseline_date",
    observation_end = gp_observation_end,
    min_lookback_days = 365,
    min_followup_days = 365,
    min_coverage_fraction = 0.8
  )
```

For pre-extracted or synthetic data, the same end-to-end function can be run
without a Spark connection:

```r
gp_result <- migraine_query |>
  run_gp_workflow(
    clinical_data = gp_clinical_raw,
    registration_data = gp_registrations_raw,
    participants = cohort,
    observation_end = gp_observation_end,
    window_start = "study_start",
    window_end = "study_end",
    index_date = "baseline_date"
  )
```

The two RAP outputs are written to the active worker filesystem. Upload them
with `dx upload` if they need to persist as project files. The observation end
must be chosen from the data release/provider information for the active RAP
project; a missing GP code alone is never sufficient evidence that a
participant is disease-free.

Official references: [UK Biobank Primary Care
Data](https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/primary_care_data.pdf),
[`gp_registrations` record
table](https://biobank.ndph.ox.ac.uk/showcase/rectab.cgi?id=1061), and [GP date
coding 819](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=819).

The original lower-level clinical interface remains available:

```r
gp_raw <- migraine_query |>
  rap_plan_gp_query(database = database_name) |>
  rap_run_gp_query(
    connection = sc,
    output = "gp_migraine_records.csv"
  )

gp_cases <- gp_raw |>
  parse_gp_clinical() |>
  summarise_gp_diagnoses()
```


## AI Agent Skills for UKB data analyses

UKBAnalytica ships a curated set of AI agent skills under
`inst/skills/UKBAnalytica_skills/`. Each skill is a self-contained prompt
document that a Claude Code agent (or any compatible AI assistant) can load to
generate RAP-ready R scripts, plan analysis workflows, and interpret aggregate
outputs — without ever seeing real participant-level data.

| Skill | Phase | Coverage |
|---|---|---|
| `ukbsci-rap-extract` | P2 | Discover UKB fields and execute extractions via `dx extract_dataset` (sync) or table-exporter (async) |
| `ukbsci-cohort` | P2 | Define disease phenotypes from ICD-10/9, self-report, death, OPCS4, cancer registry, and build Cox-ready survival datasets |
| `ukbsci-workflow` | P2 | End-to-end study orchestrator — produces a phased plan and calls sub-skills in order |
| `ukbsci-regression` | P3 | Batch linear / logistic / Cox / GLM / negative-binomial / GAM regression; unified `run_regression()` interface; PH diagnostics; lag sensitivity; Fine-Gray competing risks; trend tests |
| `ukbsci-survival` | P3 | Kaplan-Meier curves with risk tables and log-rank p-values |
| `ukbsci-baseline` | P3 | Stratified Table 1 (baseline characteristics) via `tableone` |
| `ukbsci-propensity` | P4 | Propensity scores, PSM, IPTW (ATE/ATT/ATC), balance diagnostics, weighted regression |
| `ukbsci-mediation` | P4 | Causal mediation analysis (4-way decomposition); single and multi-mediator with sensitivity |
| `ukbsci-subgroup-sensitivity` | P4 | Subgroup × interaction tests across Cox / logistic / linear / GLM / negbin; complete-case and early-event sensitivity filters |
| `ukbsci-imputation` | P4 | Multiple imputation (mice) and Rubin's-rules pooling (mitools) with FMI diagnostics |
| `ukbsci-metabolomics` | P5 | Nightingale NMR metabolite ORA — name mapping, classification, custom or MetaboAnalystR backend, dot/bar plot visualization |
| `ukbsci-proteomics` | P5 | UKB Olink / UKB-PPP: ID mapping, GO/KEGG ORA, STRING PPI, community detection |
| `ukbsci-ml` | P5 | End-to-end ML workflows (classification + survival): split, feature select, tune, fit, evaluate, SHAP |
| `ukbsci-preprocess` | P5 | Variable cleaning, composite-variable builders (BP, air pollution, diet score), variable-set catalog |
| `ukbsci-plot` | P6 | Manuscript figures: forest, volcano, calibration; shared neutral theme/palettes; multi-format save helper |

### Data privacy boundary

All skills enforce a strict **script-generation-only** boundary:

- **What the agent receives:** column names, variable roles, intended analysis design, and aggregate outputs (flow counts, coefficient tables, model metrics, enrichment results, non-identifying figures).
- **What the agent must never receive:** real UKB participant rows, `eid` values, exact dates, raw RAP fields, row-level predictions, SHAP matrices, screenshots, or log excerpts containing row-level values.

The workflow is: describe your schema → the agent generates an R script → you run the script inside RAP → you share only aggregate results back for interpretation. Real participant-level data should remain inside the approved RAP project and RAP-controlled storage at all times.


## Reuse an R environment on RAP

RStudio workers on UKB-RAP are temporary, so packages installed only on the
worker are removed when the instance is terminated. To reuse an existing
project and its installed R packages, initialize an `renv` project and back up
the complete project directory to persistent RAP project storage.

After installing the required packages, create or update the dependency
snapshot:

```r
install.packages("renv")
renv::init()      # run once for a new project
install.packages("UKBAnalytica")
renv::snapshot()  # run again after changing package dependencies
```

In the RStudio Terminal, move into the project directory and create the
backup:

```bash
cd ~/ukbanalytica_project
dx-backup-folder -d /.Backups/ukbanalytica_environment.tar.gz
```

The archive is uploaded to the `.Backups` folder in the current RAP project,
not retained only on the temporary worker. It therefore remains available
after the RStudio instance is terminated. Its presence can be checked with:

```bash
dx ls -a /.Backups/
```

After starting a new RStudio instance, restore the project from RAP storage:

```bash
dx-restore-folder \
  /.Backups/ukbanalytica_environment.tar.gz \
  ~/ukbanalytica_project
```

Then activate the restored environment in R:

```r
setwd("~/ukbanalytica_project")
renv::activate()
renv::restore()
library(UKBAnalytica)
```

The backup should contain project-local R package libraries, configuration
files, and analysis scripts only. System-level dependencies installed with
tools such as `apt` are not captured reliably and should instead be packaged
in a Docker image. UK Biobank participant-level data should remain in approved
RAP-controlled data locations and should not be included in an environment
archive.

## Supplementary Materials
Here we provide some learning materials for UK Biobank in which you may be interested:
- [UK Biobank database browser](https://biobank.ndph.ox.ac.uk/ukb/index.cgi)
- [UK Biobank RAP platform](https://ukbiobank.dnanexus.com/landing)
- [UK Biobank GitHub resources](https://github.com/UK-Biobank)
- [UK Biobank RAP platform user guide](https://dnanexus.gitbook.io/uk-biobank-rap)
- [UK Biobank GP description](https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/primary_care_data.pdf)
- [UK Biobank learning guides supported by our team](https://hinna0818.github.io/Bioinfo-SMU/Epidemiology/UK_Biobank/) 

## Star History

<p align="center">
  <a href="https://repostars.dev/?repos=Hinna0818%2FUKBAnalytica&amp;theme=light">
    <picture>
      <source
        media="(prefers-color-scheme: dark)"
        srcset="https://repostars.dev/api/embed?repo=Hinna0818%2FUKBAnalytica&amp;theme=dark"
      />
      <img
        src="https://repostars.dev/api/embed?repo=Hinna0818%2FUKBAnalytica&amp;theme=light"
        width="700"
        alt="UKBAnalytica star history (estimated timeline, exact current total)"
      />
    </picture>
  </a>
</p>
