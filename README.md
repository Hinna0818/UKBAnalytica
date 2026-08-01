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

## GWAS and PheWAS phenotype preparation

The GWAS/PheWAS module accepts outputs from the existing phenotyping workflow
without modifying them. It builds REGENIE-compatible phenotype/covariate files
and the four-column ICD-10 table used by the PheWAS R package. The builder is
in-memory only; writing participant-level files is a separate, explicit RAP
operation.

```r
icd10_long <- parse_icd10_diagnoses(ukb_raw)

gwas_phewas <- build_gwas_phewas_phenotype(
  data = surv_dt,
  phenotype_cols = "outcome_status",
  diagnoses = icd10_long,
  covariates = c("p31", "p21022"),
  categorical_covariates = "p31",
  diagnosis_date_col = "diag_date",
  sex_col = "p31",
  sample_qc = ukb_gwas_qc_profile("none")
)

gwas_phewas$trait_summary
gwas_phewas$qc_flow

# Inside RAP only; defaults to requiring a path below /mnt/project.
# files <- ukb_write_gwas_phewas_phenotype(
#   gwas_phewas,
#   "/mnt/project/Analysis/GWAS_PheWAS/phenotype"
# )
```

Use `ukb_plan_regenie()` and `ukb_plan_plink2_dosage()` to create auditable
dry-run command plans. Execution remains opt-in. See the
[GWAS and PheWAS workflow](https://hinna0818.github.io/UKBAnalytica/14-gwas-phewas.html)
for the complete RAP example.

The structured REGENIE planner accepts additional non-conflicting official
tokens through `step1_args` and `step2_args`. For complete control, use
`ukb_plan_regenie_command()`:

```r
custom_regenie <- ukb_plan_regenie_command(
  args = c(
    "--step", "2",
    "--pgen", "/mnt/project/Genotype/cohort",
    "--phenoFile", "/mnt/project/Analysis/phenotype.tsv",
    "--phenoColList", "disease",
    "--bt",
    "--pred", "/mnt/project/Analysis/step1_pred.list",
    "--bsize", "400",
    "--minMAC", "20",
    "--minINFO", "0.8",
    "--chr", "19",
    "--gz",
    "--out", "/mnt/project/Analysis/disease_chr19"
  ),
  label = "disease_chr19"
)

ukb_run_regenie(custom_regenie)  # dry run
```

Raw tokens are preserved, but their validity depends on the REGENIE version
installed in the execution environment.

Common genotype formats can be converted through the same PLINK2 execution
boundary:

```r
# BED/BIM/FAM -> BGEN 1.2; creates a plan and does not execute by default.
conversion <- convert_gwas_datatype(
  input = "/mnt/project/Genotype/array/ukb_array.bed",
  from = "auto",
  to = "bgen",
  output_prefix = "/mnt/project/Analysis/converted/ukb_array",
  bgen_bits = 8,
  plink2_args = c("--maf", "0.01", "--threads", "16")
)

conversion$commands$convert$display

# Explicit execution in a RAP environment with PLINK2 on PATH:
# conversion_run <- convert_gwas_datatype(
#   input = "/mnt/project/Genotype/array/ukb_array.bed",
#   to = "bgen",
#   output_prefix = "/mnt/project/Analysis/converted/ukb_array",
#   execute = TRUE
# )
```

For PLINK2 functionality without a dedicated wrapper, pass the official
command-line tokens to `ukb_plan_plink2()`, then execute the returned plan with
`ukb_run_plink2(..., execute = TRUE)`.

## Polygenic risk scores

Published UK Biobank PRS Release V2 fields can be inspected and loaded without
reading genotype files:

```r
ukb_prs_catalog("standard", c("HT", "T2D"))

# From a participant table already extracted in RAP:
published_prs <- load_ukb_prs(
  participant_data,
  trait = c("HT", "T2D"),
  type = "standard",
  standardize = TRUE
)

# Enhanced PRS defaults to setting scores outside field 26200 to NA.
# Use eligibility = "require" to retain testing-subgroup participants only.

# Or omit `data` to use the existing RAP phenotype extractor:
# published_prs <- load_ukb_prs(
#   trait = c("HT", "T2D"), type = "standard"
# )
```

For a custom fixed-weight PRS, first normalize a public weight table, then
create and explicitly execute a PLINK2 plan. The package passes genotype paths
to PLINK2 and never loads or caches UKB BGEN data in R.

```r
weights <- prepare_prs_weights(
  published_weights,
  variant_col = "rsid",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  chr_col = "chr",
  pos_col = "position",
  effect_allele_freq_col = "effect_allele_frequency",
  weight_col = "odds_ratio",
  weight_type = "or",       # converted to log(OR)
  genome_build = "GRCh37"
)

# ukb_variant_map is a target-map table with ID/CHR/POS/ALT/REF/ALT_FREQ.
# Builds must match; no silent liftOver is performed.
weights <- harmonize_prs_weights(
  weights,
  variant_map = ukb_variant_map,
  target_allele1_freq_col = "ALT_FREQ",
  target_build = "GRCh37",
  palindromic = "infer_by_af",
  output = "/mnt/project/Analysis/PRS/t2d_weights_harmonized.tsv"
)

attr(weights, "harmonization_qc")

prs_plan <- ukb_plan_prs(
  weights = weights,
  bgen = sprintf(
    "/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb_c%d_b0_v3.bgen",
    1:22
  ),
  sample = "/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3_s487395.sample",
  keep_file = "/mnt/project/Analysis/PRS/cohort.keep",
  frequency_file = "/mnt/project/Analysis/PRS/ukb_reference.acount",
  output_prefix = "/mnt/project/Analysis/PRS/t2d",
  score_name = "T2D_PRS",
  plink2_args = c("--threads", "16")
)

ukb_run_prs(prs_plan) # dry run
# result <- ukb_run_prs(prs_plan, execute = TRUE)
# result$scores
```

When several scores use the same target genotype map, combine their
harmonized weights before planning. This produces one PLINK2 genotype scan per
chromosome, with one score column per PRS. Splitting the weights also prevents
commands from being created for chromosomes with no contributing variants.

```r
multi_weights <- combine_prs_weights(
  list(HT_PRS = ht_weights, T2D_PRS = t2d_weights),
  output = "/mnt/project/Analysis/PRS/ht_t2d_weights.tsv"
)

weights_by_chr <- split_prs_weights(
  multi_weights,
  output_dir = "/mnt/project/Analysis/PRS/weights_by_chr"
)

multi_plan <- ukb_plan_prs(
  weights = weights_by_chr,
  genotype = sprintf("/mnt/project/Genotype/pgen/ukb_chr%d.pgen", 1:22),
  genotype_format = "pgen", # BED/BIM/FAM is also accepted with "bed"
  keep_file = "/mnt/project/Analysis/PRS/cohort.keep",
  output_prefix = "/mnt/project/Analysis/PRS/ht_t2d"
)

# Reruns only incomplete chromosomes; up to four PLINK2 processes at once.
# multi_result <- ukb_run_prs(
#   multi_plan, execute = TRUE, workers = 4, resume = TRUE
# )
# multi_result$scores[, .(IID, HT_PRS, T2D_PRS)]
```

`workers` controls independent chromosome-level PLINK2 processes. Keep the
sum of PLINK2 threads across workers within the RAP instance CPU allocation.
With `resume = TRUE`, complete `.sscore` (and `.sscore.vars`, when requested)
outputs are reused; partial outputs require `overwrite = TRUE` before rerun.

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
  <a href="https://www.star-history.com/?repos=Hinna0818%2FUKBAnalytica&type=date&legend=top-left">
    <img 
      src="https://api.star-history.com/svg?repos=Hinna0818/UKBAnalytica&type=date&legend=top-left"
      width="600"
      alt="Star History Chart"
    />
  </a>
</p>
