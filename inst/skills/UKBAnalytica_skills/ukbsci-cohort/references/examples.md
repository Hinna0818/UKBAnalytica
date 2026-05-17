# `ukbsci-cohort` — examples

Three minimal examples. All assume `dt` is the phenotype `data.table`
extracted by `ukbsci-rap-extract` and **already resident in RAP project storage**.

Standard header for every example:

```r
###############################################################################
# UKBAnalytica skill: ukbsci-cohort
# Citation: He N. UKBAnalytica R package. https://github.com/Hinna0818/UKBAnalytica_v2
###############################################################################

library(UKBAnalytica)
library(data.table)

stopifnot(data.table::is.data.table(dt))
stopifnot(all(c("eid", "p53_i0") %in% colnames(dt)))
```

---

## Example A — Aortic aneurysm cohort with adjustment exposures

A common pattern: one **primary disease** (AA) and several **comorbidities**
that will be used as covariates downstream.

```r
defs <- get_predefined_diseases()[
  c("AA", "Hypertension", "Diabetes", "CVD")
]

cohort <- build_survival_dataset(
  dt                  = dt,
  disease_definitions = defs,
  prevalent_sources   = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources     = c("ICD10", "ICD9", "Death"),
  censor_date         = as.Date("2023-10-31"),
  baseline_col        = "p53_i0",
  primary_disease     = "AA",
  output              = "wide",
  show_flow           = TRUE,
  dt_threads          = 8
)

# Inspect attrition (aggregate counts only; shareable with the local agent)
flow <- attr(cohort, "participant_flow")
print(flow)

# QC: how many incident events
cat("Incident AA events:",
    sum(cohort$outcome_status == 1, na.rm = TRUE), "\n")
cat("Censored:",
    sum(cohort$outcome_status == 0, na.rm = TRUE), "\n")
cat("Prevalent (NA outcome):",
    sum(is.na(cohort$outcome_status)), "\n")

# Save aggregate flow for the analysis report
fwrite(flow, "/mnt/project/results/aa_cohort_flow.csv")
```

---

## Example B — Custom COPD definition combining ICD + algorithm + self-report

When a curated definition is not exactly what the user wants, build a custom
one. This example combines hospital codes, UKB's algorithmic COPD endpoint,
and self-report.

```r
copd_def <- create_disease_definition(
  name              = "COPD_custom",
  icd10_pattern     = "^J4[0-4]",
  icd9_pattern      = "^49[0-2]",
  sr_codes          = c(1112, 1113, 1472),
  death_icd10       = "^J4[0-4]",
  algo_date_field   = 42016,
  algo_source_field = 42017
)

defs <- list(COPD = copd_def,
             Diabetes = get_predefined_diseases()$Diabetes)

# Before locking in source choices, see the overlap
qc <- compare_data_sources(
  dt                  = dt,
  disease_definitions = defs,
  baseline_col        = "p53_i0"
)
print(qc)

# Build the cohort: use Algorithm for outcome (UKB-validated dates)
cohort <- build_survival_dataset(
  dt                  = dt,
  disease_definitions = defs,
  prevalent_sources   = c("ICD10", "ICD9", "Self-report", "Death", "Algorithm"),
  outcome_sources     = c("Algorithm", "Death"),
  censor_date         = as.Date("2023-10-31"),
  baseline_col        = "p53_i0",
  primary_disease     = "COPD",
  show_flow           = TRUE
)

# Stratify incident cases by 5-year follow-up windows
strat <- select_incident_by_years(
  df             = cohort,
  n_years        = 5,
  status_col     = "outcome_status",
  time_col       = "outcome_surv_time",
  output         = "split",
  copy           = TRUE,
  verbose        = TRUE
)
str(strat)
#> List of 2: $ within_5_years, $ after_5_years
```

---

## Example C — Sensitivity comparison across source sets

For sensitivity analyses: prevalent-case flags built three different ways.

```r
defs <- get_predefined_diseases()[c("Hypertension", "Diabetes")]

sens <- extract_disease_history_sensitivity(
  dt                  = dt,
  diseases            = names(defs),
  disease_definitions = defs,
  baseline_col        = "p53_i0"
)

# Sens now has, per disease:
#   <Disease>_history_ICD10      (strict)
#   <Disease>_history_hospital   (ICD-10 + ICD-9)
#   <Disease>_history_all        (+ Self-report)
sens[, .N, by = .(
  Hypertension_history_ICD10,
  Hypertension_history_hospital,
  Hypertension_history_all
)]
```

These columns become the exposure / strata variables for
`ukbsci-subgroup-sensitivity`.

---

## Example D — Composite endpoint (MACE)

```r
defs_all <- get_predefined_diseases()
mace <- combine_disease_definitions(
  defs_all$MI, defs_all$Stroke, defs_all$HF,
  name = "MACE"
)

cohort <- build_survival_dataset(
  dt                  = dt,
  disease_definitions = list(MACE = mace,
                             Diabetes = defs_all$Diabetes),
  prevalent_sources   = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources     = c("ICD10", "ICD9", "Death"),
  primary_disease     = "MACE",
  show_flow           = TRUE
)
```

The cohort has columns `MACE_history`, `MACE_incident`, `Diabetes_history`,
`Diabetes_incident`, and `outcome_status` / `outcome_surv_time` driven by
MACE.

---

## Snippet — pre-flight column audit

Use before `build_survival_dataset()` to fail fast on missing source columns:

```r
need_for_source <- list(
  ICD10           = c("p41270", "p41280_a0"),
  ICD9            = c("p41271", "p41281_a0"),
  OPCS4           = c("p41272", "p41282_a0"),
  Death           = c("p40000_i0", "p40001_i0"),
  CancerRegistry  = c("p40005_i0", "p40006_i0",
                      "p40011_i0", "p40012_i0"),
  `Self-report`   = c("p20002_i0_a0", "p20008_i0_a0")
)

audit_sources <- function(dt, sources) {
  cols <- colnames(dt)
  problems <- list()
  for (s in sources) {
    needs <- need_for_source[[s]]
    if (is.null(needs)) next
    pat <- sub("_a0$", "_a", needs)
    miss <- !sapply(pat, function(p) any(startsWith(cols, p)))
    if (any(miss)) problems[[s]] <- needs[miss]
  }
  problems
}

audit <- audit_sources(dt, c("ICD10", "ICD9", "Death", "Self-report"))
if (length(audit) > 0) {
  message("Missing columns for sources: ",
          paste(names(audit), collapse = ", "))
  print(audit)
  stop("Re-run ukbsci-rap-extract to pull these fields before continuing.")
}
```
