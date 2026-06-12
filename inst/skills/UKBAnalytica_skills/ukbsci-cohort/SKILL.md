---
name: ukbsci-cohort
description: >
  Build a UK Biobank disease cohort from a phenotype table extracted on the
  Research Analysis Platform (RAP), using the UKBAnalytica R package.
  Covers disease-definition construction and catalog lookup
  (get_predefined_diseases, get_disease_catalog, get_pomegranate_diseases,
  create_disease_definition, combine_disease_definitions), source parsers
  (parse_icd10_diagnoses, parse_icd9_diagnoses, parse_opcs4_procedures,
  parse_cancer_registry, parse_death_records, parse_self_reported_illnesses),
  multi-source case extraction (extract_cases_by_source, extract_disease_diagnosis,
  extract_disease_history, extract_disease_history_sensitivity,
  compare_data_sources), reusable follow-up time skeletons
  (ukb_time_skeleton), and the Cox-ready dataset builder
  build_survival_dataset() with prevalent vs incident separation and follow-up
  time computation, plus select_incident_by_years() for time-window
  stratification. Use this skill when the user asks to define a
  UKB disease phenotype, inspect curated/Pomegranate diagnostic codes,
  separate prevalent from incident cases, build a survival dataset, compare
  diagnostic sources, or compute follow-up time.
  Triggers: UKB cohort, UKB disease definition, prevalent vs incident,
  survival dataset, build_survival_dataset, ICD10 phenotyping, UKB 队列,
  病例定义, 生存数据集, /ukbsci-cohort. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-cohort — UK Biobank disease cohort & survival dataset builder

## 0. RAP guardrails (must always hold)

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

This skill assumes its input — the phenotype `data.table` — was extracted by
**`ukbsci-rap-extract`** (or equivalent `rap_*` calls) and lives in RAP
project storage. The same rules carry over:

- Inputs and outputs default to `/mnt/project/...`. Do not emit code that
  writes the cohort table or per-disease case files to a local laptop path.
- Aggregate results (cohort flow tables, Table 1, regression coefficients)
  may be shared with the local agent. Per-participant tables are not.
- If the user asks "can I just run this locally with a small sample?", answer
  with a *simulated* toy dataset — never a sample of real participants.

---

## 1. When to load this skill

- Define a UKB disease phenotype from scratch (`create_disease_definition`)
  or pick one of the pre-validated definitions (`get_predefined_diseases`).
- Combine multiple definitions into one composite endpoint
  (`combine_disease_definitions`).
- Identify cases from one or several diagnostic sources (ICD-10, ICD-9, OPCS-4,
  death registry, cancer registry, first-occurrence, self-report,
  UKB algorithmically-defined outcomes).
- Separate **prevalent** (baseline-or-earlier) from **incident** (post-baseline)
  cases and compute follow-up time for Cox / competing-risks models.
- Run a sensitivity analysis where the definition source set varies.
- Stratify incident cases by years since enrollment (e.g. ≤5y vs >5y).

## 2. When NOT to load this skill

- The user has not extracted phenotype data yet → load **`ukbsci-rap-extract`**.
- The user wants to fit a regression / Cox model on the survival dataset →
  load **`ukbsci-regression`** or **`ukbsci-survival`**.
- The user wants a Table 1 summary of the cohort → load **`ukbsci-baseline`**.
- The user wants propensity-score matching → load **`ukbsci-propensity`**.

---

## 3. Prerequisites

Required columns in the input `data.table`:

| Purpose | UKB field IDs / RAP columns |
|---------|------------------------------|
| Participant ID | `eid` |
| Baseline date | `p53_i0` (date of attending assessment centre, instance 0) |
| ICD-10 | `p41270` (concatenated codes) + `p41280_a*` (date array) |
| ICD-9 | `p41271` + `p41281_a*` |
| OPCS-4 procedures | `p41272` or `p41272_a*` + `p41282_a*` |
| Death | `p40000_i*` (date) + `p40001_i*` (primary cause) + `p40002_i*_a*` (secondary causes) |
| Cancer registry | `p40005_i*`, `p40006_i*`, `p40011_i*`, `p40012_i*` |
| Self-report | `p20002_i*_a*` (codes) + `p20008_i*_a*` (year as fractional decimal) |
| First Occurrence | UKB Category 1712 date fields (`p131xxx`) + source fields (`p131xxx + 1`) |
| Algorithmic outcomes | UKB Category 42 date field (e.g. `p42016` for COPD) + source field (`p42017`) |

The skill must inspect `dt` and warn the user about any source they request
that is not represented in the columns.

```r
library(UKBAnalytica)
library(data.table)

# dt is the data.table produced by rap_extract_pheno() / rap_submit_extract()
stopifnot(data.table::is.data.table(dt))
stopifnot("eid" %in% colnames(dt))
stopifnot("p53_i0" %in% colnames(dt))   # baseline date
```

---

## 4. Pipeline (4 phases)

### Phase 1 — Pick or build disease definitions

Start with the curated catalogue. In UKBAnalytica 1.0.0, the curated library
contains **86** analysis-ready disease definitions. The source-aware catalog
also exposes **313** Pomegranate-derived definitions. The expanded union of
curated and Pomegranate definitions currently contains **331** usable disease
keys; the conservative intersection contains **56** shared keys.

```r
defs_all <- get_predefined_diseases()
names(defs_all)
#> "AA", "TAA", "AAA", "CVD", "MI", "HF", "Stroke", "Hypertension",
#> "Diabetes", "T1DM", "T2DM", "Atrial_Fibrillation", "Lung_Cancer", ...

# Subset to your study endpoints
diseases <- defs_all[c("AA", "Hypertension", "Diabetes")]
```

Inspect exact code evidence before committing to an endpoint:

```r
# Code-level table: source, code_system, field_id, code, validation_status
get_disease_catalog(disease = "COPD")
get_disease_catalog(disease = "diabetes", code_system = "ICD-10")

# Pomegranate-derived disease definitions that current parsers can use
pom_defs <- get_predefined_diseases(source = "pomegranate")

# Matched curated + Pomegranate definitions.
# "intersection" is conservative; "union" expands matched code sets and
# retains unmatched definitions from both sources.
shared_defs <- get_predefined_diseases(
  source = "both",
  merge_type = "intersection"
)
expanded_defs <- get_predefined_diseases(
  source = "both",
  merge_type = "union"
)
```

Use `source = "curated"` for protocol-driven analyses where stable manual
definitions are preferred. Use `source = "both", merge_type = "union"` for
endpoint discovery or sensitivity analysis. When the exact source matters,
write the selected source and merge strategy into the analysis script.

Define a custom phenotype:

```r
copd_def <- create_disease_definition(
  name             = "COPD",
  icd10_pattern    = "^J4[0-4]",        # J40-J44
  icd9_pattern     = "^49[0-2]",        # 490-492
  sr_codes         = c(1112, 1113, 1472),
  death_icd10      = "^J4[0-4]",
  opcs4_pattern    = NULL,
  algo_date_field   = 42016,            # UKB algorithmic COPD date field
  algo_source_field = 42017             # Source: hospital / death / algorithm
)
```

Custom **cancer** phenotype using histology + behaviour codes:

```r
lung_ca_def <- create_disease_definition(
  name                 = "LungCancer_strict",
  cancer_icd10_pattern = "^C34",
  cancer_histology     = c(8070, 8071, 8140, 8240),   # squamous, adenoca, etc.
  cancer_behaviour     = 3L                            # malignant only
)
```

Combine definitions when you need composite endpoints:

```r
cvd_composite <- combine_disease_definitions(
  defs_all$MI, defs_all$Stroke, defs_all$HF,
  name = "MACE"
)
```

### Phase 2 — Parse source tables (optional inspection)

You normally don't call the parsers directly — `extract_cases_by_source()`
does it for you. But for QC / debugging it is useful to know the schemas:

```r
icd10 <- parse_icd10_diagnoses(dt)
# → eid, icd10_code, diag_date (Date), source = "ICD10"

icd9  <- parse_icd9_diagnoses(dt)
opcs4 <- parse_opcs4_procedures(dt)
canc  <- parse_cancer_registry(dt)
# → cancer_icd10_code, cancer_histology, cancer_behaviour, diag_date

death <- parse_death_records(dt)
# → death_code, death_date, source = "Death", cause_type (primary/secondary)

dod   <- get_death_dates(dt)
# → eid, death_date  (earliest per participant — for censoring)

sr    <- parse_self_reported_illnesses(dt, baseline_col = "p53_i0")
# → sr_code, diag_date (month-precision, imputed to day 15), instance, array_idx
```

If a parser returns a 0-row table, the corresponding source columns were not
present — warn the user before they request that source downstream.

### Phase 3 — Extract cases (one disease, one source, or many of each)

#### 3a. Single-disease diagnosis table — `extract_disease_diagnosis()`

Use this when the user wants a compact diagnosis table before building a full
survival dataset.

```r
mi_def <- get_predefined_diseases(
  source = "both",
  merge_type = "union",
  disease = "MI"
)

mi_dx <- extract_disease_diagnosis(
  dt                  = dt,
  disease             = "MI",
  disease_definitions = mi_def,
  sources             = c("ICD10", "ICD9", "Death"),
  censor_date         = as.Date("2023-10-31"),
  baseline_col        = "p53_i0",
  include_all         = TRUE
)
# Columns: eid, disease, diagnosed, prevalent_case, incident_case,
# earliest_date, diagnosis_source, status, surv_time
```

#### 3b. The flexible primitive — `extract_cases_by_source()`

```r
cases_long <- extract_cases_by_source(
  dt                  = dt,
  disease_definitions = diseases,
  sources             = c("ICD10", "ICD9", "Self-report", "Death",
                          "OPCS4", "CancerRegistry",
                          "FirstOccurrence", "Algorithm"),
  censor_date  = as.Date("2023-10-31"),
  baseline_col = "p53_i0"
)
str(cases_long)
#> eid, disease, status, surv_time, prevalent_case,
#> earliest_date, diagnosis_source
```

`diagnosis_source` is the *winning* source per participant–disease (whichever
gave the earliest valid date). For `"Algorithm"`, it is further refined to
`"Algorithm_<source-code>"`.

#### 3c. The history-only helper — `extract_disease_history()`

For "did the participant ever have disease X (per the chosen sources)?":

```r
history <- extract_disease_history(
  dt                  = dt,
  diseases            = c("Hypertension", "Diabetes"),
  disease_definitions = defs_all,
  sources             = c("ICD10", "ICD9", "Self-report"),
  baseline_col        = "p53_i0"
)
#> Adds columns: Hypertension_history (0/1), Diabetes_history (0/1)
```

#### 3d. Sensitivity variants — `extract_disease_history_sensitivity()`

For pre-baked sensitivity strata (ICD-10-only, ICD-10+ICD-9, all sources):

```r
sens <- extract_disease_history_sensitivity(
  dt                  = dt,
  diseases            = c("Diabetes"),
  disease_definitions = defs_all,
  baseline_col        = "p53_i0"
)
#> Adds: Diabetes_history_ICD10, Diabetes_history_hospital, Diabetes_history_all
```

#### 3e. Source-by-source comparison — `compare_data_sources()`

QC tool. Returns a per-disease table of how many cases each source uncovers
and how much they overlap:

```r
qc <- compare_data_sources(
  dt                  = dt,
  disease_definitions = diseases,
  baseline_col        = "p53_i0"
)
```

### Phase 4 — Build reusable time skeleton, then Cox-ready dataset

For prospective analyses, first standardize the time backbone once and reuse
it across disease endpoints.

```r
time_skel <- ukb_time_skeleton(
  data                 = dt,
  id_col               = "eid",
  baseline_col         = "p53_i0",
  birth_year_col       = "p34",
  birth_month_col      = "p52",
  age_col              = "p21022",
  death_date_cols      = "^(participant\\.)?p40000_i[0-9]+$",
  lost_to_followup_col = "p191",
  admin_censor_date    = as.Date("2023-10-31")
)
```

`admin_censor_date` remains the database-level administrative censoring date.
`followup_end_date` is participant-specific and uses the earliest of death,
loss-to-follow-up, and administrative censoring.

The main entry point. Most user requests should resolve here.

```r
cohort <- build_survival_dataset(
  dt                  = dt,
  disease_definitions = diseases,
  prevalent_sources   = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources     = c("ICD10", "ICD9", "Death"),
  censor_date         = as.Date("2023-10-31"),
  baseline_col        = "p53_i0",
  time_skeleton       = time_skel,
  primary_disease     = "AA",          # which disease drives outcome_*
  output              = "wide",        # or "long"
  include_all         = TRUE,
  show_flow           = TRUE,
  dt_threads          = 8
)

# Inspect attrition
flow <- attr(cohort, "participant_flow")
if (!is.null(flow)) print(flow)

# Check created columns without printing participant rows
stopifnot(all(c("AA_history", "AA_incident",
                "outcome_status", "outcome_surv_time") %in% names(cohort)))
cohort[, .(
  n = .N,
  n_incident = sum(outcome_status == 1, na.rm = TRUE),
  n_censored = sum(outcome_status == 0, na.rm = TRUE),
  n_prevalent = sum(is.na(outcome_status))
)]
```

Column semantics added by `build_survival_dataset()`:

| Column | Meaning |
|--------|---------|
| `<Disease>_history` | 1 = prevalent case (≤ baseline from `prevalent_sources`); 0 otherwise |
| `<Disease>_incident` | 1 = incident case (> baseline from `outcome_sources`); 0 otherwise |
| `outcome_status` | 1 = incident `primary_disease` event; 0 = censored; **NA** = prevalent (not at risk) |
| `outcome_surv_time` | Follow-up time in years; **NA** for prevalent |

**Prevalent ≠ at risk.** Participants with `<primary>_history = 1` have
`outcome_status = NA` and `outcome_surv_time = NA` so that downstream Cox
fits drop them automatically.

### Phase 5 — Optional: stratify incident cases by years since enrollment

For lag analyses ("MI in the first 5 years vs after"):

```r
strat <- select_incident_by_years(
  df             = cohort,
  n_years        = 5,
  time_col       = "outcome_surv_time",
  status_col     = "outcome_status",
  baseline_col   = "p53_i0",
  event_date_col = "earliest_date",
  group_col      = "incident_timing",
  output         = "split",          # or "combined"
  copy           = TRUE,
  verbose        = TRUE
)
str(strat)
#> List of 2: $ within_5_years, $ after_5_years
```

`output = "combined"` returns one data.table with a `incident_timing` column
taking values `"within_5_years"` / `"after_5_years"`.

---

## 5. Decision rules the agent should ask the user about

Pause and ask before generating code if any of these is ambiguous:

1. **Which disease(s)?** Check `get_disease_catalog(disease = "<term>")`
   and `names(get_predefined_diseases())` first. Only build custom definitions
   if the requested endpoint is absent or the protocol requires different
   source/code rules.
2. **`primary_disease`?** Required for `outcome_status` / `outcome_surv_time`.
   Default to the first disease, but confirm.
3. **`prevalent_sources` vs `outcome_sources`?** Defaults are:
   - prevalent: `c("ICD10", "ICD9", "Self-report", "Death")` — broad,
     because we want to *exclude* people who already have the disease
   - outcome: `c("ICD10", "ICD9", "Death")` — Self-report excluded because
     its year-only dates corrupt follow-up time
   Confirm whether the user wants to add `OPCS4`, `CancerRegistry`,
   `FirstOccurrence`, or `Algorithm` to outcomes.
4. **`censor_date`?** Default `2023-10-31`. Should match the UKB data release
   the user has access to — confirm with them.
5. **`output = "wide"` vs `"long"`?** Wide is the right default; long is for
   per-source / per-disease audits.
6. **`show_flow = TRUE`?** Default on; the attrition table is **safe to
   export** because it contains only aggregate counts.

---

## 6. Key functions (cheat-sheet)

| Function | Purpose | Returns |
|----------|---------|---------|
| `get_predefined_diseases()` | Catalogue of curated disease definitions | named list |
| `get_disease_catalog(source, disease, code_system)` | Source-aware code evidence table | data.frame |
| `get_pomegranate_diseases(disease)` | Convert Pomegranate catalog rows to disease definitions | named list |
| `get_pomegranate_source_manifest()` | Provenance for Pomegranate YAML and portal audit sources | data.frame |
| `load_pomegranate_portal_coding(path)` | Load an optional local portal CSV audit table | data.frame |
| `create_disease_definition()` | Build a custom definition | list |
| `combine_disease_definitions(..., name)` | Composite endpoint (union) | list |
| `parse_icd10_diagnoses(dt)` | ICD-10 long table | data.table |
| `parse_icd9_diagnoses(dt)` | ICD-9 long table | data.table |
| `parse_opcs4_procedures(dt)` | OPCS-4 long table | data.table |
| `parse_cancer_registry(dt)` | Cancer registry long table | data.table |
| `parse_death_records(dt)` | Death registry long table | data.table |
| `get_death_dates(dt)` | Earliest death date per eid | data.table |
| `parse_self_reported_illnesses(dt, baseline_col)` | Self-report long table | data.table |
| `extract_cases_by_source(dt, defs, sources, ...)` | Flexible multi-source extractor | data.table |
| `extract_disease_history(dt, diseases, defs, sources, baseline_col)` | Boolean history flags | data.table |
| `extract_disease_history_sensitivity(dt, diseases, defs, baseline_col)` | Pre-baked sensitivity flags | data.table |
| `compare_data_sources(dt, defs, baseline_col)` | Source-overlap QC | data.table |
| `extract_disease_diagnosis(dt, disease, disease_definitions, sources, ...)` | Diagnosis table for selected diseases | data.table |
| `ukb_time_skeleton(data, id_col, baseline_col, ..., admin_censor_date)` | Reusable time backbone | data.table |
| `build_survival_dataset(dt, defs, prevalent_sources, outcome_sources, ..., primary_disease, ..., show_flow, dt_threads)` | **Main** Cox-ready dataset | data.table (+ `attr "participant_flow"`) |
| `select_incident_by_years(df, n_years, ...)` | Stratify by follow-up years | data.table or list |

See [`references/functions.md`](references/functions.md) for full argument
tables and return-value details.

---

## 7. Common pitfalls

1. **`primary_disease` not in `disease_definitions`.** Function errors out.
   Always pass a name from `names(disease_definitions)`.
2. **Sources requested but columns missing.** If `dt` lacks `p41270`,
   requesting `"ICD10"` returns 0 cases silently and the user gets a
   misleading "no events" cohort. Inspect `parse_icd10_diagnoses(dt)`
   first when in doubt.
3. **Self-report in outcomes.** Self-report has only fractional-year
   precision; including it in `outcome_sources` corrupts follow-up time.
   Default leaves it out; only override if you fully understand the cost.
4. **Cancer histology codes.** UKB stores them in `p40011_i*`. When using
   `cancer_histology`, codes must match exactly (UKB uses the same numeric
   codes as the WHO ICD-O-3 dictionary).
5. **First Occurrence special dates.** UKB encodes "before 1900", "1901-01-01",
   "2037-07-07", etc. as sentinel values; the parser drops them automatically
   — do not re-add them.
6. **Date-format chaos.** UKB columns may be character, R Date, or compact
   `YYYYMMDD` numerics. `extract_cases_by_source()` routes everything
   through `.safe_as_date()`. Trust the function but read its warnings.
7. **Censoring at wrong date.** Default `2023-10-31` corresponds to one
   specific UKB release. If the user is on a different release, override
   `censor_date`.
8. **Combining `Algorithm` with hospital sources.** UKB's algorithmic
   outcomes already incorporate hospital and death data. Adding `ICD10` *and*
   `Algorithm` to `outcome_sources` is usually redundant; pick one. Use
   `compare_data_sources()` to quantify the overlap.
9. **`show_flow = TRUE` but `output = "long"`.** The attrition table only
   makes sense for wide format. Confirm before mixing.

---

## 8. Related skills

| Skill | When to switch |
|-------|----------------|
| **`ukbsci-rap-extract`** | When the user has no phenotype data yet, or needs additional fields. |
| **`ukbsci-preprocess`** | After cohort build; clean missing codes, derive composite covariates. |
| **`ukbsci-baseline`** | After cohort build; produce Table 1 stratified by `<primary>_history` or exposure. |
| **`ukbsci-regression`** | Fit Cox / logistic / linear models on `outcome_status` / `outcome_surv_time`. |
| **`ukbsci-survival`** | KM curves, competing risks, time-dependent ROC. |
| **`ukbsci-subgroup-sensitivity`** | Run sensitivity analyses across source variants. |
| **`ukbsci-workflow`** | End-to-end orchestrator; calls this skill as Phase 2. |

---

## 9. References

- [`references/functions.md`](references/functions.md) — every function
  signature, argument table, return-value details
- [`references/disease-catalog.md`](references/disease-catalog.md) — the
  pre-validated phenotype catalogue from `get_predefined_diseases()`
- [`references/source-matrix.md`](references/source-matrix.md) — which
  source needs which UKB field IDs; recommended combinations
- [`references/examples.md`](references/examples.md) — three minimal
  end-to-end examples
