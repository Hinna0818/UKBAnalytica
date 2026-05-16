---
name: ukbsci-preprocess
description: >
  Variable cleaning and feature construction for UK Biobank Research
  Analysis Platform (RAP) cohorts using UKBAnalytica. Negative-code
  sanitation (ukb_clean_missing), unified preprocessing with predefined
  variable sets (preprocess_baseline, prepare_analysis_dataset), composite
  variable builders (calculate_blood_pressure, calculate_air_pollution,
  calculate_diet_score), and curated variable-set lookups
  (get_variable_set, get_variable_sets, get_variable_info). Use this skill
  when the user asks to clean UKB negative-codes (-1, -3, -7, …), derive
  blood pressure / air pollution / diet score, compose an analysis-ready
  table, or look up predefined UKB variable sets. Triggers: UKB
  preprocessing, variable cleaning, derive BP, air pollution exposure, diet
  score, negative code, UKB 变量预处理, 缺失码处理, 复合变量,
  /ukbsci-preprocess. Hard rule: cleaned participant-level tables stay inside
  RAP project storage; aggregate missingness/QC summaries and de-identified
  figures can be exported.
---

# ukbsci-preprocess — UKB variable cleaning & composite-feature builders

## 0. RAP guardrails

Shared privacy boundary: do not export UK Biobank RAP individual-level raw
data, direct identifiers (`eid`), exact dates, raw RAP fields, or row-level
source tables that can be linked back to participants. De-identified analytical
figures and aggregate summaries (curves, coefficients, metrics, feature-level
or bin-level source tables) are generally exportable when no identifying or raw
participant-level fields accompany them.

The cleaned cohort is participant-level. Stay inside `/mnt/project/...`.
Summary statistics, missingness counts, and pre-vs-post sanity tables are
aggregate and safe to export.

---

## 1. When to load

- Convert UKB negative sentinel codes (-1, -2, -3, -7, -9, -10, -11, -13,
  -15) to `NA`.
- Build composite variables: averaged blood pressure, averaged air pollution,
  diet score.
- Apply a predefined variable set (e.g. `clinical_core`) to a freshly
  extracted RAP table.
- Wrap up preprocessing + optional imputation in one call
  (`prepare_analysis_dataset`).

## 2. When NOT to load

- Pulling fields off RAP → `ukbsci-rap-extract`.
- Multiple imputation (more than the optional flag here) →
  `ukbsci-imputation`.
- Defining disease cases / survival outcome → `ukbsci-cohort`.

---

## 3. Prerequisites

```r
library(UKBAnalytica); library(data.table)
stopifnot(data.table::is.data.table(dt))
```

---

## 4. Pipeline

### Phase 1 — Browse predefined variables

```r
get_variable_info(category = "all")         # 9 categories supported
get_variable_info(category = "biomarkers")
get_variable_sets(category = "lifestyle")    # all sets, with field IDs
get_variable_set("clinical_core", output = "field_id")  # int vector
```

The catalogues cover: `demographics`, `anthropometrics`, `lifestyle`,
`socioeconomic`, `blood_pressure`, `medications`, `biomarkers`, `pollution`,
`diet`.

### Phase 2 — Sanitize negative codes

```r
dt <- ukb_clean_missing(
  dt,
  cols          = NULL,              # NULL = scan all numeric cols
  invalid_codes = c(-1, -3, -7, -9, -11, -13, -15)
)
```

This is **the first** transformation. Run it before computing composite
variables — otherwise the means / scores are silently corrupted by
sentinels.

### Phase 3 — Build composite variables

```r
# Mean SBP / DBP across two readings; falls back to manual if automated missing
dt <- calculate_blood_pressure(dt, type = "sbp", prefer = "auto")
dt <- calculate_blood_pressure(dt, type = "dbp", prefer = "auto")

# Average air pollution across time points
dt <- calculate_air_pollution(dt, pollutants = c("NO2","PM10","PM2.5","NOx"))

# Healthy diet score 0–7
dt <- calculate_diet_score(
  dt,
  components  = c("fruit","vegetable","fish","meat","cereal","milk"),
  na_handling = "strict"           # "strict" = NA in any component → NA score
)
```

### Phase 4 — Unified call (`preprocess_baseline`)

For routine workflows:

```r
dt <- preprocess_baseline(
  dt,
  variables       = c("age","sex","bmi","sbp","dbp","ldl_cholesterol",
                       "hba1c","crp","smoking_status","alcohol_intake"),
  custom_mapping  = NULL,
  missing_action  = "keep",        # "drop" filters rows with any NA
  invalid_codes   = c(-1, -3)
)
```

`preprocess_baseline()` understands the same variable names as
`get_variable_info()`; passing a name not in the catalogue raises an error.

### Phase 5 — End-to-end with optional imputation

```r
ready <- prepare_analysis_dataset(
  dt,
  preprocess_vars = c("age","sex","bmi","sbp","ldl_cholesterol","hba1c"),
  impute          = FALSE,           # TRUE to call mice via run_imputation
  impute_method   = "pmm",
  impute_m        = 5,
  seed            = 1234
)
```

If `impute = TRUE`, the return value is the list-of-imputations from
`run_imputation()`. For full MI workflows use `ukbsci-imputation`.

---

## 5. Common pitfalls

1. **Order matters.** Always `ukb_clean_missing()` before any composite
   variable; otherwise sentinel codes contaminate the means.
2. **`calculate_diet_score()` rounds to one decimal place** in some
   versions. If you need higher precision for downstream Cox, scale the
   raw components yourself.
3. **`prefer = "manual"`** in `calculate_blood_pressure()` swaps the
   precedence — handy if your cohort has more manual readings than
   automated.
4. **`na_handling = "partial"`** in `calculate_diet_score()` averages
   available components and scales — useful but introduces a different
   missingness mechanism. Document it explicitly.
5. **`missing_action = "drop"`** in `preprocess_baseline()` removes any
   row with any NA in the listed `variables`. Confirm this is what the
   user wants, especially for high-dimensional feature sets.
6. **Custom mappings.** `custom_mapping = list(my_alias = "p123_i0")`
   lets you rename without touching the catalogue. Useful for non-standard
   fields (e.g. derived imaging metrics).
7. **Predefined variables can drift.** When UKB releases new fields, the
   catalogue is updated. Run `get_variable_info()` interactively before
   freezing scripts.

---

## 6. Key functions

| Function | Purpose |
|----------|---------|
| `ukb_clean_missing(df, cols = NULL, invalid_codes = c(-1,-3,-7,-9,-11,-13,-15))` | Replace sentinel codes with NA |
| `calculate_blood_pressure(df, type = c("sbp","dbp"), prefer = c("auto","manual"))` | Mean BP across readings |
| `calculate_air_pollution(df, pollutants = c("NO2","PM10","PM2.5","NOx"))` | Time-averaged exposures |
| `calculate_diet_score(df, components = c("fruit","vegetable","fish","meat","cereal","milk"), na_handling = c("strict","partial"))` | 0–7 healthy diet score |
| `preprocess_baseline(df, variables, custom_mapping = NULL, missing_action = c("keep","drop"), invalid_codes = c(-1,-3))` | Unified pipeline using predefined variable map |
| `prepare_analysis_dataset(df, preprocess_vars = NULL, impute = FALSE, impute_method = "pmm", impute_m = 5, seed = NULL)` | Preprocessing + optional MI |
| `get_variable_info(category = "all")` | Catalogue rows |
| `get_variable_set(set, output = c("data.frame","field_id","ukb_col"))` | One set |
| `get_variable_sets(set = NULL, category = NULL)` | All sets with filter |

See [`references/variable-catalog.md`](references/variable-catalog.md) for
the 9 categories and example variables.

---

## 7. Related skills

| Skill | When |
|-------|------|
| `ukbsci-rap-extract` | Get the raw fields first. |
| `ukbsci-cohort` | Build outcomes from the cleaned table. |
| `ukbsci-baseline` | Table 1 of the cleaned cohort. |
| `ukbsci-imputation` | Full MI pooling (vs. the one-shot option here). |

---

## 8. References

- [`references/functions.md`](references/functions.md)
- [`references/variable-catalog.md`](references/variable-catalog.md) — the
  9 categories and what's in each
- [`references/examples.md`](references/examples.md)
