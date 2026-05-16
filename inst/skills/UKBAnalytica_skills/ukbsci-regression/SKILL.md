---
name: ukbsci-regression
description: >
  Fit linear, logistic, and Cox regression models (including extensions for
  proportional-hazards diagnostics, lagged sensitivity, competing risks,
  and trend tests) on UK Biobank Research Analysis Platform (RAP) cohorts
  built with UKBAnalytica. Wraps runmulti_lm, runmulti_logit, runmulti_cox,
  runmulti_cox_lag, runmulti_cox_zph, runmulti_competing, runmulti_trend, and
  ukb_cox_diagnostics. Use this skill when the user asks to fit Cox /
  logistic / linear models on the cohort, run batch regressions across
  multiple exposures, check proportional-hazards assumptions, run a lagged
  sensitivity analysis, perform competing-risks Fine-Gray modeling, or test a
  dose-response trend across ordered exposure categories. Triggers: UKB
  regression, Cox model, logistic regression, batch regression, proportional
  hazards test, competing risks, Fine-Gray, p_trend, UKB 回归, Cox 模型,
  /ukbsci-regression. Hard rule: cohort data with participant-level rows stays
  inside RAP project storage; aggregate coefficient tables and de-identified
  diagnostic figures can be exported.
---

# ukbsci-regression — UK Biobank batch regression

## 0. RAP guardrails

Shared privacy boundary: do not export UK Biobank RAP individual-level raw
data, direct identifiers (`eid`), exact dates, raw RAP fields, or row-level
source tables that can be linked back to participants. De-identified analytical
figures and aggregate summaries (curves, coefficients, metrics, feature-level
or bin-level source tables) are generally exportable when no identifying or raw
participant-level fields accompany them.

Input: a cohort `data.table` produced by `ukbsci-cohort` (or any
preprocessed UKB analysis dataset) that **already lives in RAP project
storage**. Output: aggregate coefficient tables (estimate / CI / p-value),
which are safe to export.

Forbidden: writing per-participant predicted values / residuals / influence
diagnostics to local disk. Keep those in `/mnt/project/...` if you need them
for QC.

---

## 1. When to load

- Fit one or several univariate / multivariable Cox, logistic, or linear
  models on a UKB cohort.
- Batch-run regression across **many exposures at once** (each gets its own
  model).
- Check the proportional-hazards assumption for Cox models.
- Run a lagged sensitivity analysis (drop events in the first N years).
- Fine-Gray competing-risks modeling.
- Linear-trend ("p-trend") testing across ordered exposure categories.

## 2. When NOT to load

- Building the cohort or survival dataset → `ukbsci-cohort`.
- Subgroup × interaction tests → `ukbsci-subgroup-sensitivity`.
- IPTW / propensity-adjusted regression → `ukbsci-propensity`.
- Multiple-imputation pooling → `ukbsci-imputation`.
- KM curves / cumulative incidence / time-dependent ROC → `ukbsci-survival`.
- Forest plots of the output → `ukbsci-plot`.

---

## 3. Prerequisites

```r
library(UKBAnalytica)
library(data.table)
library(survival)         # Cox

stopifnot(data.table::is.data.table(cohort))
# For Cox: outcome columns from build_survival_dataset()
stopifnot(all(c("outcome_status", "outcome_surv_time") %in% colnames(cohort)))
```

For competing-risks models the user needs a separate competing-event
indicator (e.g. non-cardiovascular death dates), or a single status column
encoded as 0 = censored, 1 = primary event, 2 = competing event.

---

## 4. Pipeline (4 phases)

### Phase 1 — Pick the right model family

| Question | Function | Outcome |
|----------|----------|---------|
| Continuous outcome | `runmulti_lm()` | numeric (e.g. SBP, BMI, biomarker level) |
| Binary outcome at fixed time | `runmulti_logit()` | 0/1 |
| Time-to-event with censoring | `runmulti_cox()` | `outcome_status` + `outcome_surv_time` |
| Time-to-event + competing risk | `runmulti_competing()` | status with `{0,1,2}` coding |
| Cox with PH diagnostics | `runmulti_cox_zph()` | same as Cox |
| Cox with early-event exclusion | `runmulti_cox_lag()` | same as Cox |
| Ordered exposure trend | `runmulti_trend()` | depends on `model_type` |

### Phase 2 — Choose exposures & covariates

The `*_multi_*` family runs **one model per `main_var`** but shares the same
covariate set. Build the covariate set from the analysis plan, e.g.

```r
covars <- c("age", "sex", "bmi", "smoking_status",
            "Hypertension_history", "Diabetes_history")
exposures <- c("sbp", "ldl_cholesterol", "hba1c", "crp")
```

### Phase 3 — Fit

```r
# Linear: continuous outcome, one model per exposure
lm_res <- runmulti_lm(
  data      = cohort,
  main_var  = exposures,
  covariates = covars,
  outcome   = "carotid_imt"
)

# Logistic: binary outcome at baseline
logit_res <- runmulti_logit(
  data      = cohort,
  main_var  = exposures,
  covariates = covars,
  outcome   = "Hypertension_history"
)

# Cox: time-to-event
cox_res <- runmulti_cox(
  data      = cohort,
  main_var  = exposures,
  covariates = covars,
  endpoint  = c("outcome_surv_time", "outcome_status")
)

# Save aggregate results
fwrite(cox_res, "/mnt/project/<area>/04-results/02-cox_main.csv")
```

### Phase 4 — Diagnostics & extensions

#### 4a. PH diagnostics (Schoenfeld residual test)

```r
zph <- runmulti_cox_zph(
  data       = cohort,
  main_var   = exposures,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status"),
  transform  = "km",
  alpha      = 0.05
)
zph[, .(variable, HR, lower95, upper95, pvalue,
        zph_pvalue, ph_violation)]
```

A `TRUE` in `ph_violation` means the exposure may have a time-varying effect;
consider stratifying the offending covariate or adding a `tt()` term in a
follow-up model.

#### 4b. Lagged sensitivity

Drops participants whose event occurred within `lag_years` of baseline —
guards against reverse causation.

```r
lag_res <- runmulti_cox_lag(
  data      = cohort,
  main_var  = exposures,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status"),
  lag_years  = c(0, 1, 2, 5),
  verbose    = TRUE
)
fwrite(lag_res, "/mnt/project/<area>/04-results/04-cox_lag.csv")
```

The output is *long*: one row per `(lag_years, variable, level)` triple.

#### 4c. Competing risks (Fine-Gray subdistribution hazard)

Two coding modes:

```r
# Single-column mode: status_col has {0 = censored, 1 = primary, 2 = competing}
fg1 <- runmulti_competing(
  data           = cohort,
  main_var       = exposures,
  covariates     = covars,
  time_col       = "outcome_surv_time",
  event_col      = "outcome_status_full",
  event_value    = 1,
  compete_value  = 2
)

# Dual-column mode: separate columns for primary and competing events
fg2 <- runmulti_competing(
  data        = cohort,
  main_var    = exposures,
  covariates  = covars,
  time_col    = "outcome_surv_time",
  event_col   = "primary_event",
  compete_col = "competing_event"
)
```

Output: subdistribution hazard ratios (`SHR`) instead of `HR`.

#### 4d. Trend across ordered exposure categories

```r
trend_res <- runmulti_trend(
  data            = cohort,
  main_var        = c("ldl_cholesterol_q", "sbp_q"),  # quartile factors
  covariates      = covars,
  model_type      = "cox",
  endpoint        = c("outcome_surv_time", "outcome_status"),
  ref_level       = "Q1",
  score_method    = "integer",                        # 1, 2, 3, 4
  include_level_estimates = TRUE
)
trend_res[, .(variable, level, estimate, lower95, upper95,
              pvalue, p_trend, trend_estimate)]
```

Output: per-level estimates **plus** an overall `p_trend` per exposure.

#### 4e. Bespoke PH diagnostics for a single fitted Cox object

```r
fit <- survival::coxph(
  survival::Surv(outcome_surv_time, outcome_status) ~ sbp + age + sex,
  data = cohort
)
diag <- ukb_cox_diagnostics(fit, transform = "km", alpha = 0.05)
diag$table
diag$global_pvalue
```

---

## 5. Common pitfalls

1. **Prevalent NA in `outcome_status`.** `build_survival_dataset()`
   intentionally codes prevalent cases as `NA` so Cox drops them. Do **not**
   impute these to 0.
2. **Factor reference level.** Both `runmulti_logit()` and `runmulti_cox()`
   honor the factor's first level as reference. Use `relevel()` (or set
   `ref_level` in `runmulti_trend`) before fitting if you want a different
   reference.
3. **Multicollinearity.** Batch regressions don't auto-check VIF. If
   `covars` contains correlated variables (e.g. `sbp` and `hypertension`),
   confirm with the user.
4. **Wrong endpoint order.** `endpoint = c("time", "status")` — time first.
   Reversing it makes the Cox model meaningless.
5. **Competing-risks status coding.** In single-column mode, the *value*
   `1` defaults to the primary event and `2` to the competing event. Confirm
   with the user; mis-coding inverts the result.
6. **`runmulti_cox_lag()` drops events with `surv_time ≤ lag_years`.** It
   does NOT shift the time origin; you are running the same regression on a
   trimmed cohort.
7. **`runmulti_trend()` requires factor exposures.** If `main_var` is
   continuous, the function errors out. Quartile-categorise first.
8. **`p_trend` interpretation.** `score_method = "integer"` is "score as
   linear" — assumes equal spacing between levels. Use `"custom"` with
   explicit scores when category midpoints are uneven.

---

## 6. Key functions (cheat-sheet)

| Function | Returns | Notable columns |
|----------|---------|-----------------|
| `runmulti_lm(data, main_var, covariates, outcome)` | data.frame | `variable, beta, lower95, upper95, pvalue` |
| `runmulti_logit(data, main_var, covariates, outcome)` | data.frame | `variable, OR, lower95, upper95, pvalue` |
| `runmulti_cox(data, main_var, covariates, endpoint)` | data.frame | `variable, HR, lower95, upper95, pvalue, n, n_event` |
| `runmulti_cox_zph(..., transform, alpha)` | data.frame | + `zph_pvalue, ph_violation, ph_global_violation` |
| `runmulti_cox_lag(..., lag_years)` | data.frame (long) | + `lag_years, n_removed` |
| `runmulti_competing(data, main_var, time_col, event_col, compete_col, event_value, compete_value)` | data.frame | `SHR, lower95, upper95, pvalue, n_compete` |
| `runmulti_trend(..., model_type, score_method, custom_scores)` | data.frame | per-level + `p_trend, trend_estimate` |
| `ukb_cox_diagnostics(model, transform, alpha)` | list | `table, global_pvalue, cox_zph` |

See [`references/functions.md`](references/functions.md) for full argument
tables.

---

## 7. Related skills

| Skill | When to switch |
|-------|----------------|
| `ukbsci-cohort` | Build / re-build the survival dataset before regressing. |
| `ukbsci-survival` | KM curves, time-dependent ROC, post-Cox visualization. |
| `ukbsci-propensity` | IPTW-weighted regression; PS matching analysis. |
| `ukbsci-subgroup-sensitivity` | Subgroup × interaction tests. |
| `ukbsci-mediation` | Mediation modeling. |
| `ukbsci-imputation` | Pool regression coefficients across imputations. |
| `ukbsci-plot` | Forest plot, volcano plot of regression output. |
| `ukbsci-workflow` | End-to-end orchestrator; calls this skill in Phase 6. |

---

## 8. References

- [`references/functions.md`](references/functions.md) — full signatures
- [`references/examples.md`](references/examples.md) — three end-to-end
  examples (Cox main, Cox with lag, Fine-Gray competing risks)
