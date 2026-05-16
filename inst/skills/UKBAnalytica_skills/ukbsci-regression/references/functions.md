# `ukbsci-regression` — function reference

All signatures verbatim from `UKBAnalytica` ≥ 0.6.2.2 (`R/regression.R` +
`R/regression_extensions.R`).

---

## `runmulti_lm()`

```r
runmulti_lm(data, main_var, covariates = NULL, outcome, ...)
```

| Arg | Type | Meaning |
|-----|------|---------|
| `data` | data.frame / data.table | analysis dataset |
| `main_var` | char vector | exposure column names (one model per var) |
| `outcome` | char | continuous outcome column |
| `covariates` | char vector or `NULL` | adjustment variables |

**Returns:** `data.frame` with `variable`, `beta`, `lower95`, `upper95`,
`pvalue`.

**Errors when:** `outcome` / `main_var` / `covariates` missing from `data`,
or `main_var` empty.

---

## `runmulti_logit()`

```r
runmulti_logit(data, main_var, covariates = NULL, outcome, ...)
```

Same shape as `runmulti_lm` but for binary (0/1) outcomes.

**Returns:** `variable, OR, lower95, upper95, pvalue`.

**Errors when:** `outcome` is not binary 0/1.

---

## `runmulti_cox()`

```r
runmulti_cox(data, main_var, covariates = NULL,
             endpoint = c("time", "status"), ...)
```

| Arg | Meaning |
|-----|---------|
| `endpoint` | length-2 char: `c(time_col, status_col)` |

**Returns:** `variable, coef, se, z, HR, lower95, upper95, pvalue, n, n_event`.

Requires the `survival` package.

---

## `runmulti_cox_lag()`

```r
runmulti_cox_lag(data, main_var, covariates = NULL,
                 endpoint = c("time", "status"),
                 lag_years = c(0, 1, 2, 5),
                 verbose = TRUE, ...)
```

For each value in `lag_years`, drops participants with `surv_time ≤ lag` and
refits. Output is **long**:

`lag_years, variable, level, n_input, n_removed, n, n_event, HR, lower95,
upper95, pvalue`.

Internally calls `sensitivity_exclude_early_events()`.

---

## `runmulti_cox_zph()`

```r
runmulti_cox_zph(data, main_var, covariates = NULL,
                 endpoint = c("time", "status"),
                 transform = c("km", "rank", "identity"),
                 alpha = 0.05,
                 keep_models = FALSE, ...)
```

Schoenfeld-residual PH test (via `survival::cox.zph`). Output adds
`zph_pvalue, zph_global_pvalue, ph_violation, ph_global_violation`.

If `keep_models = TRUE`, the fitted `coxph` objects are attached as the
attribute `"cox_models"`.

---

## `runmulti_competing()`

```r
runmulti_competing(data, main_var, covariates = NULL,
                   time_col, event_col,
                   compete_col = NULL,
                   event_value = 1, compete_value = 2,
                   conf_level = 0.95, ...)
```

Fine-Gray subdistribution hazards.

| Mode | Spec |
|------|------|
| Single-column | `compete_col = NULL`; `event_col` is `{0, event_value, compete_value}` |
| Dual-column | `event_col` and `compete_col` are separate 0/1 indicators |

**Returns:** `variable, level, n, n_event, n_compete, SHR, lower95, upper95,
pvalue`.

---

## `runmulti_trend()`

```r
runmulti_trend(data, main_var, outcome = NULL, covariates = NULL,
               model_type = c("cox", "logistic", "linear"),
               endpoint = NULL,
               ref_level = NULL,
               score_method = c("integer", "median", "custom"),
               custom_scores = NULL,
               include_level_estimates = TRUE, ...)
```

Tests dose-response trend across an **ordered factor** exposure.

| Arg | Meaning |
|-----|---------|
| `model_type` | which underlying model |
| `score_method` | `"integer"` (1,2,3…) or `"custom"` |
| `custom_scores` | named list: `list(exposure_name = c("Q1"=0, "Q2"=1, ...))` |

**Returns:** per-level rows + `p_trend, trend_estimate, trend_lower95,
trend_upper95, score_method`. `"median"` is reserved for future release.

---

## `ukb_cox_diagnostics()`

```r
ukb_cox_diagnostics(model,
                    transform = c("km", "rank", "identity"),
                    terms = TRUE, global = TRUE,
                    alpha = 0.05,
                    return_object = TRUE)
```

Run PH diagnostics on a single fitted `coxph` model.

**Returns:** `list(table, global_pvalue, cox_zph)` where `cox_zph` is the raw
`survival::cox.zph` object (only when `return_object = TRUE`).

---

## Dependencies

| Need | Suggests-package |
|------|------------------|
| Cox / Fine-Gray | `survival` |
| Robust SE (downstream of `run_weighted_analysis`) | `sandwich`, `lmtest` |

No participant-level data leaves the project — outputs are aggregate
coefficient tables.
