# `ukbsci-propensity` — function reference

All signatures reflect `UKBAnalytica` 1.0.0.

---

## `estimate_propensity_score()`

```r
estimate_propensity_score(data, treatment, covariates,
                          method = c("logistic", "gbm"),
                          formula = NULL)
```

| Arg | Type | Meaning |
|-----|------|---------|
| `treatment` | char | binary (0/1) treatment column |
| `covariates` | char vector | adjustment variables |
| `method` | char | `"logistic"` (default) or `"gbm"` |
| `formula` | formula | optional override; auto-built if `NULL` |

**Returns:** `data.table` with new column `ps` ∈ `[0.001, 0.999]`.

Requires `gbm` for `method = "gbm"`.

---

## `match_propensity()`

```r
match_propensity(data, ps_col = "ps", treatment,
                 ratio = 1, caliper = 0.2,
                 method = c("nearest", "optimal"),
                 replace = FALSE,
                 exact_match = NULL)
```

| Arg | Meaning |
|-----|---------|
| `ratio` | matching ratio (1 = 1:1) |
| `caliper` | caliper width in SD of PS |
| `method` | nearest-neighbor or optimal |
| `replace` | sampling with/without replacement |
| `exact_match` | char vector of variables to match exactly |

Wraps `MatchIt::matchit()`. **Returns:** matched `data.table` with
`match_id`, `match_distance`.

---

## `calculate_weights()`

```r
calculate_weights(data, ps_col = "ps", treatment,
                  weight_type = c("ATE", "ATT", "ATC"),
                  stabilized = TRUE,
                  truncate = c(0.01, 0.99))
```

IPTW weights:

| `weight_type` | Estimand |
|---------------|----------|
| `"ATE"` | Average treatment effect over the population |
| `"ATT"` | Effect on the treated |
| `"ATC"` | Effect on the controls |

Stabilized weights multiply by marginal P(T=1). Truncation clips at the
specified quantiles.

**Returns:** `data.table` with new column `weight`.

---

## `assess_balance()`

```r
assess_balance(data, treatment, covariates,
               method = c("unmatched", "matched", "weighted"),
               weight_col = NULL,
               threshold = 0.1)
```

**Returns:** `data.frame` with `variable, mean_treated, mean_control, smd,
variance_ratio, balanced`.

`balanced` is `TRUE` when `|smd| ≤ threshold`. Requires `weight_col` for
`method = "weighted"`.

---

## `plot_balance()`

```r
plot_balance(balance_before, balance_after,
             threshold = 0.1,
             title = "Covariate Balance",
             xlab  = "Standardized Mean Difference")
```

Love plot: pre vs post points joined by horizontal lines; vertical reference
lines at ±`threshold`.

**Returns:** `ggplot2` object.

---

## `plot_ps_distribution()`

```r
plot_ps_distribution(data, ps_col = "ps", treatment,
                     type = c("histogram", "density", "mirror"),
                     matched = FALSE,
                     match_col = NULL)
```

| `type` | Output |
|--------|--------|
| `"histogram"` | overlapping histograms by group |
| `"density"` | kernel density curves |
| `"mirror"` | back-to-back histograms (clinical-paper style) |

**Returns:** `ggplot2` object.

---

## `run_weighted_analysis()`

```r
run_weighted_analysis(data, exposure, outcome = NULL,
                      covariates = NULL,
                      weight_col = "weight",
                      model_type = c("cox", "logistic", "linear"),
                      endpoint   = NULL,
                      robust_se  = TRUE)
```

IPTW (or matched-weighted) regression with `sandwich` + `lmtest` robust SE
(HC0 for linear/logit; coxph `robust = TRUE` for Cox).

| `model_type` | Required |
|--------------|----------|
| `"cox"` | `endpoint = c(time_col, status_col)` |
| `"logistic"` / `"linear"` | `outcome` |

**Returns:** `data.frame` with `variable`, effect (HR/OR/β), `lower95`,
`upper95`, `pvalue`.

---

## Suggests-package dependencies

| Need | Package |
|------|---------|
| Logistic PS | base R |
| GBM PS | `gbm` |
| Matching | `MatchIt` |
| Balance | (built-in) — `cobalt` optional for cross-checking |
| Robust SE | `sandwich`, `lmtest` |
| Weighted Cox | `survival` |
