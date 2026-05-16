# `ukbsci-imputation` — function reference

## `run_imputation()`

```r
run_imputation(data, id_col = "eid", vars,
               factor_vars     = NULL,
               method          = "pmm",
               m               = 5,
               maxit           = 10,
               seed            = 1234,
               print           = TRUE,
               additional_data = NULL,
               additional_join = c("inner", "left"))
```

| Arg | Type | Meaning |
|-----|------|---------|
| `vars` | char vector | columns to impute |
| `factor_vars` | subset of `vars` | columns to coerce to factor |
| `method` | char or named vector | mice method (default `pmm`); name per-var to override |
| `m` | int | number of imputed datasets |
| `additional_data` | named list of data.frames | external tables merged into each imputed copy by `id_col` |

**Returns:** list:

- `imp` — `mice::mids` (use `mice::plot(imp$imp)` for trace diagnostics)
- `data_list` — list of `m` completed `data.tables` with all original columns
  + merged additional data

Empty strings are converted to `NA` before imputing.

---

## `create_imputation_list()`

```r
create_imputation_list(datasets, validate = TRUE)
```

Wraps a list of completed `data.frame`s into an `mitools::imputationList`.
Validates column-name and row-count consistency.

---

## `fit_mi_models()`

```r
fit_mi_models(datasets, formula,
              model_type = c("lm", "logistic", "poisson", "cox", "negbin"),
              family     = NULL,
              ...)
```

Per-imputation model fitting. Survival models require `survival`; negbin
requires `MASS`. Failed fits dropped with a warning; errors when fewer than
2 succeed.

**Returns:** list of fitted model objects.

---

## `pool_mi_models()`

```r
pool_mi_models(models       = NULL,
               datasets     = NULL,
               formula      = NULL,
               model_type   = c("lm","logistic","poisson","cox","negbin"),
               family       = NULL,
               df.complete  = Inf,
               conf.level   = 0.95,
               exponentiate = NULL)
```

Two calling modes: pass pre-fitted `models` (preferred when you need
unusual weighting), or pass `datasets` + `formula` and let the function
refit. Pooling uses `mitools::MIcombine` (Rubin's rules + FMI).

**Returns:** S3 `mi_pooled_result`:

| Slot | Meaning |
|------|---------|
| `pooled` | data.frame: `term, estimate, std.error, statistic, df, p.value, conf.low, conf.high, fmi` |
| `mi_result` | raw `MIresult` |
| `n_imputations` | int |
| `model_type` | char |
| `formula` | model formula |
| `exponentiated` | logical |
| `conf.level` | numeric |
| `call` | function call |

S3 methods on `mi_pooled_result`: `print`, `summary`, `coef`, `confint`,
`vcov`, `tidy` (broom-style with `conf.int`, `exponentiate`).

---

## `pool_custom_estimates()`

```r
pool_custom_estimates(estimates, variances,
                      df.complete = Inf,
                      conf.level  = 0.95,
                      labels      = NULL)
```

Pool an arbitrary vector of point estimates (per imputation) with their
variance-covariance matrices. Useful for mediation `TNIE`, calibration
intercepts, ML-derived risk scores, etc.

**Returns:** S3 `mi_pooled_result`.

---

## `plot_mi_pooled()` / `plot_mi_diagnostics()`

```r
plot_mi_pooled(mi_result, terms = NULL,
               exponentiate = NULL,
               null_value   = NULL,
               title        = "Pooled Estimates (Multiple Imputation)",
               colors       = NULL,
               show_fmi     = TRUE)

plot_mi_diagnostics(mi_result, type = c("fmi","variance_ratio","df"),
                    title = NULL)
```

- `plot_mi_pooled`: forest plot of pooled estimates; point size = FMI.
- `plot_mi_diagnostics`: per-term bar of FMI / between-within variance
  ratio / effective df, intercept excluded.

Both return `ggplot2` objects.

---

## Rubin's-rules summary (what the pool does)

- Pooled point estimate = mean of `m` per-imputation estimates.
- Pooled variance = mean within-imputation variance + (1 + 1/m) × between-imputation variance.
- FMI ≈ (1 + 1/m) × between / total variance.
- DF correction shrinks ∞ → small-sample value when between-variance dominates.
- CIs use `t(df)`; p-values via the same approximation.

---

## Dependency summary

| Need | Suggests |
|------|----------|
| Imputation | `mice` |
| Pooling | `mitools` |
| Cox | `survival` |
| Negative-binomial | `MASS` |
| Robust SE (per-imputation) | `sandwich`, `lmtest` |
