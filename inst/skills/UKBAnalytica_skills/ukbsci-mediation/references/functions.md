# `ukbsci-mediation` — function reference

Requires `regmedint`. Decomposition follows Valeri & VanderWeele (2013):
4-way decomposition of TE = CDE + PNIE + (interaction reference) +
(interaction mediated).

---

## `run_mediation()`

```r
run_mediation(data, exposure, mediator, outcome,
              covariates       = NULL,
              exposure_levels  = c(0, 1),
              mediator_value   = 0,
              covariate_values = NULL,
              mediator_type    = c("continuous", "binary"),
              outcome_type     = c("linear", "logistic", "cox"),
              endpoint         = NULL,
              interaction      = TRUE,
              boot             = FALSE,
              boot_n           = 1000,
              conf_level       = 0.95)
```

| Arg | Meaning |
|-----|---------|
| `exposure_levels` | comparison contrast (treated vs. reference) |
| `mediator_value` | value at which CDE is evaluated |
| `covariate_values` | reference covariate vector; `NULL` → covariate means |
| `endpoint` | for Cox, `c(time_col, status_col)` |
| `interaction` | include E × M term in outcome model (default `TRUE`) |
| `boot` | bootstrap CIs (vs. delta method) |

**Returns:** S3 `mediation_result` list:

| Slot | Contents |
|------|----------|
| `effects` | data.frame: CDE, PNDE, TNIE, TNDE, PNIE, TE, PM with SE / CI / p |
| `mediator_model` | fitted mediator model |
| `outcome_model` | fitted outcome model |
| `regmedint_obj` | raw `regmedint` object |
| `params` | analysis parameters |

S3 methods: `print`, `summary` (with `exponentiate` arg), `coef`, `confint`.

---

## `run_multi_mediator()`

```r
run_multi_mediator(data, exposure, mediators, outcome,
                   covariates    = NULL,
                   mediator_type = "continuous",
                   outcome_type  = "linear",
                   endpoint      = NULL,
                   ...)
```

Iterates `run_mediation()` across `mediators`. **Returns:** `data.frame` with
one row per mediator: `mediator, tnie, tnie_se, tnie_lower, tnie_upper,
tnie_p, pnde, pnde_se, pnde_p, te, te_se, te_p, pm, pm_se, pm_p`.

Failures for a single mediator yield an `NA` row + warning.

No multiple-comparisons correction — apply externally (e.g. BH).

---

## `plot_mediation()`

```r
plot_mediation(mediation_result,
               type = c("effects", "path", "decomposition"),
               show_ci       = TRUE,
               show_pvalue   = TRUE,
               exponentiate  = FALSE,
               title         = NULL,
               colors        = NULL)
```

| `type` | Output |
|--------|--------|
| `"effects"` | bar chart of CDE / PNDE / TNIE / TE |
| `"path"` | path diagram (a × b vs. c′) |
| `"decomposition"` | stacked bar of TE = NDE + NIE |

**Returns:** `ggplot2` object.

---

## `plot_mediation_forest()`

```r
plot_mediation_forest(multi_mediation_result,
                      effect_type  = c("tnie", "pnde", "te", "pm"),
                      exponentiate = FALSE,
                      null_value   = 0,
                      title        = "Mediation Analysis: Multiple Mediators")
```

Forest plot from `run_multi_mediator()` output, ordered by effect size.
**Returns:** `ggplot2` object.
