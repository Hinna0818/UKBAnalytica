# `ukbsci-subgroup-sensitivity` — function reference

## `run_subgroup_analysis()`

```r
run_subgroup_analysis(data, exposure, outcome = NULL,
                      subgroup_var, covariates = NULL,
                      model_type = c("cox","logistic","linear"),
                      endpoint  = NULL,
                      ref_level = NULL)
```

| Arg | Meaning |
|-----|---------|
| `exposure` | exposure column |
| `subgroup_var` | categorical modifier |
| `model_type` | one of `cox`, `logistic`, `linear` |
| `endpoint` | length-2 for Cox |
| `ref_level` | reference level for subgroup factor |

**Returns:** data.frame:

`subgroup_var, subgroup, n, n_event, estimate (HR/OR/β), lower95, upper95,
pvalue, p_interaction`.

`p_interaction` is the interaction term p-value from the full-cohort model,
repeated across rows of the same `subgroup_var`.

---

## `run_multi_subgroup()`

```r
run_multi_subgroup(data, exposure, outcome = NULL,
                   subgroup_vars,
                   covariates = NULL,
                   model_type = c("cox","logistic","linear"),
                   endpoint   = NULL)
```

Iterates `run_subgroup_analysis()` over `subgroup_vars`; returns a long
concatenated data.frame. No multiple-comparisons correction.

---

## `sensitivity_exclude_early_events()`

```r
sensitivity_exclude_early_events(data,
                                 endpoint = c("outcome_surv_time","outcome_status"),
                                 n_years,
                                 copy = TRUE,
                                 verbose = TRUE)
```

Drops rows with `surv_time ≤ n_years`. **Does not shift the time origin.**

**Returns:** filtered data + `attr(., "sensitivity_info")`:
`method, endpoint, n_years, n_input, n_removed, n_output`.

Errors when `endpoint` not length 2, time not numeric, or status not 0/1.

---

## `sensitivity_exclude_missing_covariates()`

```r
sensitivity_exclude_missing_covariates(data, covariates,
                                       copy = TRUE,
                                       stepwise = FALSE,
                                       verbose = TRUE)
```

Complete-case filter. `stepwise = TRUE` records per-variable attrition.

**Returns:** filtered data + `attr(., "sensitivity_info")`, optionally
`attr(., "complete_case_flow")`.

Prefer multiple-imputation (`ukbsci-imputation`) over complete-case for
primary analysis; this function is a sensitivity check.
