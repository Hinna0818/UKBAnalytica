# `ukbsci-subgroup-sensitivity` — function reference

## `run_subgroup_analysis()`

```r
run_subgroup_analysis(data, exposure, outcome = NULL,
                      subgroup_var, covariates = NULL,
                      model_type = c("cox","logistic","linear","glm","negbin"),
                      family = "poisson",
                      endpoint  = NULL,
                      ref_level = NULL)
```

| Arg | Meaning |
|-----|---------|
| `exposure` | exposure column |
| `subgroup_var` | categorical modifier |
| `model_type` | one of `cox`, `logistic`, `linear`, `glm`, `negbin` |
| `family` | GLM family for `model_type = "glm"` |
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
                   model_type = c("cox","logistic","linear","glm","negbin"),
                   family = "poisson",
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

---

## `ukb_participant_flow()` and `plot_participant_flow()`

```r
ukb_participant_flow(data, steps, id_col = NULL,
                     outcome_col = NULL, event_value = 1,
                     start_label = "Initial population")

plot_participant_flow(flow, show_removed = TRUE,
                      show_events = TRUE, fill = "#2F6FA3")
```

`steps` is a named list of one-sided formulas, functions, logical vectors, or
character vectors of complete-case variables. Returns an aggregate flow table
with retained/removed counts and event counts. The plotter returns a ggplot.

---

## `ukb_sensitivity_suite()`

```r
ukb_sensitivity_suite(data,
                      exposure,
                      covariates = NULL,
                      endpoint = c("outcome_surv_time", "outcome_status"),
                      early_event_years = c(2, 4, 6),
                      complete_case_covariates = NULL,
                      additional_covariate_sets = NULL,
                      conf_level = 0.95,
                      verbose = TRUE)
```

Fits a primary Cox model and common sensitivity models with the same endpoint
and exposure. Returns a `ukb_sensitivity_suite` list with `summary`, `models`,
`flows`, and `settings`.
