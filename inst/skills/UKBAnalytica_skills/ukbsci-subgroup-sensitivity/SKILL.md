---
name: ukbsci-subgroup-sensitivity
description: >
  Subgroup interaction tests and sensitivity-analysis filters for UK Biobank
  Research Analysis Platform (RAP) cohorts built with UKBAnalytica. Covers
  run_subgroup_analysis and run_multi_subgroup (stratum-specific Cox / logistic /
  linear / GLM / negative-binomial effects + interaction p-value),
  sensitivity_exclude_early_events (drop events within N years of baseline
  to guard against reverse causation), and sensitivity_exclude_missing_covariates
  (complete-case analysis with optional flow tracking), ukb_participant_flow /
  plot_participant_flow for cohort attrition, and ukb_sensitivity_suite for
  common Cox sensitivity runs. Use this skill when
  the user asks for subgroup analysis, interaction tests, effect modification,
  complete-case sensitivity, or early-event exclusion sensitivity. Subgroup
  analysis now supports all model families: cox, logistic, linear, glm
  (any GLM family including Poisson, quasi-Poisson, Gamma), and negbin
  (negative-binomial). Triggers: subgroup analysis, interaction p-value,
  effect modification, sensitivity analysis, complete-case analysis, lag analysis,
  亚组分析, 敏感性分析, /ukbsci-subgroup-sensitivity. Hard rule: local agents must
  not read or inspect real UKB RAP participant-level data; generate scripts for
  RAP execution and interpret aggregate outputs only.
---

# ukbsci-subgroup-sensitivity — Subgroup & sensitivity analyses

## 0. RAP guardrails

Input: cohort `data.table`. Output: aggregate subgroup-effect tables and
sensitivity flow tables that may be shared with the local agent.

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

---

## 1. When to load

- Run effect estimates **within strata** of a categorical modifier.
- Report a formal **interaction p-value** for the modifier.
- Sweep over multiple potential modifiers in one go.
- Drop events in the first N years (lag analysis for reverse-causation).
- Complete-case re-run, optionally with a participant-flow trace.
- Record sequential inclusion/exclusion counts for a flow diagram.
- Run a compact Cox sensitivity suite with early-event and covariate-adjustment variants.

## 2. When NOT to load

- Main-effect models without subgroups → `ukbsci-regression`.
- Mediation analysis → `ukbsci-mediation`.
- Propensity-score weighting → `ukbsci-propensity`.

---

## 3. Prerequisites

```r
library(UKBAnalytica); library(survival); library(data.table)
stopifnot(all(c("outcome_status", "outcome_surv_time") %in% colnames(cohort)))
```

---

## 4. Pipeline

### Phase 0 — Participant flow table

Use `ukb_participant_flow()` when the user needs a transparent attrition table
before modeling. The output is aggregate and can be shared with the local
agent.

```r
flow <- ukb_participant_flow(
  data = cohort,
  steps = list(
    "Complete primary covariates" = c("age", "sex", "bmi", "smoking_status"),
    "Positive follow-up time"     = ~ outcome_surv_time > 0
  ),
  id_col      = "eid",
  outcome_col = "outcome_status"
)

plot_participant_flow(flow)
```

### Phase 1 — Single-modifier subgroup analysis

`run_subgroup_analysis()` estimates stratum-specific effects **and** an interaction
p-value for any of five model families: `"cox"`, `"logistic"`, `"linear"`, `"glm"`,
`"negbin"`.

```r
# Cox subgroup (default use-case)
sub_cox <- run_subgroup_analysis(
  data         = cohort,
  exposure     = "sbp",
  subgroup_var = "sex",
  covariates   = c("age", "bmi", "smoking_status"),
  model_type   = "cox",
  endpoint     = c("outcome_surv_time", "outcome_status"),
  ref_level    = "Male"
)
# data.frame: one row per subgroup level
#   subgroup_var, subgroup, n, n_event, estimate (HR / OR / β / IRR),
#   lower95, upper95, pvalue, p_interaction

# GLM subgroup — Poisson (incidence rate ratio)
sub_glm <- run_subgroup_analysis(
  data         = cohort,
  exposure     = "sbp",
  outcome      = "hospitalisation_count",
  subgroup_var = "sex",
  covariates   = c("age", "bmi", "smoking_status"),
  model_type   = "glm",
  family       = "poisson"
)
# family accepts character, function, or family object
# quasi-Poisson → Wald CI; log/logit link → estimate column is exponentiated

# Negative-binomial subgroup (IRR = exp(beta))
sub_nb <- run_subgroup_analysis(
  data         = cohort,
  exposure     = "sbp",
  outcome      = "hospitalisation_count",
  subgroup_var = "sex",
  covariates   = c("age", "bmi", "smoking_status"),
  model_type   = "negbin"
)
```

**Interaction p-value method:**
- Binary subgroup or continuous exposure × binary factor → Wald p on single interaction term.
- Multi-level factor → likelihood-ratio test (LRT) via `anova()`; quasi-likelihood falls back
  to an F-test.

### Phase 2 — Multiple modifiers

```r
multi_sub <- run_multi_subgroup(
  data           = cohort,
  exposure       = "sbp",
  subgroup_vars  = c("sex", "age_band", "smoking_status", "Diabetes_history"),
  covariates     = c("age", "bmi"),
  model_type     = "cox",           # or "glm", "negbin", "logistic", "linear"
  endpoint       = c("outcome_surv_time", "outcome_status")
)
fwrite(multi_sub, "/mnt/project/<area>/04-results/03-cox_subgroup.csv")
```

The output is long; `p_interaction` is repeated for each level of a given subgroup variable.
`model_type` and `family` are forwarded to `run_subgroup_analysis()` — all five families apply.

### Phase 3 — Sensitivity analysis: drop early events

```r
trim <- sensitivity_exclude_early_events(
  data     = cohort,
  endpoint = c("outcome_surv_time", "outcome_status"),
  n_years  = 2,
  copy     = TRUE,
  verbose  = TRUE
)
attr(trim, "sensitivity_info")
# list: method, endpoint, n_years, n_input, n_removed, n_output

# Then re-fit:
cox_lag2 <- runmulti_cox(
  data = trim,
  main_var = "sbp",
  covariates = c("age", "sex", "bmi", "smoking_status"),
  endpoint = c("outcome_surv_time", "outcome_status")
)
```

Use `runmulti_cox_lag()` instead when you want **the same regression
fit at multiple lag values in one call** (see `ukbsci-regression`).

### Phase 4 — Sensitivity analysis: complete-case

```r
cc <- sensitivity_exclude_missing_covariates(
  data       = cohort,
  covariates = c("age", "sex", "bmi", "smoking_status",
                 "ldl_cholesterol", "hba1c"),
  copy       = TRUE,
  stepwise   = TRUE,                         # records per-variable drop
  verbose    = TRUE
)
attr(cc, "sensitivity_info")
attr(cc, "complete_case_flow")  # only when stepwise = TRUE

cox_cc <- runmulti_cox(
  data = cc,
  main_var = "sbp",
  covariates = c("age", "sex", "bmi", "smoking_status",
                 "ldl_cholesterol", "hba1c"),
  endpoint = c("outcome_surv_time", "outcome_status")
)
fwrite(cox_cc, "/mnt/project/<area>/04-results/04-cox_sensitivity.csv")
```

When the user has multiply-imputed data, prefer pooled MI analysis
(`ukbsci-imputation`) over complete-case. Complete-case is a sensitivity
check, not the main analysis.

### Phase 5 — Compact Cox sensitivity suite

Use `ukb_sensitivity_suite()` when the user wants the standard Cox sensitivity
bundle with one exposure and shared covariates.

```r
sens <- ukb_sensitivity_suite(
  data = cohort,
  exposure = "exposure",
  covariates = c("age", "sex", "bmi", "smoking_status"),
  endpoint = c("outcome_surv_time", "outcome_status"),
  early_event_years = c(2, 4, 6),
  complete_case_covariates = c("age", "sex", "bmi", "smoking_status"),
  additional_covariate_sets = list(extra = c("education", "TDI")),
  verbose = FALSE
)

sens$summary
```

---

## 5. Common pitfalls

1. **Subgroup variable with NA.** Rows with `NA` on `subgroup_var` are
   dropped silently by `run_subgroup_analysis()`. The reported `n` reflects
   the post-drop sample.
2. **No within-stratum variation.** A subgroup level may have zero events
   for a binary exposure — the function warns and returns `NA` for that
   row. Filter / merge upstream.
3. **`p_interaction` interpretation.** It is the p-value for the
   exposure × subgroup interaction term in the full-cohort model — NOT a
   pairwise stratum comparison.
4. **Multiple-comparisons correction.** `run_multi_subgroup()` does NOT
   adjust. With 4+ modifiers, apply BH externally.
5. **`sensitivity_exclude_early_events()` does not shift origin.** It
   simply drops events with `surv_time ≤ n_years`; the time variable is
   unchanged. Mention this in the manuscript footnote.
6. **`copy = TRUE` is the safe default** for `data.table` inputs —
   it preserves the original. Setting `FALSE` mutates `cohort` in place.
7. **`stepwise = TRUE` is informative but slow** for wide cohorts. Use it
   only when you need a participant-flow figure for the manuscript.

---

## 6. Key functions

| Function | Returns |
|----------|---------|
| `ukb_participant_flow(data, steps, id_col, outcome_col)` | aggregate attrition table |
| `plot_participant_flow(flow)` | ggplot |
| `run_subgroup_analysis(data, exposure, outcome = NULL, subgroup_var, covariates = NULL, model_type = c("cox","logistic","linear","glm","negbin"), family = "poisson", endpoint = NULL, ref_level = NULL)` | data.frame: per-level estimates + `p_interaction` |
| `run_multi_subgroup(..., subgroup_vars, model_type, family)` | concatenated long table |
| `sensitivity_exclude_early_events(data, endpoint, n_years, copy = TRUE, verbose = TRUE)` | filtered data + `sensitivity_info` attr |
| `sensitivity_exclude_missing_covariates(data, covariates, copy = TRUE, stepwise = FALSE, verbose = TRUE)` | filtered data + `sensitivity_info` (+ `complete_case_flow` if stepwise) |
| `ukb_sensitivity_suite(data, exposure, covariates, endpoint, early_event_years, complete_case_covariates, additional_covariate_sets)` | list with `summary`, `models`, `flows`, `settings` |

**`model_type` effect estimate column mapping:**

| `model_type` | `estimate` column | CI method |
|---|---|---|
| `"cox"` | Hazard ratio (HR) | `summary(coxph)$conf.int` |
| `"logistic"` | Odds ratio (OR) | Profile-likelihood |
| `"linear"` | Beta (mean difference) | `confint(lm)` |
| `"glm"` log/logit link | Ratio `exp(beta)` | Profile-likelihood (Wald for quasi) |
| `"glm"` identity link | Beta | Profile-likelihood |
| `"negbin"` | IRR `exp(beta)` | Profile-likelihood |

---

## 7. Related skills

| Skill | When |
|-------|------|
| `ukbsci-cohort` | Build cohort with `outcome_status` / `outcome_surv_time`. |
| `ukbsci-regression` | Main-effect models for comparison. |
| `ukbsci-imputation` | Pooled MI analysis (preferred over complete-case). |
| `ukbsci-plot` | Subgroup forest plot of the output. |

---

## 8. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md)
