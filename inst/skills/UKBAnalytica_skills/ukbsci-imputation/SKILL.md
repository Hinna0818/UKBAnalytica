---
name: ukbsci-imputation
description: >
  Multiple imputation by chained equations and Rubin's-rules pooling for UK
  Biobank Research Analysis Platform (RAP) cohorts using UKBAnalytica.
  Wraps mice (run_imputation, create_imputation_list) and mitools::MIcombine
  (pool_mi_models, pool_custom_estimates) to produce m imputed datasets, fit
  regression models on each (fit_mi_models), and pool estimates with FMI
  diagnostics (S3 mi_pooled_result + plot_mi_pooled, plot_mi_diagnostics).
  Supports merging additional datasets (e.g. proteomics) into each imputed
  copy via the additional_data argument. Use this skill when the user asks
  for multiple imputation, MI pooling, Rubin's rules, FMI, or wants to merge
  per-imputation analyses across covariate / proteomics datasets. Triggers:
  multiple imputation, MI, mice, Rubin's rules, FMI, pool imputations,
  imputed dataset, 多重插补, 插补合并, /ukbsci-imputation. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-imputation — Multiple imputation and pooling for UKB cohorts

## 0. RAP guardrails

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

The list of `m` imputed `data.tables` is **per-participant** and must remain
inside RAP project storage. Pooled coefficient tables and `plot_mi_*` figures
are aggregate and may be shared with the local agent. Recommended on-disk layout:

```
/mnt/project/<area>/03-cohort/
├── imp_list.rds        # list of m completed data.tables (RAP-only)
└── imp_mids.rds        # the mice mids object (optional)
```

---

## 1. When to load

- Cohort has non-trivial missingness on important covariates
  (`> ~5–10%` per variable).
- The user wants Cox / logistic / linear estimates **pooled** across `m`
  imputations using Rubin's rules.
- Multi-mediator or multi-exposure analyses to be repeated under MI.
- Merging a *separate* RAP-resident table (proteomics, metabolomics) into
  each imputed copy.

## 2. When NOT to load

- Missingness is < 5% on every covariate and complete-case is adequate →
  `ukbsci-subgroup-sensitivity` (complete-case helper).
- The user wants a quick exploratory model, no pooling → `ukbsci-regression`.

---

## 3. Prerequisites

```r
library(UKBAnalytica); library(mice); library(mitools); library(data.table)

# Optional model-specific suggestions
library(survival)        # Cox
library(MASS)            # neg-binomial
```

`run_imputation()` and `pool_mi_models()` error out with clear messages if
`mice` / `mitools` are not installed.

---

## 4. Pipeline (4 phases)

### Phase 1 — Impute

```r
imp <- run_imputation(
  data            = cohort,
  id_col          = "eid",
  vars            = c("bmi","ldl_cholesterol","hba1c","crp","sbp","dbp",
                       "smoking_status","alcohol_intake"),
  factor_vars     = c("smoking_status","alcohol_intake"),
  method          = "pmm",          # default predictive mean matching
  m               = 5,
  maxit           = 10,
  seed            = 1234,
  print           = TRUE,
  additional_data = list(proteomics = protein_dt),    # optional
  additional_join = "inner"
)

# imp$imp       → mice mids object (use mice::plot(imp$imp) for trace plots)
# imp$data_list → list of m completed data.tables (each = full row schema)
saveRDS(imp$data_list, "/mnt/project/<area>/03-cohort/imp_list.rds")
```

### Phase 2 — Fit models per imputation

```r
models <- fit_mi_models(
  datasets   = imp$data_list,                    # or imp_list directly
  formula    = survival::Surv(outcome_surv_time, outcome_status) ~
                 sbp + age + sex + bmi + smoking_status,
  model_type = "cox"
)
# models is a list of fitted coxph objects; failed fits are dropped with a warning
```

Alternative: pass an `imputationList` and have `pool_mi_models()` refit:

```r
imp_list <- create_imputation_list(imp$data_list)
class(imp_list)  # "imputationList"
```

### Phase 3 — Pool via Rubin's rules

```r
pooled <- pool_mi_models(
  models       = models,
  formula      = survival::Surv(outcome_surv_time, outcome_status) ~
                   sbp + age + sex + bmi + smoking_status,
  model_type   = "cox",
  df.complete  = Inf,                              # large-sample
  conf.level   = 0.95,
  exponentiate = NULL                              # auto: TRUE for cox
)
class(pooled)  # "mi_pooled_result"

# Quick views
summary(pooled, exponentiate = TRUE)
tidy(pooled, conf.int = TRUE, exponentiate = TRUE)
coef(pooled);  confint(pooled);  vcov(pooled)

fwrite(tidy(pooled, exponentiate = TRUE),
       "/mnt/project/<area>/04-results/02e-cox_pooled.csv")
```

The `pooled$pooled` slot contains: `term, estimate, std.error, statistic,
df, p.value, conf.low, conf.high, fmi`.

### Phase 4 — Diagnostics & visualization

```r
p_pool <- plot_mi_pooled(
  pooled,
  terms        = c("sbp","bmi","smoking_status"),
  exponentiate = TRUE,
  null_value   = 1,
  show_fmi     = TRUE,
  title        = "Pooled Cox HRs across 5 imputations"
)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig10-mi_pooled_forest.pdf",
                p_pool, width = 8, height = 5)

p_fmi <- plot_mi_diagnostics(pooled, type = "fmi")
p_vr  <- plot_mi_diagnostics(pooled, type = "variance_ratio")
p_df  <- plot_mi_diagnostics(pooled, type = "df")
```

FMI ≥ ~0.3 suggests the variable's pooled SE is heavily inflated by
between-imputation variance — consider increasing `m`.

### Phase 5 — Pool custom (non-model) estimates

For derived quantities like mediation `TNIE` across imputations, give
`pool_custom_estimates()` the list of per-imputation point estimates and
variance-covariance matrices:

```r
pooled_med <- pool_custom_estimates(
  estimates = lapply(imp$data_list, function(d) {
    r <- run_mediation(d, exposure = "sbp_above_median",
                       mediator = "ldl_cholesterol",
                       outcome  = "outcome_surv_time",
                       endpoint = c("outcome_surv_time","outcome_status"),
                       mediator_type = "continuous",
                       outcome_type  = "cox")
    coef(r)
  }),
  variances = lapply(imp$data_list, function(d) {
    r <- run_mediation(d, exposure = "sbp_above_median",
                       mediator = "ldl_cholesterol",
                       outcome  = "outcome_surv_time",
                       endpoint = c("outcome_surv_time","outcome_status"),
                       mediator_type = "continuous",
                       outcome_type  = "cox")
    vcov(r$regmedint_obj)
  }),
  df.complete = Inf,
  conf.level  = 0.95,
  labels      = c("CDE","PNDE","TNIE","TNDE","PNIE","TE","PM")
)
summary(pooled_med)
```

---

## 5. Common pitfalls

1. **`run_imputation()` re-coerces types.** Variables in `vars` are coerced
   to numeric (or factor if listed in `factor_vars`). Variables outside
   `vars` are left alone and merged back on `id_col`.
2. **Empty strings → NA.** `run_imputation()` converts `""` to `NA` before
   imputing — confirm this is the desired behavior, especially if you
   encoded "Prefer not to answer" as `""` upstream.
3. **`m = 5` is the default**, fine for `< 30%` missingness. For heavier
   missingness, use `m = 20` or more.
4. **Convergence.** Inspect `mice::plot(imp$imp)` — if trace plots have not
   stabilised, increase `maxit`.
5. **Cox + small subgroups.** Some imputations may fail to converge.
   `fit_mi_models()` drops failures and warns; if fewer than 2 succeed it
   errors out.
6. **`exponentiate = NULL` auto-decides.** It is `TRUE` for `cox` and
   `logistic`, `FALSE` for `lm`. Override when needed.
7. **Robust SEs are not directly supported in the pooled output.** If you
   need IPTW + MI, weight inside each imputation by hand, then pass the
   fitted weighted models to `pool_mi_models(models = ...)`.
8. **Cluster the merge.** When using `additional_data`, every imputed copy
   gets the same external table merged on `id_col`. If the external table
   itself has missingness, impute those columns inside `vars` instead.

---

## 6. Key functions

| Function | Purpose | Returns |
|----------|---------|---------|
| `run_imputation(data, id_col, vars, factor_vars, method, m, maxit, seed, print, additional_data, additional_join)` | mice MI + merge back | list(`imp`, `data_list`) |
| `create_imputation_list(datasets, validate = TRUE)` | wrap as `mitools::imputationList` | imputationList |
| `fit_mi_models(datasets, formula, model_type, family = NULL, ...)` | per-imputation fits | list of model objects |
| `pool_mi_models(models = NULL, datasets = NULL, formula = NULL, model_type, family = NULL, df.complete = Inf, conf.level = 0.95, exponentiate = NULL)` | Rubin's rules pool | S3 `mi_pooled_result` |
| `pool_custom_estimates(estimates, variances, df.complete = Inf, conf.level = 0.95, labels = NULL)` | pool arbitrary effects | S3 `mi_pooled_result` |
| `plot_mi_pooled(mi_result, terms, exponentiate, null_value, title, colors, show_fmi)` | forest with FMI | ggplot |
| `plot_mi_diagnostics(mi_result, type = c("fmi","variance_ratio","df"))` | diagnostic bar | ggplot |
| S3 on `mi_pooled_result` | `print, summary, coef, confint, vcov, tidy` | varies |

---

## 7. Related skills

| Skill | When |
|-------|------|
| `ukbsci-cohort` | Build cohort + outcomes first. |
| `ukbsci-preprocess` | Negative-code sanitation before imputation. |
| `ukbsci-regression` | Single-dataset model as comparison anchor. |
| `ukbsci-mediation` | Repeat across imputations + `pool_custom_estimates()`. |
| `ukbsci-propensity` | Weighted models per imputation, then pool. |
| `ukbsci-plot` | Polished forest of pooled estimates. |

---

## 8. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md)
