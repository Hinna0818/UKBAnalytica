---
name: ukbsci-ml
description: >
  Machine learning workflows (binary / multiclass classification and
  right-censored survival) for UK Biobank Research Analysis Platform (RAP)
  cohorts using UKBAnalytica. End-to-end pipelines (ukb_ml_workflow,
  ukb_ml_flow, ukb_ml_survival_workflow) cover frozen train / validation /
  test split, feature selection (Boruta, filter, glmnet),
  hyperparameter tuning (grid / random / Bayesian), final refit,
  test-set evaluation, model comparison, and model-agnostic interpretation
  via SHAP. Supports RF / XGBoost / glmnet / SVM / NNet / RPart / Naive
  Bayes / Logistic for classification, and RSF / GBM / coxnet for survival.
  Output evaluation: ROC, PR, calibration, DCA, KS, gain & lift, confusion
  matrix, threshold optimization, SHAP summary / dependence / force, plus
  plot_ml_* and plot_shap_* visualizations. Use this skill when the user
  asks for a UKB-based predictive model, ML pipeline, SHAP interpretation,
  C-index for a survival model, or comparison between several ML algorithms.
  Triggers: UKB ML, machine learning, XGBoost, random forest, SHAP, C-index,
  survival ML, calibration, decision curve, AUC, /ukbsci-ml. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-ml — Machine learning workflows for UKB cohorts

## 0. RAP guardrails

**Core rule:** do not export UK Biobank participant-level raw data, direct
identifiers (`eid`, exact dates, raw RAP fields), or row-level source tables
that can be linked back to participants.

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

| Output | Sharing rule |
|--------|-------------|
| UK Biobank RAP raw participant-level data | Do not export. |
| Direct identifiers (`eid`, raw dates, RAP field values) | Do not export. |
| Train / validation / test datasets with raw fields or `eid` | Do not export. |
| Trained model `.rds` | Keep in RAP project workspace unless sharing is explicitly permitted by project governance. |
| Per-row predictions containing `eid` or re-identifiable row index | Do not export. |
| Per-row predictions without any identifier | Keep in RAP; do not share with the local agent. Export only aggregate metrics or rendered figures. |
| Aggregate metrics table (AUC, C-index, calibration intercept, RMSE…) | Shareable with the local agent. |
| ROC / PR / calibration / DCA / KS / gain plots | Shareable as rendered figures. |
| SHAP global importance table (mean \|SHAP\| per feature) | Shareable — aggregate only. |
| SHAP beeswarm / dependence / summary plots | Shareable as rendered figures only, provided no row-level source data, row indices, or identifiers accompany the file. |
| SHAP force plot for a real participant row | Do not share with the local agent. Use a synthetic or representative prototype instead. |

---

## 1. When to load

- User wants a predictive model for a binary outcome (e.g. 10-year MACE
  risk), a multiclass outcome (diabetes subtype), or a survival endpoint.
- User wants to compare several algorithms on the same train / test split.
- User wants SHAP-based interpretation.
- User wants Decision Curve Analysis or threshold optimization.

## 2. When NOT to load

- Effect estimation / inferential modeling → `ukbsci-regression`.
- Mediation, propensity scoring → respective skills.
- Just a Cox model without ML → `ukbsci-regression`.

---

## 3. Prerequisites

```r
library(UKBAnalytica); library(data.table); library(survival)

# Suggests (install on RAP only what your chosen model needs):
# ranger xgboost glmnet e1071 nnet rpart randomForestSRC gbm MASS
# Boruta rBayesianOptimization fastshap pROC
```

Before using high-level workflow helpers, verify that the latest
UKBAnalytica code is installed or loaded in the RAP session. If
`ukb_ml_flow()` or `ukb_ml_compare_flows()` is not found, reinstall the
current package version using the project-approved method or load the current
source tree with `pkgload::load_all()`.

After ML interface changes, run the shared skill smoke test:

```r
Rscript inst/skills/UKBAnalytica_skills/evals/smoke_test_core.R
```

---

## 4. Recommended entry points

**Use the right function for your task:**

| Task | Recommended function |
|------|----------------------|
| Single model, most analyses | `ukb_ml_flow()` |
| Compare several algorithms or feature sets | `ukb_ml_compare_flows()` |
| Fully staged workflow with explicit feature selection + per-stage artefacts | `ukb_ml_workflow()` |
| Manual step-by-step control | lower-level functions (`ukb_ml_split_data`, `ukb_ml_tune`, …) |

### Return-object cheat-sheet

| Function | Key slots in the returned object |
|----------|----------------------------------|
| `ukb_ml_workflow()` | `$split`, `$feature_result`, `$tune_result`, `$threshold_result`, `$final_model`, `$final_test_metrics`, `$final_test_predictions` |
| `ukb_ml_flow()` | `$model`, `$metrics`, `$predictions`, `$roc`, `$shap` (if `compute_shap = TRUE`) |
| `ukb_ml_compare_flows()` | `$flows`, `$metrics`, `$comparison`, `$predictions`, `$roc`, `$thresholds` |

Never mix components across workflow objects (e.g. `workflow_obj$roc` does not
exist; use `flow_obj$roc`).

---

### 4a. `ukb_ml_flow()` — single model (recommended for most analyses)

```r
flow <- ukb_ml_flow(
  formula      = outcome_status ~ age + sex + bmi + sbp + ldl_cholesterol +
                                   hba1c + crp + smoking_status,
  data         = cohort,
  model        = "xgboost",
  outcome_type = "binary",
  split_params = list(split = "train_valid_test",
                       train_ratio = 0.7,
                       validation_ratio = 0.1,
                       test_ratio = 0.2,
                       stratify_by = "outcome"),
  tune    = TRUE,
  tune_params = list(search = "random", n_iter = 30,
                      resampling = "validation",
                      metric = "auc", maximize = TRUE),
  threshold_method = "youden",
  compute_shap = TRUE,
  shap_params  = list(nsim = 100, sample_n = 2000),
  seed    = 1234,
  verbose = TRUE
)

flow$metrics                   # aggregate — shareable with the local agent
# flow$predictions              # keep on RAP (may contain row indices)

ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig13-roc.pdf",
                plot_ml_roc(flow))
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig14-calib.pdf",
                plot_ml_calibration(flow))
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig15-dca.pdf",
                plot_ml_dca(flow))
```

### 4b. `ukb_ml_workflow()` — staged workflow with explicit feature selection

Use when you need per-stage artefacts (Boruta importance table, tuning grid,
held-out validation metrics).

```r
flow <- ukb_ml_workflow(
  formula        = outcome_status ~ .,
  data           = cohort,
  model          = "xgboost",
  outcome_type   = "binary",
  split_params   = list(split = "train_valid_test",
                         train_ratio = 0.7, validation_ratio = 0.1,
                         test_ratio = 0.2, stratify_by = "outcome"),
  feature_select = "boruta",
  tune           = TRUE,
  tune_params    = list(search = "random", n_iter = 30,
                         resampling = "validation",
                         metric = "auc", maximize = TRUE),
  threshold_method = "youden",
  fit_final      = TRUE,
  evaluate_test  = TRUE,
  seed           = 1234
)

flow$final_test_metrics                           # aggregate — shareable
# flow$final_test_predictions                     # keep on RAP

ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig13-roc.pdf",
                plot_ml_roc(flow$final_model, newdata = flow$split$test))
```

**Note on `use_validation_in_refit`:** `ukb_ml_fit_final()` uses the
`use_validation_in_refit` argument to control whether validation-set rows
are folded into the final refit. This defaults to `TRUE` in
`ukb_ml_fit_final()` itself. In high-level frozen-test workflows such as
`ukb_ml_flow()` and `ukb_ml_compare_flows()`, check the function default and
pass this argument explicitly when you need to include or exclude validation
rows from the final model.

### 4b. Survival — `ukb_ml_survival_workflow()`

```r
sflow <- ukb_ml_survival_workflow(
  formula = survival::Surv(outcome_surv_time, outcome_status) ~
              age + sex + bmi + sbp + ldl_cholesterol + hba1c,
  data    = cohort,
  model   = "rsf",                          # "rsf", "gbm_surv", "coxnet"
  split_params = list(split = "train_valid_test",
                       train_ratio = 0.7, validation_ratio = 0.1,
                       test_ratio = 0.2,
                       stratify_by = "event"),
  tune    = TRUE,
  tune_params = list(search = "grid", metric = "c_index"),
  fit_final     = TRUE,
  evaluate_test = TRUE,
  seed   = 1234
)
sflow$final_test_metrics      # c_index + (optional) time-dependent IBS
sflow$final_test_predictions  # per-row survival probabilities — RAP-only
```

Predict at specific horizons:

```r
preds <- ukb_ml_survival_predict(
  sflow,
  newdata = sflow$split$test,
  times   = c(1, 3, 5, 10),
  type    = "survival"        # "survival", "risk", "chf"
)
```

---

## 5. Sub-steps (when you don't want the all-in-one workflow)

```r
# 1. Frozen split
split <- ukb_ml_split_data(
  cohort, outcome = "outcome_status",
  outcome_type = "binary",
  split = "train_valid_test",
  train_ratio = 0.7, validation_ratio = 0.1, test_ratio = 0.2,
  stratify_by = "outcome", seed = 1234
)

# 2. Feature selection
fs <- ukb_ml_feature_select(split, formula = outcome_status ~ .,
                             method = "boruta", boruta_params = list(),
                             seed = 1234)
formula_fs <- fs$formula

# 3. Tuning
tune <- ukb_ml_tune(split, formula = formula_fs, model = "xgboost",
                     search = "random", n_iter = 30,
                     resampling = "validation", metric = "auc",
                     seed = 1234)
best <- tune$best_params

# 4. Final fit
final <- ukb_ml_fit_final(split, formula = formula_fs, model = "xgboost",
                           best_params = best, use_validation_in_refit = TRUE,
                           seed = 1234)

# 5. Test evaluation
eval <- ukb_ml_evaluate_test(final, split, threshold = NULL)
eval$metrics
```

Use this sub-stepped flow when the user needs **per-stage artefacts**
(e.g. Boruta status table, tuning grid history) for the manuscript.

---

## 6. Comparison across models

```r
cmp <- ukb_ml_compare_flows(
  formula   = outcome_status ~ .,
  data      = cohort,
  models    = c("rf","xgboost","glmnet","logistic"),
  outcome_type = "binary",
  split_params = list(split = "train_test", train_ratio = 0.7,
                       stratify_by = "outcome"),
  tune_params  = list(search = "random", n_iter = 30),
  feature_params = list(method = "none"),
  n_repeats = 1, seed = 1234
)
# Plot overlay
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig16-roc_compare.pdf",
                plot_ml_roc_compare(cmp))
```

Alternative: compare *feature sets* with one model:

```r
cmp_fs <- ukb_ml_compare_feature_sets(
  split, formula = outcome_status ~ .,
  models = "xgboost",
  feature_selections = list(none = NULL, boruta = list(),
                             filter = list(top_k = 20))
)
```

---

## 7. SHAP interpretation

```r
shap <- ukb_shap(
  flow,                       # ukb_ml_flow or ukb_ml_workflow object
  data    = flow$split$test,
  nsim    = 100,
  sample_n = 2000,            # random subsample for speed
  seed    = 1234,
  method  = "auto"            # "auto" picks native XGBoost SHAP for tree models
)
```

**Sharing rules for SHAP outputs:**

| SHAP output | Sharing rule |
|-------------|-------------|
| Importance table (mean \|SHAP\| per feature) | Shareable — aggregate summary, no identifiers. |
| Summary / beeswarm / dependence plots | Shareable as rendered figures only; do **not** attach the row-level SHAP source matrix or feature-value table. |
| Force plot for a real participant row | Do not share with the local agent. Use a synthetic or representative prototype — a row constructed from group medians, not an actual participant. |
| Raw `shap$shap_values` matrix | Keep on RAP and out of the agent context; it has one row per observation, even if `eid` is absent. |

```r
# Summary table is aggregate; beeswarm is a rendered figure. Do not export
# the row-level SHAP matrix or feature-value table.
fwrite(ukb_shap_summary(shap, top_n = 20, plot = FALSE),
       "/mnt/project/<area>/04-results/08-ml_shap_summary.csv")
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig17-shap_summary.pdf",
                plot_shap_summary(shap, top_n = 15))
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig18-shap_beeswarm.pdf",
                plot_shap_beeswarm(shap, top_n = 15))

# Dependence plot — figure only, no raw source CSV export
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig18b-shap_dep_sbp.pdf",
                ukb_shap_dependence(shap, feature = "sbp", color_by = "age"))

# Force plot — use a SYNTHETIC prototype, not a real participant row
prototype <- cohort[, lapply(.SD, median, na.rm = TRUE),
                    .SDcols = names(shap$feature_values)]
shap_proto <- ukb_shap(flow, data = prototype, nsim = 50, seed = 1234)
ukb_shap_force(shap_proto, obs_id = 1)
```

For survival ML, use `ukb_ml_survival_shap()` (same return structure,
computed via `fastshap` internally).

---

## 8. Threshold optimization & decision-curve analysis

```r
thr <- ukb_ml_threshold(
  truth = flow$split$test[[flow$split$outcome]],
  prob  = predict(flow$final_model, newdata = flow$split$test, type = "prob")[,2],
  method = "youden"
)
thr$threshold

dca <- ukb_ml_dca(flow$final_model, newdata = flow$split$test,
                   threshold = NULL, positive_class = NULL)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig19-dca.pdf",
                plot_ml_dca(flow$final_model, newdata = flow$split$test))
```

---

## 9. Common pitfalls

1. **Leakage at split.** Always use the same `ukb_ml_split` object across
   feature selection / tuning / final fit / evaluation. Re-splitting between
   stages leaks the test set.
2. **`outcome_type = "auto"`** misfires when the target is an integer-coded
   binary with three distinct values (`0, 1, NA`). Set explicitly.
3. **Class imbalance.** XGBoost ignores it by default. Either set
   `scale_pos_weight` via `param_grid` or stratify the split via
   `stratify_by = "outcome"`.
4. **Threshold of 0.5 is rarely optimal.** Use `threshold_method = "youden"`
   in the workflow or `ukb_ml_threshold()` post-hoc.
5. **Survival outcome formula.** Must use `survival::Surv(time, event) ~ …`.
   Bare `~ outcome` is a classification model.
6. **`use_validation_in_refit` differs by layer.** Low-level
   `ukb_ml_fit_final()` defaults to refitting on train + validation, but
   high-level frozen-test workflows may expose a different default. Check the
   function formals and set the argument explicitly when validation rows must
   remain a separate evaluation set.
7. **SHAP cost.** `nsim = 100` is fine; pushing to 1000 multiplies cost.
   For ~50k test rows use `sample_n = 2000` (random subsample is honest for
   the global summary).
8. **Force plots:** do not share a force plot built from a real participant's
   raw feature values with the local agent. Construct a synthetic prototype
   (e.g. group-median row) for any illustrative force-plot figure. The SHAP
   source matrix (`shap$shap_values`) has one row per observation and stays
   inside RAP.
9. **`fastshap` vs native XGBoost.** With `method = "auto"`, the function
   picks the native XGBoost SHAP for tree models (fast, exact). `fastshap`
   handles everything else but is slower.

---

## 10. Key functions (selected)

| Layer | Function |
|-------|----------|
| Models supported | `ukb_ml_supported_models(outcome_type)` |
| End-to-end (classification) | `ukb_ml_workflow(formula, data, ..., model, outcome_type, split_params, feature_select, tune, threshold_method, fit_final, evaluate_test, seed, verbose)` |
| End-to-end flexible | `ukb_ml_flow(formula = NULL, data = NULL, split = NULL, train_data = NULL, test_data = NULL, ..., model, model_id, outcome_type, tune, threshold_method, compute_shap, shap_data, shap_params, seed, verbose)` |
| End-to-end (survival) | `ukb_ml_survival_workflow(formula, data, ..., model, split_params, tune_params, fit_final, evaluate_test, seed, verbose)` |
| Split | `ukb_ml_split_data(...)`, `ukb_ml_as_split(train_data, test_data, validation_data = NULL, id_col)` |
| Survival split | `ukb_ml_survival_split_data(...)`, `ukb_ml_survival_as_split(...)` |
| Feature selection | `ukb_ml_feature_select(split, formula, method, ...)` / survival: `ukb_ml_survival_feature_select(...)` |
| Tuning | `ukb_ml_tune(split, formula, model, search, param_grid, ...)` / survival: `ukb_ml_survival_tune(...)` |
| Fit final | `ukb_ml_fit_final(...)` / survival: `ukb_ml_survival_fit_final(...)` |
| Predict | `ukb_ml_predict(object, newdata, type)` / survival: `ukb_ml_survival_predict(object, newdata, times, type)` |
| Evaluate | `ukb_ml_evaluate_test(object, split, metrics, threshold, positive_class)` / survival: `ukb_ml_survival_evaluate_test(...)` |
| Single-metric helpers | `ukb_ml_roc, ukb_ml_pr, ukb_ml_calibration, ukb_ml_dca, ukb_ml_confusion, ukb_ml_threshold, ukb_ml_gain_lift, ukb_ml_ks` |
| Compare | `ukb_ml_compare(object1, object2, ...)`, `ukb_ml_compare_flows(...)`, `ukb_ml_compare_feature_sets(...)` |
| Importance | `ukb_ml_importance(object, type)` / survival: `ukb_ml_survival_importance(object)` |
| SHAP | `ukb_shap(object, data, nsim, sample_n, seed, method)` / `ukb_shap_summary`, `ukb_shap_dependence`, `ukb_shap_force`; survival: `ukb_ml_survival_shap(...)` |
| Plot ML | `plot_ml_roc`, `plot_ml_roc_compare`, `plot_ml_pr`, `plot_ml_calibration`, `plot_ml_confusion`, `plot_ml_dca`, `plot_ml_gain`, `plot_ml_lift`, `plot_ml_ks`, `plot_ml_importance`, `plot_ml_compare` |
| Plot SHAP | `plot_shap_summary`, `plot_shap_beeswarm`, `plot_shap_dependence`, `plot_shap_force` |

`ukb_ml_supported_models()` returns: `logistic, linear, rf (ranger), xgboost,
glmnet, svm (e1071), nnet, rpart, naive_bayes`. For survival:
`rsf (randomForestSRC), gbm_surv (gbm), coxnet (glmnet)`.

---

## 11. Related skills

| Skill | When |
|-------|------|
| `ukbsci-cohort` | Build outcomes first. |
| `ukbsci-preprocess` | Feature engineering, missing-code handling. |
| `ukbsci-imputation` | Pool ML metrics across imputations (custom). |
| `ukbsci-plot` | Manuscript-ready figure assembly. |
| `ukbsci-workflow` | Phase 7 — ML model construction & evaluation. |

---

## 12. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md) — binary XGBoost flow,
  RSF survival flow, multi-model comparison + SHAP summary
- [`references/models-and-tuning.md`](references/models-and-tuning.md) —
  per-model tunable parameter cheat-sheet
