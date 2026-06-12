# `ukbsci-ml` — function reference

Reflects `UKBAnalytica` 1.0.0 (`R/ml_workflow.R`, `R/ml_model.R`,
`R/ml_evaluate.R`, `R/ml_shap.R`, `R/ml_survival.R`).

---

## End-to-end workflows

```r
ukb_ml_workflow(formula, data = NULL, split = NULL, model,
                outcome_type     = c("auto","binary","multiclass","continuous"),
                split_params     = list(),
                feature_select   = c("none","boruta","filter","glmnet"),
                feature_params   = list(),
                tune             = TRUE, tune_params = list(),
                threshold_method = c("none","fixed","youden"),
                threshold_params = list(),
                fit_final        = TRUE,
                evaluate_test    = TRUE,
                seed             = NULL,
                verbose          = TRUE, ...)

ukb_ml_flow(formula = NULL, data = NULL, split = NULL,
            train_data = NULL, test_data = NULL, validation_data = NULL,
            id_col = NULL, outcome = NULL, features = NULL,
            model = "xgboost", model_id = "model", model_label = NULL,
            outcome_type = c("auto","binary","multiclass","continuous"),
            split_params = list(), param_grid = NULL,
            tune = TRUE, tune_params = list(), best_params = NULL,
            threshold_method = c("none","fixed","youden"),
            threshold_params = list(),
            metrics = NULL, positive_class = NULL,
            use_validation_in_refit = FALSE,
            compute_shap = FALSE, shap_data = NULL, shap_params = list(),
            seed = NULL, verbose = TRUE)

ukb_ml_survival_workflow(formula, data = NULL, split = NULL,
                         model = c("rsf","gbm_surv","coxnet"),
                         outcome_type = "survival",
                         split_params = list(),
                         tune = TRUE, tune_params = list(),
                         fit_final = TRUE, evaluate_test = TRUE,
                         seed = NULL, verbose = TRUE, ...)
```

**Returns:** S3 `ukb_ml_workflow` / `ukb_ml_flow` / `ukb_ml_survival_workflow`
containing `$split, $feature_result, $tune_result, $threshold_result,
$final_model, $final_test_metrics, $final_test_predictions` (and `$shap`
for `ukb_ml_flow(compute_shap = TRUE)`).

---

## Split helpers

```r
ukb_ml_split_data(df, outcome = NULL,
                  outcome_type = c("auto","binary","multiclass","continuous"),
                  split = c("train_test","train_valid_test"),
                  train_ratio = 0.7, validation_ratio = 0.1,
                  test_ratio = 0.2, split_ratio = NULL,
                  stratify_by = c("auto","outcome","custom","none"),
                  stratify_col = NULL, regression_bins = 5,
                  seed = NULL, verbose = TRUE)

ukb_ml_as_split(train_data, test_data, validation_data = NULL,
                id_col = NULL, check_overlap = TRUE,
                outcome = NULL,
                outcome_type = c("auto","binary","multiclass","continuous"))

ukb_ml_survival_split_data(df, time_var, event_var,
                           split = c("train_test","train_valid_test"),
                           train_ratio = 0.7, validation_ratio = 0.1,
                           test_ratio = 0.2,
                           stratify_by = c("none","event"),
                           seed = NULL, verbose = TRUE)

ukb_ml_survival_as_split(train_data, test_data, validation_data = NULL,
                         time_var, event_var, check_overlap = FALSE)
```

`ukb_ml_split` is the canonical S3 token threaded through every other step.

---

## Feature selection

```r
ukb_ml_feature_select(split, formula,
                      method = c("none","boruta","filter","glmnet"),
                      outcome_type = c("auto","binary","multiclass","continuous"),
                      max_features = NULL,
                      boruta_params = list(), keep_tentative = TRUE,
                      seed = NULL, verbose = TRUE)

ukb_ml_survival_feature_select(split, formula,
                               method = c("none","boruta","filter","glmnet"),
                               boruta_params = list(),
                               seed = NULL, verbose = TRUE)
```

**Returns:** `list(formula, selected_features, status, info)`.

---

## Tuning

```r
ukb_ml_tune(split, formula, model,
            outcome_type = c("auto","binary","multiclass","continuous"),
            search       = c("grid","random","bayes"),
            param_grid   = NULL, param_space = NULL,
            n_iter       = NULL,
            resampling   = c("cv","validation"), folds = 5,
            metric = NULL, maximize = NULL,
            seed = NULL, verbose = TRUE, ...)

ukb_ml_survival_tune(split, formula,
                     model = c("rsf","gbm_surv","coxnet"),
                     search = c("grid","random"),
                     param_grid = NULL,
                     metric = "c_index",
                     seed = NULL, verbose = TRUE, ...)
```

`search = "bayes"` requires `rBayesianOptimization` (falls back to random
when missing).

---

## Final fit + prediction

```r
ukb_ml_fit_final(split, formula, model, best_params = list(),
                 outcome_type = c("auto","binary","multiclass","continuous"),
                 feature_spec = NULL, threshold = NULL,
                 use_validation_in_refit = TRUE,
                 seed = NULL, verbose = TRUE, ...)

ukb_ml_survival_fit_final(split, formula, model, best_params = list(),
                          seed = NULL, verbose = TRUE, ...)

ukb_ml_predict(object, newdata = NULL, type = c("response","prob"))
ukb_ml_survival_predict(object, newdata = NULL,
                        times = c(1, 3, 5, 10),
                        type  = c("survival","risk","chf"), ...)
```

---

## Test-set evaluation

```r
ukb_ml_evaluate_test(object, split, metrics = NULL,
                     threshold = NULL, positive_class = NULL,
                     verbose = TRUE)

ukb_ml_survival_evaluate_test(object, split,
                              times = c(1, 3, 5, 10),
                              verbose = TRUE)
```

Returns `list(metrics, predictions)` (predictions are participant-level —
RAP-only).

---

## Single-metric helpers

```r
ukb_ml_roc(object, newdata = NULL)
ukb_ml_roc_data(truth, prob, model_id = NULL, model_label = NULL,
                positive_class = NULL, ci = TRUE,
                ci_method = c("delong","bootstrap"), quiet = TRUE)
ukb_ml_pr(object, newdata = NULL)
ukb_ml_calibration(object, newdata = NULL, n_bins = 10)
ukb_ml_confusion(object, newdata = NULL, threshold = NULL,
                 positive_class = NULL)
ukb_ml_dca(object, newdata = NULL, threshold = NULL,
           positive_class = NULL)
ukb_ml_gain_lift(object, newdata = NULL, threshold = NULL,
                 positive_class = NULL, n_bins = 10)
ukb_ml_ks(object, newdata = NULL)
ukb_ml_threshold(truth, prob,
                 method = c("fixed","youden"),
                 fixed_threshold = 0.5, positive_class = NULL)
ukb_ml_metrics(object, newdata = NULL, metrics = NULL, ci = FALSE,
               ci_method = c("bootstrap","delong"),
               n_boot = 1000, verbose = TRUE, ...)
```

---

## Comparison

```r
ukb_ml_compare(object1, object2, ..., metric = NULL,
               positive_class = NULL)

ukb_ml_compare_flows(formula = NULL, data = NULL,
                     models = c("rf","xgboost"),
                     outcome_type = c("auto","binary","multiclass","continuous"),
                     split_params = list(), tune_params = list(),
                     feature_params = list(),
                     n_repeats = 1, seed = NULL, verbose = TRUE)

ukb_ml_compare_feature_sets(split, formula,
                             models = c("rf","xgboost"),
                             feature_selections = list(none = NULL,
                                                        boruta = list(),
                                                        filter = list()),
                             outcome_type = c("auto","binary","multiclass","continuous"),
                             seed = NULL, verbose = TRUE)
```

---

## Importance + SHAP

```r
ukb_ml_importance(object, type = c("model","shap"),
                   newdata = NULL, shap_params = list())
ukb_ml_survival_importance(object, ...)

ukb_shap(object, data = NULL, nsim = 100, sample_n = NULL,
         seed = NULL, verbose = TRUE,
         class_level = NULL,
         method = c("auto","fastshap","xgboost"), ...)

ukb_shap_summary(shap_object, top_n = 10, plot = TRUE)
ukb_shap_dependence(shap_object, feature, color_by = NULL, plot = TRUE)
ukb_shap_force(shap_object, obs_id = 1, plot = TRUE)

ukb_ml_survival_shap(object, data = NULL, nsim = 100, sample_n = NULL,
                     seed = NULL, verbose = TRUE, ...)
```

---

## Visualization functions (return ggplot)

`plot_ml_roc, plot_ml_roc_compare, plot_ml_pr, plot_ml_calibration,
plot_ml_confusion, plot_ml_dca, plot_ml_gain, plot_ml_lift, plot_ml_ks,
plot_ml_importance, plot_ml_compare`.

`plot_shap_summary, plot_shap_beeswarm, plot_shap_dependence,
plot_shap_force`.

---

## `ukb_ml_supported_models()`

```r
ukb_ml_supported_models(outcome_type = c("all","binary","multiclass","continuous"))
```

Returns: data.frame with columns `model, label, outcome_types, package,
default_tuning`. See `models-and-tuning.md` for the table.
