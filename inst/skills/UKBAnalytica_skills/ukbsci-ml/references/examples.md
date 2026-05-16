# `ukbsci-ml` — examples

```r
library(UKBAnalytica); library(data.table); library(survival)
```

## A. Binary XGBoost workflow with SHAP

```r
flow <- ukb_ml_workflow(
  formula = outcome_status ~ age + sex + bmi + sbp + ldl_cholesterol +
                              hba1c + crp + smoking_status,
  data    = cohort,
  model   = "xgboost", outcome_type = "binary",
  split_params = list(split = "train_valid_test",
                       train_ratio = 0.7, validation_ratio = 0.1,
                       test_ratio = 0.2, stratify_by = "outcome"),
  feature_select = "boruta",
  tune    = TRUE,
  tune_params = list(search = "random", n_iter = 30,
                      resampling = "validation",
                      metric = "auc", maximize = TRUE),
  threshold_method = "youden",
  fit_final = TRUE, evaluate_test = TRUE,
  seed = 1234
)

# Aggregate metrics — safe to export
fwrite(as.data.frame(t(flow$final_test_metrics)),
       "/mnt/project/<area>/04-results/07-ml_metrics_xgb.csv")

ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig13-roc.pdf",
                plot_ml_roc(flow$final_model, newdata = flow$split$test))
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig14-calib.pdf",
                plot_ml_calibration(flow$final_model, newdata = flow$split$test))
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig15-dca.pdf",
                plot_ml_dca(flow$final_model, newdata = flow$split$test))

# Global SHAP
shap <- ukb_shap(flow, data = flow$split$test, nsim = 100, sample_n = 2000,
                  seed = 1234, method = "auto")
imp <- ukb_shap_summary(shap, top_n = 15, plot = FALSE)
fwrite(imp, "/mnt/project/<area>/04-results/08-ml_shap_summary.csv")
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig17-shap_summary.pdf",
                plot_shap_summary(shap, top_n = 15))
```

## B. Random Survival Forest

```r
sflow <- ukb_ml_survival_workflow(
  formula = Surv(outcome_surv_time, outcome_status) ~ age + sex + bmi + sbp +
              ldl_cholesterol + hba1c + crp + smoking_status,
  data = cohort, model = "rsf",
  split_params = list(split = "train_valid_test",
                       train_ratio = 0.7, validation_ratio = 0.1,
                       test_ratio = 0.2, stratify_by = "event"),
  tune        = TRUE,
  tune_params = list(search = "grid",
                      param_grid = expand.grid(
                        ntree = c(500, 1000),
                        mtry  = c(3, 5, 8),
                        nodesize = c(5, 15, 30)),
                      metric = "c_index"),
  fit_final = TRUE, evaluate_test = TRUE,
  seed = 1234
)
sflow$final_test_metrics$c_index

preds <- ukb_ml_survival_predict(sflow,
            newdata = sflow$split$test,
            times   = c(1, 3, 5, 10), type = "survival")
```

## C. Multi-model comparison

```r
cmp <- ukb_ml_compare_flows(
  formula = outcome_status ~ .,
  data    = cohort,
  models  = c("logistic","glmnet","rf","xgboost"),
  outcome_type = "binary",
  split_params = list(split = "train_test", train_ratio = 0.7,
                       stratify_by = "outcome"),
  tune_params  = list(search = "random", n_iter = 20),
  feature_params = list(method = "none"),
  seed = 1234
)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig16-roc_compare.pdf",
                plot_ml_roc_compare(cmp))
```
