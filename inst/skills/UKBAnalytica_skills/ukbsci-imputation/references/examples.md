# `ukbsci-imputation` — examples

```r
library(UKBAnalytica); library(mice); library(mitools); library(survival); library(data.table)
```

## A. Standard 5-imputation pooled Cox

```r
imp <- run_imputation(
  data   = cohort, id_col = "eid",
  vars   = c("bmi","ldl_cholesterol","hba1c","crp","sbp","smoking_status"),
  factor_vars = "smoking_status",
  method = "pmm", m = 5, maxit = 10, seed = 1234
)

models <- fit_mi_models(
  datasets = imp$data_list,
  formula = Surv(outcome_surv_time, outcome_status) ~ sbp + age + sex + bmi + smoking_status,
  model_type = "cox"
)

pooled <- pool_mi_models(
  models = models,
  formula = Surv(outcome_surv_time, outcome_status) ~ sbp + age + sex + bmi + smoking_status,
  model_type = "cox", exponentiate = TRUE
)
print(pooled)
fwrite(tidy(pooled, exponentiate = TRUE),
       "/mnt/project/<area>/04-results/02e-cox_pooled.csv")
```

## B. Pool a mediation TNIE across imputations

```r
ests <- list(); vcovs <- list()
for (d in imp$data_list) {
  r <- run_mediation(d,
    exposure = "sbp_above_median", mediator = "ldl_cholesterol",
    outcome  = "outcome_surv_time",
    endpoint = c("outcome_surv_time","outcome_status"),
    mediator_type = "continuous", outcome_type = "cox")
  ests[[length(ests)+1]]  <- coef(r)
  vcovs[[length(vcovs)+1]] <- vcov(r$regmedint_obj)
}
pooled_med <- pool_custom_estimates(
  estimates = ests, variances = vcovs,
  labels    = c("CDE","PNDE","TNIE","TNDE","PNIE","TE","PM")
)
summary(pooled_med)
```

## C. Diagnostic figures

```r
p_pool <- plot_mi_pooled(pooled, terms = c("sbp","bmi","smoking_status"),
                         exponentiate = TRUE, null_value = 1, show_fmi = TRUE,
                         title = "Pooled Cox HRs (5 imputations)")
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig10-mi_pooled_forest.pdf",
                p_pool, width = 8, height = 5)

ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig10b-mi_fmi.pdf",
                plot_mi_diagnostics(pooled, type = "fmi"), width = 7, height = 4)
```
