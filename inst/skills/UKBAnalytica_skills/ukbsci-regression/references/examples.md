# `ukbsci-regression` — examples

Standard header for every example:

```r
library(UKBAnalytica)
library(data.table)
library(survival)

stopifnot(data.table::is.data.table(cohort))
stopifnot(all(c("outcome_status", "outcome_surv_time") %in% colnames(cohort)))

covars <- c("age", "sex", "bmi", "smoking_status",
            "Hypertension_history", "Diabetes_history")
```

---

## A. Cox model — one outcome, many exposures

```r
exposures <- c("sbp", "ldl_cholesterol", "hba1c", "crp")

cox_res <- runmulti_cox(
  data       = cohort,
  main_var   = exposures,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status")
)
print(cox_res)
fwrite(cox_res, "/mnt/project/<area>/04-results/02-cox_main.csv")
```

## B. PH diagnostics + lagged sensitivity

```r
zph <- runmulti_cox_zph(
  data       = cohort,
  main_var   = exposures,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status"),
  transform  = "km", alpha = 0.05
)
problematic <- zph[ph_violation == TRUE, variable]

# Re-run a lagged-sensitivity Cox on the problematic exposures only
lag_res <- runmulti_cox_lag(
  data       = cohort,
  main_var   = problematic,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status"),
  lag_years  = c(0, 2, 5)
)
fwrite(lag_res, "/mnt/project/<area>/04-results/04-cox_lag.csv")
```

## C. Fine-Gray competing risks

```r
# Suppose cohort has `outcome_status_full` with {0=censored, 1=AA, 2=non-AA death}
fg <- runmulti_competing(
  data           = cohort,
  main_var       = exposures,
  covariates     = covars,
  time_col       = "outcome_surv_time",
  event_col      = "outcome_status_full",
  event_value    = 1,
  compete_value  = 2
)
fwrite(fg, "/mnt/project/<area>/04-results/02c-fg_competing.csv")
```

## D. Ordered-trend across quartiles

```r
cohort[, sbp_q := cut(sbp,
                      quantile(sbp, c(0, .25, .5, .75, 1), na.rm = TRUE),
                      labels = paste0("Q", 1:4),
                      include.lowest = TRUE)]

trend_res <- runmulti_trend(
  data            = cohort,
  main_var        = "sbp_q",
  covariates      = covars,
  model_type      = "cox",
  endpoint        = c("outcome_surv_time", "outcome_status"),
  ref_level       = "Q1",
  score_method    = "integer",
  include_level_estimates = TRUE
)
trend_res
```

## E. Single-model PH diagnostics

```r
fit <- coxph(Surv(outcome_surv_time, outcome_status) ~ sbp + age + sex,
             data = cohort)
diag <- ukb_cox_diagnostics(fit, transform = "km", alpha = 0.05)
diag$table         # per-term zph
diag$global_pvalue # global PH test
```
