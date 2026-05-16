# `ukbsci-subgroup-sensitivity` — examples

```r
library(UKBAnalytica); library(survival); library(data.table)
covars <- c("age","sex","bmi","smoking_status")
```

## A. Single-modifier subgroup analysis

```r
sub <- run_subgroup_analysis(
  data = cohort, exposure = "sbp",
  subgroup_var = "sex", covariates = setdiff(covars, "sex"),
  model_type = "cox",
  endpoint = c("outcome_surv_time","outcome_status"),
  ref_level = "Male"
)
print(sub)
```

## B. Many modifiers at once

```r
multi_sub <- run_multi_subgroup(
  data = cohort, exposure = "sbp",
  subgroup_vars = c("sex","age_band","smoking_status","Diabetes_history"),
  covariates = c("age","bmi"),
  model_type = "cox",
  endpoint = c("outcome_surv_time","outcome_status")
)
fwrite(multi_sub, "/mnt/project/<area>/04-results/03-cox_subgroup.csv")
```

## C. Lag sensitivity (drop early events)

```r
trim <- sensitivity_exclude_early_events(
  data = cohort,
  endpoint = c("outcome_surv_time","outcome_status"),
  n_years = 2
)
attr(trim, "sensitivity_info")

cox_lag2 <- runmulti_cox(trim, "sbp", covars,
                         endpoint = c("outcome_surv_time","outcome_status"))
fwrite(cox_lag2,
       "/mnt/project/<area>/04-results/04a-cox_lag2.csv")
```

## D. Complete-case + flow trace

```r
cc <- sensitivity_exclude_missing_covariates(
  data = cohort,
  covariates = c(covars, "ldl_cholesterol", "hba1c"),
  stepwise = TRUE
)
flow <- attr(cc, "complete_case_flow")
print(flow)

cox_cc <- runmulti_cox(cc, "sbp", covars,
                       endpoint = c("outcome_surv_time","outcome_status"))
fwrite(cox_cc, "/mnt/project/<area>/04-results/04b-cox_complete_case.csv")
```
