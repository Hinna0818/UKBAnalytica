# `ukbsci-propensity` — examples

```r
library(UKBAnalytica); library(data.table); library(survival)
covars_ps <- c("age","sex","bmi","smoking_status",
               "Hypertension_history","Diabetes_history",
               "ldl_cholesterol","sbp")
```

## A. PSM with caliper + exact-match on sex

```r
ps_data <- estimate_propensity_score(cohort, "treatment", covars_ps,
                                     method = "logistic")

matched <- match_propensity(
  data = ps_data, ps_col = "ps", treatment = "treatment",
  ratio = 1, caliper = 0.2, method = "nearest",
  replace = FALSE, exact_match = "sex"
)

b_pre  <- assess_balance(cohort,  "treatment", covars_ps, method = "unmatched")
b_post <- assess_balance(matched, "treatment", covars_ps, method = "matched")
plot_balance(b_pre, b_post, title = "Love plot — PSM")
```

## B. IPTW with stabilized ATE weights, truncation

```r
weighted <- calculate_weights(
  data = ps_data, ps_col = "ps", treatment = "treatment",
  weight_type = "ATE", stabilized = TRUE,
  truncate = c(0.01, 0.99)
)
b_w <- assess_balance(weighted, "treatment", covars_ps,
                      method = "weighted", weight_col = "weight")
plot_ps_distribution(ps_data, "ps", "treatment", type = "mirror")
```

## C. Weighted Cox with robust SE

```r
wt_cox <- run_weighted_analysis(
  data = weighted, exposure = "treatment",
  weight_col = "weight",
  model_type = "cox",
  endpoint   = c("outcome_surv_time", "outcome_status"),
  robust_se  = TRUE
)
fwrite(wt_cox, "/mnt/project/<area>/04-results/02d-cox_iptw.csv")
```
