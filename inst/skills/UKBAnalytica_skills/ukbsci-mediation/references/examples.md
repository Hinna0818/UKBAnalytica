# `ukbsci-mediation` — examples

```r
library(UKBAnalytica); library(regmedint); library(data.table)
covars <- c("age","sex","bmi","smoking_status")
```

## A. SBP → LDL → incident AA (Cox)

```r
res <- run_mediation(
  data = cohort,
  exposure = "sbp_above_median", mediator = "ldl_cholesterol",
  outcome  = "outcome_surv_time", covariates = covars,
  mediator_type = "continuous", outcome_type = "cox",
  endpoint = c("outcome_surv_time","outcome_status"),
  interaction = TRUE, boot = FALSE
)
summary(res, exponentiate = TRUE)
res$effects
```

## B. Screen 4 candidate mediators

```r
multi <- run_multi_mediator(
  data = cohort,
  exposure = "sbp_above_median",
  mediators = c("ldl_cholesterol","hba1c","crp","bmi"),
  outcome = "outcome_surv_time", covariates = covars,
  mediator_type = "continuous", outcome_type = "cox",
  endpoint = c("outcome_surv_time","outcome_status")
)
multi[order(multi$tnie_p), ]
fwrite(multi, "/mnt/project/<area>/04-results/06-mediation_screen.csv")
```

## C. Sensitivity + forest plot

```r
sens <- run_sensitivity_mediation(res, rho_values = seq(-0.9, 0.9, 0.1))
plot_mediation(res, type = "effects", exponentiate = TRUE)

p_forest <- plot_mediation_forest(multi, effect_type = "tnie", null_value = 0)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig09-mediation_forest.pdf",
                p_forest, width = 8, height = 5)
```
