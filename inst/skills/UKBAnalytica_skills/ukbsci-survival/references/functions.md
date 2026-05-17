# `ukbsci-survival` — function reference

## `plot_km_curve()`

```r
plot_km_curve(data, time_col, status_col,
              group_col    = NULL,
              conf_int     = TRUE,
              risk_table   = TRUE,
              censor_marks = TRUE,
              palette      = "jco",
              title        = NULL,
              xlab         = "Time (years)",
              ylab         = "Survival Probability",
              legend_title = "Group",
              median_line  = TRUE,
              pvalue       = TRUE,
              xlim         = NULL,
              break_time   = NULL)
```

| Arg | Type | Meaning |
|-----|------|---------|
| `data` | data.frame / data.table | survival dataset (from `build_survival_dataset()`) |
| `time_col` | char | follow-up time column (years) |
| `status_col` | char | event indicator (1 = event, 0 = censored) |
| `group_col` | char or `NULL` | optional stratifier |
| `palette` | char | `"jco"`, `"nejm"`, `"lancet"`, or hex vector |
| `pvalue` | logical | log-rank p when `group_col` present |
| `risk_table` | logical | append numbers-at-risk panel |

**Returns:**

| `risk_table = TRUE` | `risk_table = FALSE` |
|---|---|
| `list($km_plot, $risk_table)` from `survminer::ggsurvplot` | bare `ggplot2` object |

Requires `survival` + `survminer`.

**Behavior:**

- Drops `NA` time or status rows automatically (via `survfit`).
- Drops `NA` group rows when `group_col` is supplied.
- Calculates log-rank p only when `group_col` is provided.
- Median line drawn only when curves cross 0.5.

---

## Building the figure-source CSV (the skill always emits this alongside the figure)

```r
fit <- survival::survfit(
  survival::Surv(outcome_surv_time, outcome_status) ~ group_var,
  data = cohort
)
km_data <- data.table::data.table(
  time = fit$time, n.risk = fit$n.risk, n.event = fit$n.event,
  n.censor = fit$n.censor, surv = fit$surv,
  lower = fit$lower, upper = fit$upper,
  strata = rep(names(fit$strata), fit$strata)
)
```

The CSV is aggregate (counts + Kaplan-Meier point estimates) and may be shared
with the local agent.
