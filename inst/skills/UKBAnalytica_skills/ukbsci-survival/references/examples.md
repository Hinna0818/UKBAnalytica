# `ukbsci-survival` — examples

```r
library(UKBAnalytica); library(survival); library(survminer); library(data.table)
```

## A. Overall KM

```r
p <- plot_km_curve(
  data = cohort, time_col = "outcome_surv_time",
  status_col = "outcome_status",
  conf_int = TRUE, risk_table = TRUE,
  title = "Overall survival",
  xlab  = "Follow-up (years)",
  ylab  = "Event-free probability"
)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig03-km_overall.pdf", p, width=7, height=6)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig03-km_overall.png", p, width=7, height=6, dpi=300)
```

## B. Stratified KM with log-rank p

```r
p <- plot_km_curve(
  data         = cohort,
  time_col     = "outcome_surv_time",
  status_col   = "outcome_status",
  group_col    = "Hypertension_history",
  conf_int     = TRUE,
  risk_table   = TRUE,
  palette      = "jco",
  pvalue       = TRUE,
  median_line  = TRUE,
  xlim         = c(0, 15),
  break_time   = 3,
  legend_title = "Baseline HTN"
)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig03b-km_htn.pdf", p, width=7, height=6)
```

## C. Save figure source (always)

```r
fit <- survfit(Surv(outcome_surv_time, outcome_status) ~ Hypertension_history,
               data = cohort)
km_data <- data.table::data.table(
  time = fit$time, n.risk = fit$n.risk, n.event = fit$n.event,
  n.censor = fit$n.censor, surv = fit$surv,
  lower = fit$lower, upper = fit$upper,
  strata = rep(names(fit$strata), fit$strata)
)
fwrite(km_data, "/mnt/project/<area>/05-figs/data/Fig03b-km_htn.csv")
```

## D. Custom palette + manuscript layout

```r
my_colors <- c("#2F6FA3", "#C74732")    # ukbsci_clinical baseline + accent

p <- plot_km_curve(
  data = cohort, time_col = "outcome_surv_time",
  status_col = "outcome_status",
  group_col = "exposure_group",
  palette = my_colors,
  conf_int = TRUE, risk_table = TRUE,
  pvalue = TRUE,
  xlim = c(0, 12), break_time = 2
)
```
