---
name: ukbsci-survival
description: >
  Survival visualization and post-Cox summary helpers for UK Biobank
  Research Analysis Platform (RAP) cohorts built with UKBAnalytica.
  Covers Kaplan-Meier curves with risk tables and log-rank p-values
  (plot_km_curve), and pairs with ukbsci-regression for time-to-event
  analyses. Use this skill when the user asks to draw a KM curve, show
  cumulative incidence, build a risk-table-augmented survival figure, or
  visualize stratified survival across exposure / treatment groups.
  Triggers: UKB KM curve, Kaplan-Meier, survival curve, log-rank, risk
  table, cumulative incidence, UKB 生存曲线, KM 曲线, /ukbsci-survival.
  Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-survival — Survival visualization for UK Biobank cohorts

## 0. RAP guardrails

This skill consumes a cohort `data.table` from `ukbsci-cohort` and produces
PNG / PDF figures plus an aggregate "data behind the figure" CSV. Figures
contain curves + numbers-at-risk only and may be shared with the local agent.
Per-participant intermediate tables (e.g. expanded long-format survival
sets) stay in `/mnt/project/...`.

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

---

## 1. When to load

- Draw a Kaplan-Meier curve, stratified or unstratified.
- Add a risk table beneath the KM curve.
- Show log-rank p-value comparing arms.
- Annotate median survival lines / censoring marks.

## 2. When NOT to load

- Fit / batch-fit Cox or Fine-Gray models → `ukbsci-regression`.
- Build the survival dataset or define cases → `ukbsci-cohort`.
- Subgroup × interaction tests → `ukbsci-subgroup-sensitivity`.
- Time-dependent ROC for ML predictions → `ukbsci-ml`.
- Forest / volcano / calibration plots → `ukbsci-plot`.

---

## 3. Prerequisites

```r
library(UKBAnalytica)
library(survival)
library(survminer)            # plot_km_curve uses ggsurvplot under the hood

stopifnot(all(c("outcome_status", "outcome_surv_time") %in% colnames(cohort)))
```

---

## 4. Pipeline

### Phase 1 — One-arm KM (no grouping)

```r
p <- plot_km_curve(
  data       = cohort,
  time_col   = "outcome_surv_time",
  status_col = "outcome_status",
  conf_int   = TRUE,
  risk_table = TRUE,
  title      = "Overall survival",
  xlab       = "Follow-up (years)",
  ylab       = "Event-free probability"
)
```

`plot_km_curve()` returns a `ggplot2` object (or a list `$km_plot + $risk_table`
when `risk_table = TRUE`). Save with the standard helper:

```r
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig03-km_overall.pdf",
                p, width = 7, height = 6)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig03-km_overall.png",
                p, width = 7, height = 6, dpi = 300)
```

### Phase 2 — Stratified KM

```r
p <- plot_km_curve(
  data         = cohort,
  time_col     = "outcome_surv_time",
  status_col   = "outcome_status",
  group_col    = "Hypertension_history",
  conf_int     = TRUE,
  risk_table   = TRUE,
  censor_marks = TRUE,
  palette      = "jco",              # "jco", "nejm", "lancet", or custom vector
  median_line  = TRUE,
  pvalue       = TRUE,                # log-rank p
  legend_title = "Baseline HTN",
  xlim         = c(0, 15),
  break_time   = 3
)
```

Log-rank p-value is computed **only when `group_col` is supplied**.

### Phase 3 — Save aggregate figure source data

`plot_km_curve()` does not write a source CSV automatically. Build it from
`survival::survfit()` and save alongside the figure for transparency:

```r
fit <- survfit(Surv(outcome_surv_time, outcome_status) ~ Hypertension_history,
               data = cohort)
km_data <- data.table::data.table(
  time     = fit$time,
  n.risk   = fit$n.risk,
  n.event  = fit$n.event,
  n.censor = fit$n.censor,
  surv     = fit$surv,
  lower    = fit$lower,
  upper    = fit$upper,
  strata   = rep(names(fit$strata), fit$strata)
)
fwrite(km_data,
       "/mnt/project/<area>/05-figs/data/Fig03-km_overall.csv")
```

---

## 5. Common pitfalls

1. **`plot_km_curve()` requires `survminer`** for the risk-table layout. If
   it is not installed on the RAP node, the user will see a missing-package
   error; install before running.
2. **`palette` accepts presets ("jco", "nejm", "lancet") or a custom vector
   of hex codes.** Mixing the two will silently fall back to a default.
3. **`xlim` units must match `time_col`.** `outcome_surv_time` in
   UKBAnalytica is in **years**; passing `xlim = c(0, 365)` would clip almost
   the whole curve.
4. **Prevalent participants.** They have `outcome_status = NA` and
   `outcome_surv_time = NA`; `survfit()` drops them automatically. Don't
   filter manually unless you intend to also drop censored cases.
5. **Group with `NA` levels.** `group_col` rows where the stratifier is `NA`
   are dropped silently. If the group sizes look smaller than the cohort,
   that is why.
6. **Risk-table fonts.** Default font sizes target a 7-inch-wide figure. For
   manuscript single-column (3.5 inch), shrink `fontsize = 2.5` in the
   ggsurvplot defaults (edit the returned object).
7. **Median-line for low event-rate cohorts.** If the curve never crosses
   0.5, the median line is invisible — this is correct behaviour, not a bug.

---

## 6. Key functions

| Function | Returns |
|----------|---------|
| `plot_km_curve(data, time_col, status_col, group_col = NULL, conf_int = TRUE, risk_table = TRUE, censor_marks = TRUE, palette = "jco", title = NULL, xlab = "Time (years)", ylab = "Survival Probability", legend_title = "Group", median_line = TRUE, pvalue = TRUE, xlim = NULL, break_time = NULL)` | `ggplot2` object (or list `$km_plot + $risk_table` when `risk_table = TRUE`) |

---

## 7. Related skills

| Skill | When to switch |
|-------|----------------|
| `ukbsci-cohort` | Re-build `outcome_status` / `outcome_surv_time` if missing. |
| `ukbsci-regression` | Fit Cox / Fine-Gray models on the same outcome. |
| `ukbsci-subgroup-sensitivity` | KM by interaction strata + interaction tests. |
| `ukbsci-plot` | All other figure families (forest, volcano, calibration). |
| `ukbsci-workflow` | Phase 6b — KM curves as the headline survival figure. |

---

## 8. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md) — overall KM, stratified
  KM, and "two-panel" KM (incidence + risk-table) example
