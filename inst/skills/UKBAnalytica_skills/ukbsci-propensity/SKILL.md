---
name: ukbsci-propensity
description: >
  Propensity-score analyses for UK Biobank Research Analysis Platform (RAP)
  cohorts built with UKBAnalytica. Covers PS estimation via logistic
  regression or gradient-boosted models (estimate_propensity_score), 1:N
  nearest-neighbor or optimal matching (match_propensity), inverse-probability
  of treatment weighting (calculate_weights, ATE / ATT / ATC, stabilized,
  truncated), covariate balance diagnostics (assess_balance, plot_balance,
  plot_ps_distribution), and weighted regression with robust standard errors
  (run_weighted_analysis). Use this skill when the user asks for propensity
  score matching, PSM, IPTW, ATE / ATT, covariate balance, Love plot, PS
  distribution check, or a weighted Cox / logistic / linear analysis on a
  UKB cohort. Triggers: propensity score, PSM, IPTW, ATE, ATT, Love plot,
  covariate balance, 倾向评分, 倾向得分匹配, /ukbsci-propensity. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-propensity — Propensity scores, matching, and IPTW on UKB cohorts

## 0. RAP guardrails

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

The matched / weighted dataset is **participant-level** and stays in RAP
project storage. Shareable artefacts are aggregate: balance tables, Love
plots, PS distribution figures, and final weighted-effect estimates.

---

## 1. When to load

- Estimate a propensity score for a binary treatment / exposure.
- Build a 1:N PSM-matched sub-cohort.
- Construct IPTW weights (ATE / ATT / ATC, stabilized, optionally truncated).
- Assess covariate balance before vs after matching / weighting.
- Run a weighted Cox / logistic / linear analysis with robust SEs.

## 2. When NOT to load

- Vanilla unadjusted / multivariable Cox → `ukbsci-regression`.
- Multiple-imputation pooling of weighted estimates → use this skill **plus**
  `ukbsci-imputation`.
- Mediation analysis → `ukbsci-mediation`.

---

## 3. Prerequisites

```r
library(UKBAnalytica)
library(data.table)

# Suggested packages
# - MatchIt for matching
# - gbm for gradient-boosted PS
# - sandwich, lmtest for HC0 robust SEs
# - survival for weighted Cox

stopifnot("treatment" %in% colnames(cohort))
stopifnot(all(cohort$treatment %in% c(0, 1, NA)))
```

---

## 4. Pipeline (5 phases)

### Phase 1 — Estimate PS

```r
covars_ps <- c("age", "sex", "bmi", "smoking_status",
               "Hypertension_history", "Diabetes_history",
               "ldl_cholesterol", "sbp")

ps_data <- estimate_propensity_score(
  data       = cohort,
  treatment  = "treatment",
  covariates = covars_ps,
  method     = "logistic"             # "logistic" or "gbm"
)
# Adds column `ps`, bounded to [0.001, 0.999]
```

For non-linear / non-additive relationships, prefer `method = "gbm"`.

### Phase 2 — Choose matching OR weighting

**Matching path:**

```r
matched <- match_propensity(
  data        = ps_data,
  ps_col      = "ps",
  treatment   = "treatment",
  ratio       = 1,
  caliper     = 0.2,
  method      = "nearest",            # or "optimal"
  replace     = FALSE,
  exact_match = c("sex")              # optional exact-match strata
)
# Adds `match_id`, `match_distance`
```

**Weighting path:**

```r
weighted <- calculate_weights(
  data         = ps_data,
  ps_col       = "ps",
  treatment    = "treatment",
  weight_type  = "ATE",               # "ATE", "ATT", or "ATC"
  stabilized   = TRUE,
  truncate     = c(0.01, 0.99)        # weight quantile clipping
)
# Adds `weight`
```

### Phase 3 — Assess balance (before / after)

```r
b_pre  <- assess_balance(
  data       = cohort,
  treatment  = "treatment",
  covariates = covars_ps,
  method     = "unmatched"
)

b_post_m <- assess_balance(
  data       = matched,
  treatment  = "treatment",
  covariates = covars_ps,
  method     = "matched"
)

b_post_w <- assess_balance(
  data       = weighted,
  treatment  = "treatment",
  covariates = covars_ps,
  method     = "weighted",
  weight_col = "weight",
  threshold  = 0.1
)

p_love <- plot_balance(
  balance_before = b_pre,
  balance_after  = b_post_m,             # or b_post_w for IPTW
  threshold      = 0.1,
  title          = "Love plot — before vs. after matching"
)

p_ps <- plot_ps_distribution(
  data      = ps_data,
  ps_col    = "ps",
  treatment = "treatment",
  type      = "mirror"               # "histogram", "density", "mirror"
)
```

The standard target: **SMD ≤ 0.1** on every covariate post-adjustment.

### Phase 4 — Weighted / matched analysis

```r
# Cox under matched cohort (no weighting; matching already balanced)
cox_matched <- runmulti_cox(
  data       = matched,
  main_var   = "treatment",
  covariates = NULL,
  endpoint   = c("outcome_surv_time", "outcome_status")
)

# Or IPTW-weighted Cox with robust SE
wt_cox <- run_weighted_analysis(
  data       = weighted,
  exposure   = "treatment",
  covariates = NULL,
  weight_col = "weight",
  model_type = "cox",
  endpoint   = c("outcome_surv_time", "outcome_status"),
  robust_se  = TRUE
)
```

`run_weighted_analysis()` also supports `model_type = "logistic"` /
`"linear"`. Robust standard errors via `sandwich` + `lmtest`.

### Phase 5 — Export

All output goes to `/mnt/project/<area>/04-results/`:

```r
fwrite(b_pre,    "04-results/05a-balance_unmatched.csv")
fwrite(b_post_m, "04-results/05b-balance_matched.csv")
fwrite(b_post_w, "04-results/05c-balance_weighted.csv")
fwrite(wt_cox,   "04-results/02d-cox_iptw.csv")

ggplot2::ggsave("05-figs/Fig07-love_plot.pdf", p_love, width = 7, height = 6)
ggplot2::ggsave("05-figs/Fig08-ps_distribution.pdf", p_ps, width = 7, height = 5)
```

---

## 5. Common pitfalls

1. **Bounded PS at `[0.001, 0.999]`.** Prevents extreme weights but
   shrinks the tails. Verify with `plot_ps_distribution()` that no group
   has piled-up boundary mass.
2. **Caliper too tight.** Default `0.2 × SD(PS)` is generous; tightening
   below `0.05 × SD` can drop most of the treated arm.
3. **Replace = TRUE.** Increases effective sample but inflates variance;
   warn the user and confirm.
4. **`weight_type = "ATE"` vs `"ATT"`** changes the estimand. PSM by default
   targets ATT (the matched cohort). Mixing in IPTW with `ATE` and reporting
   them side-by-side without flagging is misleading.
5. **Truncation defaults at 1st/99th percentile.** Aggressive — if the user
   wants to keep extremes, set `truncate = c(0, 1)`.
6. **SMD on factors.** `assess_balance()` converts factors to numeric for
   SMD; multi-level factors are summarized as a single SMD. For per-level
   SMD, expand the factor into dummies first.
7. **Robust SE for Cox.** `run_weighted_analysis(model_type = "cox",
   robust_se = TRUE)` uses `survival::coxph(..., robust = TRUE)` semantics —
   appropriate for IPTW. For matched Cox, prefer **cluster-robust SE on
   `match_id`** (build the formula manually).
8. **Missing `weight` column for `method = "weighted"`** triggers a clear
   error from `assess_balance()`. Pass `weight_col = "weight"` explicitly.

---

## 6. Key functions

| Function | Signature highlights | Returns |
|----------|---------------------|---------|
| `estimate_propensity_score(data, treatment, covariates, method = c("logistic","gbm"), formula = NULL)` | adds `ps` ∈ [0.001, 0.999] | data.table |
| `match_propensity(data, ps_col, treatment, ratio, caliper, method, replace, exact_match)` | wraps MatchIt::matchit | data.table + `match_id`, `match_distance` |
| `calculate_weights(data, ps_col, treatment, weight_type, stabilized, truncate)` | IPTW weights | data.table + `weight` |
| `assess_balance(data, treatment, covariates, method, weight_col, threshold)` | SMD + variance ratio | data.frame with `balanced` flag |
| `plot_balance(balance_before, balance_after, threshold, title, xlab)` | Love plot | ggplot |
| `plot_ps_distribution(data, ps_col, treatment, type, matched, match_col)` | PS overlap viz | ggplot |
| `run_weighted_analysis(data, exposure, outcome, covariates, weight_col, model_type, endpoint, robust_se)` | IPTW model | data.frame: HR / OR / β + 95% CI |

See [`references/functions.md`](references/functions.md) for full argument tables.

---

## 7. Related skills

| Skill | When |
|-------|------|
| `ukbsci-cohort` | Build the cohort + `outcome_status` first. |
| `ukbsci-baseline` | Post-PSM Table 1 with SMD column. |
| `ukbsci-regression` | Unweighted multivariable adjustment for comparison. |
| `ukbsci-plot` | Forest plot of weighted results. |

---

## 8. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md) — three examples
  (PSM with caliper + exact match, IPTW with stabilized ATE weights,
  weighted Cox with robust SE)
