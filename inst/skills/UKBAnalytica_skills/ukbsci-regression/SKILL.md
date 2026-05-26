---
name: ukbsci-regression
description: >
  Fit linear, logistic, Cox, GLM, negative-binomial, GAM, Fine-Gray competing-risk,
  lagged Cox, trend, PH-diagnostic, and restricted cubic spline models
  on UK Biobank Research Analysis Platform (RAP) cohorts built with UKBAnalytica.
  Wraps runmulti_lm, runmulti_logit, runmulti_cox, runmulti_glm,
  runmulti_negbin, runmulti_gam, run_regression, runmulti_cox_lag,
  runmulti_cox_zph, runmulti_competing, runmulti_trend, ukb_cox_diagnostics,
  run_rcs, and plot_rcs. Use this skill when the user asks to fit any regression model
  on the cohort, run batch regressions across multiple exposures, fit count/rate
  outcomes with Poisson/quasi-Poisson/Gamma GLM, fit negative-binomial for
  overdispersed counts, screen for non-linearity with GAM, check proportional-hazards
  assumptions, run a lagged sensitivity analysis, perform competing-risks Fine-Gray
  modeling, test a dose-response trend, or draw an exposure-response RCS curve.
  Triggers: UKB regression, Cox model,
  logistic regression, GLM, Poisson, negative binomial, GAM, batch regression,
  proportional hazards, competing risks, Fine-Gray, p_trend, restricted cubic spline,
  RCS, 非线性, UKB 回归, Cox 模型, 广义线性模型, 广义加性模型, /ukbsci-regression.
  Hard rule: local agents must not
  read or inspect real UKB RAP participant-level data; generate scripts for RAP
  execution and interpret aggregate outputs only.
---

# ukbsci-regression — UK Biobank batch regression

## 0. RAP guardrails

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

Input: a cohort `data.table` produced by `ukbsci-cohort` (or any
preprocessed UKB analysis dataset) that **already lives in RAP project
storage**. Output: aggregate coefficient tables (estimate / CI / p-value),
which may be shared with the local agent.

Forbidden: writing per-participant predicted values / residuals / influence
diagnostics to local disk. Keep those in `/mnt/project/...` if you need them
for QC.

---

## 1. When to load

- Fit one or several univariate / multivariable Cox, logistic, linear, GLM,
  negative-binomial, or GAM models on a UKB cohort.
- Batch-run regression across **many exposures at once** (each gets its own model).
- Count or rate outcomes with Poisson, quasi-Poisson, Gamma, or other GLM families.
- Overdispersed count outcomes (negative-binomial as a principled alternative to
  quasi-Poisson).
- Screen for non-linear dose-response associations with GAM splines.
- Fit and visualize restricted cubic spline exposure-response curves.
- Check the proportional-hazards assumption for Cox models.
- Run a lagged sensitivity analysis (drop events in the first N years).
- Fine-Gray competing-risks modeling.
- Linear-trend ("p-trend") testing across ordered exposure categories.

## 2. When NOT to load

- Building the cohort or survival dataset → `ukbsci-cohort`.
- Subgroup × interaction tests → `ukbsci-subgroup-sensitivity`.
- IPTW / propensity-adjusted regression → `ukbsci-propensity`.
- Multiple-imputation pooling → `ukbsci-imputation`.
- KM curves / cumulative incidence / time-dependent ROC → `ukbsci-survival`.
- Forest plots of the output → `ukbsci-plot`.

---

## 3. Prerequisites

```r
library(UKBAnalytica)
library(data.table)
library(survival)         # Cox

stopifnot(data.table::is.data.table(cohort))
# For Cox: outcome columns from build_survival_dataset()
stopifnot(all(c("outcome_status", "outcome_surv_time") %in% colnames(cohort)))
```

For competing-risks models the user needs a separate competing-event
indicator (e.g. non-cardiovascular death dates), or a single status column
encoded as 0 = censored, 1 = primary event, 2 = competing event.

---

## 4. Pipeline (4 phases)

### Phase 1 — Pick the right model family

| Question | Function | Outcome |
|----------|----------|---------|
| Continuous outcome | `runmulti_lm()` | numeric (e.g. SBP, BMI, biomarker) |
| Binary outcome at fixed time | `runmulti_logit()` | 0/1 |
| Time-to-event with censoring | `runmulti_cox()` | `outcome_status` + `outcome_surv_time` |
| Count / rate outcome (Poisson, Gamma…) | `runmulti_glm()` | any GLM family |
| Overdispersed count outcome | `runmulti_negbin()` | non-negative integer |
| Non-linear dose-response screening | `runmulti_gam()` | any smooth GAM family |
| Time-to-event + competing risk | `runmulti_competing()` | status with `{0,1,2}` coding |
| Cox with PH diagnostics | `runmulti_cox_zph()` | same as Cox |
| Cox with early-event exclusion | `runmulti_cox_lag()` | same as Cox |
| Ordered exposure trend | `runmulti_trend()` | depends on `model_type` |
| Restricted cubic spline | `run_rcs()` + `plot_rcs()` | Cox, logistic, or linear |

**Unified interface:** `run_regression(type = c("lm","logit","cox","glm","negbin","gam"))` wraps
all six base model families — use it when the model type is determined programmatically or you
want a single consistent call site.

### Phase 2 — Choose exposures & covariates

The `*_multi_*` family runs **one model per `main_var`** but shares the same
covariate set. Build the covariate set from the analysis plan, e.g.

```r
covars <- c("age", "sex", "bmi", "smoking_status",
            "Hypertension_history", "Diabetes_history")
exposures <- c("sbp", "ldl_cholesterol", "hba1c", "crp")
```

### Phase 3 — Fit

```r
# Linear: continuous outcome, one model per exposure
lm_res <- runmulti_lm(
  data       = cohort,
  main_var   = exposures,
  covariates = covars,
  outcome    = "carotid_imt"
)

# Logistic: binary outcome at baseline
logit_res <- runmulti_logit(
  data       = cohort,
  main_var   = exposures,
  covariates = covars,
  outcome    = "Hypertension_history"
)

# Cox: time-to-event
cox_res <- runmulti_cox(
  data       = cohort,
  main_var   = exposures,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status")
)

# GLM: Poisson for count / rate outcome (IRR = exp(beta))
poisson_res <- runmulti_glm(
  data       = cohort,
  main_var   = exposures,
  family     = "poisson",
  outcome    = "hospitalisation_count",
  covariates = covars
)
# exp(poisson_res$beta) gives the incidence rate ratio
# quasi-Poisson uses Wald CI automatically: family = "quasipoisson"
# Gamma with log link for positive-continuous biomarkers:
#   family = stats::Gamma(link = "log")

# Negative-binomial: overdispersed count outcome
negbin_res <- runmulti_negbin(
  data       = cohort,
  main_var   = exposures,
  outcome    = "hospitalisation_count",
  covariates = covars
)
# Output: IRR (= exp(beta)), lower95, upper95, pvalue, theta (dispersion), n

# GAM: non-linearity screening
gam_res <- runmulti_gam(
  data       = cohort,
  main_var   = exposures,
  outcome    = "carotid_imt",
  covariates = covars
)
# edf ≈ 1 → near-linear; edf > 1 → non-linear association
# smooth = FALSE for a parametric linear term inside GAM framework

# Save aggregate results
fwrite(cox_res,     "/mnt/project/<area>/04-results/02-cox_main.csv")
fwrite(negbin_res,  "/mnt/project/<area>/04-results/02-negbin_main.csv")
```

#### Restricted cubic spline exposure-response

Use RCS when the user asks for an interpretable nonlinear exposure-response
curve, especially for continuous environmental, biomarker, or omics exposures.
For quick screening across many exposures, `runmulti_gam()` is lighter; for a
publication-ready single-exposure curve, use `run_rcs()` and `plot_rcs()`.

```r
rcs_fit <- run_rcs(
  data       = cohort,
  exposure   = "pm25",
  covariates = covars,
  model_type = "cox",
  endpoint   = c("outcome_surv_time", "outcome_status"),
  knots      = 4,
  ref_value  = median(cohort$pm25, na.rm = TRUE),
  backend    = "ns"
)

rcs_plot <- plot_rcs(
  rcs_fit,
  xlab = "PM2.5",
  ylab = "HR (95% CI)",
  show_p = TRUE,
  distribution = "rug"
)
ggplot2::ggsave(
  "/mnt/project/<area>/04-results/02-rcs_pm25.pdf",
  rcs_plot, width = 5.2, height = 4.2
)
```

Interpretation: `p_overall` tests whether the exposure is associated with the
outcome; `p_nonlinear` tests whether the spline improves over a linear term.
Do not export prediction grids with participant IDs; aggregate curve data and
figures are safe to share.

#### Unified interface

```r
# run_regression() — single call site for all 6 families
cox_res2 <- run_regression(
  data       = cohort,
  main_var   = exposures,
  type       = "cox",
  endpoint   = c("outcome_surv_time", "outcome_status"),
  covariates = covars
)

run_regression(data = cohort, main_var = exposures,
               type = "glm", family = "poisson",
               outcome = "hospitalisation_count", covariates = covars)

run_regression(data = cohort, main_var = exposures,
               type = "gam", outcome = "carotid_imt",
               smooth = TRUE, covariates = covars)
```

### Phase 4 — Diagnostics & extensions

#### 4a. PH diagnostics (Schoenfeld residual test)

```r
zph <- runmulti_cox_zph(
  data       = cohort,
  main_var   = exposures,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status"),
  transform  = "km",
  alpha      = 0.05
)
zph[, .(variable, HR, lower95, upper95, pvalue,
        zph_pvalue, ph_violation)]
```

A `TRUE` in `ph_violation` means the exposure may have a time-varying effect;
consider stratifying the offending covariate or adding a `tt()` term in a
follow-up model.

#### 4b. Lagged sensitivity

Drops participants whose event occurred within `lag_years` of baseline —
guards against reverse causation.

```r
lag_res <- runmulti_cox_lag(
  data      = cohort,
  main_var  = exposures,
  covariates = covars,
  endpoint   = c("outcome_surv_time", "outcome_status"),
  lag_years  = c(0, 1, 2, 5),
  verbose    = TRUE
)
fwrite(lag_res, "/mnt/project/<area>/04-results/04-cox_lag.csv")
```

The output is *long*: one row per `(lag_years, variable, level)` triple.

#### 4c. Competing risks (Fine-Gray subdistribution hazard)

Two coding modes:

```r
# Single-column mode: status_col has {0 = censored, 1 = primary, 2 = competing}
fg1 <- runmulti_competing(
  data           = cohort,
  main_var       = exposures,
  covariates     = covars,
  time_col       = "outcome_surv_time",
  event_col      = "outcome_status_full",
  event_value    = 1,
  compete_value  = 2
)

# Dual-column mode: separate columns for primary and competing events
fg2 <- runmulti_competing(
  data        = cohort,
  main_var    = exposures,
  covariates  = covars,
  time_col    = "outcome_surv_time",
  event_col   = "primary_event",
  compete_col = "competing_event"
)
```

Output: subdistribution hazard ratios (`SHR`) instead of `HR`.

#### 4d. Trend across ordered exposure categories

```r
trend_res <- runmulti_trend(
  data            = cohort,
  main_var        = c("ldl_cholesterol_q", "sbp_q"),  # quartile factors
  covariates      = covars,
  model_type      = "cox",
  endpoint        = c("outcome_surv_time", "outcome_status"),
  ref_level       = "Q1",
  score_method    = "integer",                        # 1, 2, 3, 4
  include_level_estimates = TRUE
)
trend_res[, .(variable, level, estimate, lower95, upper95,
              pvalue, p_trend, trend_estimate)]
```

Output: per-level estimates **plus** an overall `p_trend` per exposure.

#### 4e. Bespoke PH diagnostics for a single fitted Cox object

```r
fit <- survival::coxph(
  survival::Surv(outcome_surv_time, outcome_status) ~ sbp + age + sex,
  data = cohort
)
diag <- ukb_cox_diagnostics(fit, transform = "km", alpha = 0.05)
diag$table
diag$global_pvalue
```

---

## 5. Common pitfalls

1. **Prevalent NA in `outcome_status`.** `build_survival_dataset()` codes prevalent cases as
   `NA` so Cox drops them. Do **not** impute these to 0.
2. **Factor reference level.** All regression functions honor the factor's first level as
   reference. Use `relevel()` (or set `ref_level` in `runmulti_trend`) before fitting.
3. **Multicollinearity.** Batch regressions don't auto-check VIF. Confirm `covars` doesn't
   contain correlated variables (e.g. `sbp` and `hypertension`).
4. **Wrong endpoint order.** `endpoint = c("time", "status")` — time first.
5. **GLM coefficient scale.** `runmulti_glm()` returns coefficients on the **link scale**.
   For log-link families, exponentiate: `exp(beta)`, `exp(lower95)`, `exp(upper95)`.
6. **quasi-Poisson vs negative-binomial.** quasi-Poisson adjusts SEs only; negative-binomial
   is a proper likelihood model. Prefer `runmulti_negbin()` when reporting AIC or under heavy
   overdispersion.
7. **GAM non-linearity interpretation.** `edf` close to 1 does not prove linearity — it
   reflects the penalised fit. Visual inspection of the smooth is also recommended.
8. **RCS overfitting.** Use sensible knot counts (usually 3–5), prespecify or
   report the reference value, and avoid interpreting unstable tails with sparse data.
9. **Competing-risks status coding.** Value `1` = primary event, `2` = competing event.
   Confirm with the user; mis-coding inverts the result.
10. **`runmulti_cox_lag()` drops events with `surv_time ≤ lag_years`.** Does not shift
   the time origin.
11. **`runmulti_trend()` requires factor exposures.** Quartile-categorise first if continuous.

---

## 6. Key functions (cheat-sheet)

| Function | Returns | Notable columns |
|----------|---------|-----------------|
| `runmulti_lm(data, main_var, covariates, outcome)` | data.frame | `variable, beta, lower95, upper95, pvalue` |
| `runmulti_logit(data, main_var, covariates, outcome)` | data.frame | `variable, OR, lower95, upper95, pvalue` |
| `runmulti_cox(data, main_var, covariates, endpoint)` | data.frame | `variable, HR, lower95, upper95, pvalue, n, n_event` |
| `runmulti_glm(data, main_var, family, outcome, covariates)` | data.frame | `variable, family, link, beta, lower95, upper95, pvalue, n` |
| `runmulti_negbin(data, main_var, outcome, covariates)` | data.frame | `variable, IRR, lower95, upper95, pvalue, theta, n` |
| `runmulti_gam(data, main_var, outcome, covariates, smooth, family, k)` | data.frame | smooth=TRUE: `edf, ref_df, stat, stat_type, pvalue, family, link, n`; smooth=FALSE: `beta, lower95, upper95` |
| `run_regression(data, main_var, type, outcome, endpoint, family, smooth, covariates)` | data.frame | as per underlying function |
| `runmulti_cox_zph(..., transform, alpha)` | data.frame | + `zph_pvalue, ph_violation, ph_global_violation` |
| `runmulti_cox_lag(..., lag_years)` | data.frame (long) | + `lag_years, n_removed` |
| `runmulti_competing(data, main_var, time_col, event_col, compete_col, event_value, compete_value)` | data.frame | `SHR, lower95, upper95, pvalue, n_compete` |
| `runmulti_trend(..., model_type, score_method, custom_scores)` | data.frame | per-level + `p_trend, trend_estimate` |
| `ukb_cox_diagnostics(model, transform, alpha)` | list | `table, global_pvalue, cox_zph` |
| `run_rcs(data, exposure, covariates, model_type, endpoint, outcome, knots, backend)` | `ukb_rcs` list | `prediction, p_overall, p_nonlinear, knots, ref` |
| `plot_rcs(x, ...)` | ggplot | RCS curve with 95% CI and optional exposure distribution |

**`runmulti_glm()` family notes:**
- Accepts character (`"poisson"`, `"quasipoisson"`, `"gaussian"`, `"binomial"`, `"Gamma"`…),
  a function (`stats::poisson`), or a family object (`stats::Gamma(link = "log")`).
- Quasi-families use Wald CIs; proper families use profile-likelihood CIs.
- Log/logit/cloglog link: coefficients are on the link scale; exponentiate for ratio estimates.

**`runmulti_negbin()` vs quasi-Poisson:**
- Negative-binomial is a proper probability model; enables AIC comparison, more reliable under
  heavy overdispersion.
- Quasi-Poisson only adjusts standard errors; coefficients are identical to Poisson.

**`runmulti_gam()` interpretation:**
- `edf` (estimated degrees of freedom): ≈ 1 → near-linear; > 1 → non-linear.
- A common workflow: screen with `smooth = TRUE`, then switch to `runmulti_lm()` for linear
  variables (`edf < 1.5`) and retain GAM smooth for non-linear ones.

See [`references/functions.md`](references/functions.md) for full argument
tables.

---

## 7. Related skills

| Skill | When to switch |
|-------|----------------|
| `ukbsci-cohort` | Build / re-build the survival dataset before regressing. |
| `ukbsci-survival` | KM curves, time-dependent ROC, post-Cox visualization. |
| `ukbsci-propensity` | IPTW-weighted regression; PS matching analysis. |
| `ukbsci-subgroup-sensitivity` | Subgroup × interaction tests. |
| `ukbsci-mediation` | Mediation modeling. |
| `ukbsci-imputation` | Pool regression coefficients across imputations. |
| `ukbsci-plot` | Forest plot, volcano plot of regression output. |
| `ukbsci-workflow` | End-to-end orchestrator; calls this skill in Phase 6. |

---

## 8. References

- [`references/functions.md`](references/functions.md) — full signatures
- [`references/examples.md`](references/examples.md) — three end-to-end
  examples (Cox main, Cox with lag, Fine-Gray competing risks)
