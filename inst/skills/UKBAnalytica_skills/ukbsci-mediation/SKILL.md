---
name: ukbsci-mediation
description: >
  Causal mediation analysis for UK Biobank Research Analysis Platform (RAP)
  cohorts using UKBAnalytica. Wraps regmedint (Valeri & VanderWeele 4-way
  decomposition: CDE, PNDE, TNIE, TNDE, PNIE, TE, PM) via run_mediation;
  supports multi-mediator screening (run_multi_mediator),
  sensitivity-to-unmeasured-confounding (run_sensitivity_mediation), and
  visualization (plot_mediation, plot_mediation_forest) for linear,
  logistic, and Cox outcome models with continuous or binary mediators.
  Use this skill when the user asks for mediation analysis, indirect
  effects, natural direct / indirect effects, controlled direct effects,
  proportion mediated, or mediator screening across multiple candidates.
  Triggers: mediation analysis, indirect effect, natural direct effect,
  TNIE, PNDE, proportion mediated, 中介分析, /ukbsci-mediation. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-mediation — Causal mediation on UKB cohorts

## 0. RAP guardrails

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

Inputs: cohort `data.table` with exposure, mediator(s), outcome,
covariates. Outputs: aggregate effect tables + figures that may be shared
with the local agent.
Bootstrap CI sampling stays in RAP memory; the resulting estimates and
CIs are aggregate.

---

## 1. When to load

- The user has an exposure → mediator → outcome triple and wants to
  decompose total effect into direct + indirect components.
- Screening multiple candidate mediators for the strongest indirect effect.
- Reporting CDE / PNDE / TNIE / TNDE / PNIE / TE / proportion-mediated.
- Sensitivity assessment to unmeasured exposure–mediator confounders.

## 2. When NOT to load

- Simple confounder adjustment (no mediator) → `ukbsci-regression`.
- Effect-modification / subgroup tests → `ukbsci-subgroup-sensitivity`.
- IPTW / treatment-effect estimation → `ukbsci-propensity`.

---

## 3. Prerequisites

```r
library(UKBAnalytica)
library(regmedint)         # Required — implements 4-way decomposition
library(data.table)
```

If `regmedint` is missing, `run_mediation()` errors out with a clear
message.

---

## 4. Pipeline

### Phase 1 — Specify the triple + model types

```r
exposure   <- "sbp_above_median"          # binary or continuous
mediator   <- "ldl_cholesterol"           # continuous
outcome    <- "outcome_status"            # binary for logistic, ignored for Cox time arg
covars     <- c("age", "sex", "bmi", "smoking_status")

mediator_type <- "continuous"             # "continuous" or "binary"
outcome_type  <- "cox"                    # "linear", "logistic", "cox"
```

Decision matrix:

| Mediator | Outcome | `mediator_type` × `outcome_type` |
|----------|---------|----------------------------------|
| Continuous biomarker | Time-to-event | `("continuous", "cox")` |
| Binary indicator | Binary outcome | `("binary", "logistic")` |
| Continuous biomarker | Continuous outcome | `("continuous", "linear")` |

### Phase 2 — Single mediator: `run_mediation()`

```r
res <- run_mediation(
  data             = cohort,
  exposure         = exposure,
  mediator         = mediator,
  outcome          = "outcome_surv_time",            # time variable for Cox
  covariates       = covars,
  exposure_levels  = c(0, 1),
  mediator_value   = 0,                              # CDE evaluated at m = 0
  covariate_values = NULL,                           # NULL = covariate means
  mediator_type    = "continuous",
  outcome_type     = "cox",
  endpoint         = c("outcome_surv_time", "outcome_status"),
  interaction      = TRUE,                           # include E×M
  boot             = FALSE,
  boot_n           = 1000,
  conf_level       = 0.95
)

class(res)
# "mediation_result"

res$effects      # data.frame: estimate, SE, CI, p-value for CDE, PNDE, TNIE, TNDE, PNIE, TE, PM
summary(res)     # S3 summary with optional `exponentiate = TRUE`
coef(res)        # named vector of point estimates
confint(res)     # matrix [lower, upper]
```

### Phase 3 — Multi-mediator screening

```r
candidate_mediators <- c("ldl_cholesterol", "hba1c", "crp", "bmi")

multi <- run_multi_mediator(
  data           = cohort,
  exposure       = exposure,
  mediators      = candidate_mediators,
  outcome        = "outcome_surv_time",
  covariates     = covars,
  mediator_type  = "continuous",
  outcome_type   = "cox",
  endpoint       = c("outcome_surv_time", "outcome_status")
)
# data.frame: one row per mediator with TNIE/PNDE/TE/PM + SE/CI/p
```

### Phase 4 — Sensitivity to unmeasured confounding

```r
sens <- run_sensitivity_mediation(
  mediation_result = res,
  rho_values       = seq(-0.9, 0.9, by = 0.1)
)
# data.frame: rho, tnie_adjusted, pnde_adjusted, te_adjusted
```

Caveat: the implementation is a simplified ρ-based approximation, not a
fully rigorous bounds method. State this caveat to the user when reporting
the sensitivity analysis.

### Phase 5 — Visualize

```r
p_bar  <- plot_mediation(res, type = "effects",
                         show_ci = TRUE, show_pvalue = TRUE,
                         exponentiate = FALSE,
                         title = "Mediation of SBP→AA via LDL")
p_path <- plot_mediation(res, type = "path")
p_dec  <- plot_mediation(res, type = "decomposition")

p_forest <- plot_mediation_forest(
  multi,
  effect_type = "tnie",
  exponentiate = FALSE,
  null_value   = 0
)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig09-mediation_forest.pdf",
                p_forest, width = 8, height = 5)
```

---

## 5. Key assumptions encoded

`regmedint` (and therefore this skill) assumes:

1. No unmeasured exposure-outcome confounders given `covariates`.
2. No unmeasured mediator-outcome confounders given exposure + covariates.
3. No unmeasured exposure-mediator confounders given covariates.
4. No mediator-outcome confounder affected by exposure (cross-world
   assumption for natural effects).
5. Correct functional form of mediator and outcome models.

Setting `interaction = TRUE` allows for an E × M interaction in the
outcome model (recommended; the default).

---

## 6. Common pitfalls

1. **`regmedint` not installed.** Triggers a stop with explicit message —
   install before continuing.
2. **Survival outcomes** require `outcome_type = "cox"` *and* `endpoint`.
   Passing only `outcome` is for linear / logistic.
3. **Binary mediator + Cox outcome.** Supported but doubles model
   complexity; expect convergence warnings on small subgroups.
4. **Bootstrap is slow.** `boot = TRUE, boot_n = 1000` can take many
   minutes on 200k+ cohorts. Start with `boot = FALSE` (delta-method CIs)
   and only bootstrap a confirmatory run.
5. **No multiple-comparisons correction in `run_multi_mediator()`.**
   When screening ≥ 5 mediators, apply Benjamini-Hochberg externally before
   highlighting "significant" indirect effects.
6. **Sensitivity analysis is approximate.** State this in the manuscript;
   for rigorous bounds use the original VanderWeele / Imai approach.
7. **Exponentiated estimates** apply naturally to logistic / Cox outcomes
   (OR / HR scale). For linear outcomes, leave `exponentiate = FALSE`.

---

## 7. Key functions

| Function | Returns |
|----------|---------|
| `run_mediation(...)` | S3 `mediation_result` (`effects`, `mediator_model`, `outcome_model`, `regmedint_obj`, `params`) |
| `run_multi_mediator(...)` | data.frame: per-mediator effects |
| `run_sensitivity_mediation(mediation_result, rho_values)` | data.frame: ρ-adjusted effects |
| `plot_mediation(mediation_result, type = c("effects","path","decomposition"))` | ggplot |
| `plot_mediation_forest(multi_mediation_result, effect_type)` | ggplot |
| S3 methods | `print`, `summary`, `coef`, `confint` (on `mediation_result`) |

---

## 8. Related skills

| Skill | When |
|-------|------|
| `ukbsci-cohort` | Build cohort first. |
| `ukbsci-regression` | Compare total vs. mediator-adjusted estimates. |
| `ukbsci-imputation` | Pool mediation effects across imputations (use `pool_custom_estimates`). |
| `ukbsci-plot` | Polished forest figure for the manuscript. |

---

## 9. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md)
