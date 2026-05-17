# Hand-off checklist

What the `ukbsci-workflow` skill should pass to each specialist skill at the
start of its phase. Every hand-off is a structured prompt: study context +
inputs + decisions already locked + outputs expected.

When the workflow skill activates a specialist, it should send a single
message of the form:

```
Phase <n> hand-off from ukbsci-workflow.

Study context:
  - Title: <...>
  - Hypothesis: <...>
  - Working directory: /mnt/project/<area>/
Inputs:
  - <path-1>
  - <path-2>
Locked decisions:
  - <decision-key>: <value>
  - <decision-key>: <value>
Open decisions for you:
  - <question-1>
  - <question-2>
Expected outputs:
  - <path-1>
  - <path-2>
```

Below are the per-phase hand-off contracts.

---

## Phase 1 → `ukbsci-rap-extract`

**Inputs from workflow:**

- The full field-ID list (numeric IDs + any predefined variable sets).
- The chosen `dataset` (if not auto-discovered).
- Output path: `/mnt/project/<area>/02-extract/pheno.csv`.

**Locked decisions:**

- Whether to use sync (`rap_extract_pheno`) or async (`rap_submit_extract`).
- `manifest` path: `/mnt/project/<area>/02-extract/manifest.csv`.

**Open decisions for the specialist:**

- The exact instance / array column set after `rap_plan_extract()` expansion.

**Expected outputs:**

- `02-extract/pheno.csv` (RAP-resident)
- `02-extract/manifest.csv` (field-level only; shareable)

**Post-condition the workflow must check:**

```r
stopifnot(file.exists("/mnt/project/<area>/02-extract/pheno.csv"))
```

---

## Phase 2 → `ukbsci-cohort`

**Inputs:**

- `02-extract/pheno.csv`.
- Disease list from the plan.
- Defaults for `prevalent_sources` and `outcome_sources`.

**Locked decisions:**

- `primary_disease`.
- `censor_date`.
- `baseline_col = "p53_i0"`.
- `output = "wide"` (unless user overrides).

**Open decisions for the specialist:**

- Whether to add `OPCS4` / `CancerRegistry` / `FirstOccurrence` / `Algorithm`
  to either source set.

**Expected outputs:**

- `03-cohort/cohort.csv`
- `03-cohort/cohort_flow.csv` (aggregate; shareable)
- `03-cohort/source_compare.csv` (optional aggregate; shareable)

**Post-condition:**

- `attr(cohort, "participant_flow")` must be non-NULL and presented to the
  user.

---

## Phase 3 → `ukbsci-preprocess`

**Inputs:**

- `03-cohort/cohort.csv`.
- Variable specifications from the plan (which columns are exposure,
  covariates, etc.).

**Locked decisions:**

- Negative-code sanitation (`-1`, `-2`, `-3`, `-7`, `-10` → `NA`).

**Expected outputs:**

- `03-cohort/cohort_clean.csv`.

---

## Phase 4 → `ukbsci-baseline`

**Inputs:**

- `03-cohort/cohort_clean.csv`.
- Strata column (from D4).
- Variable lists for `factor_vars` and continuous vars.

**Expected outputs:**

- `04-results/01-baseline_table.csv` (aggregate; shareable).

---

## Phase 5 → `ukbsci-imputation` (optional)

**Inputs:**

- `03-cohort/cohort_clean.csv`.
- Imputed-variable list and `m` (from D5).

**Expected outputs:**

- `03-cohort/imp_list.rds` (RAP-resident).

---

## Phase 6 → `ukbsci-regression` / `ukbsci-survival` / `ukbsci-propensity` / `ukbsci-mediation` / `ukbsci-subgroup-sensitivity`

The workflow skill may invoke several of these in sequence. Each gets:

**Inputs:**

- `03-cohort/cohort_clean.csv` (or `imp_list.rds` when Phase 5 ran).
- Outcome columns: `outcome_status`, `outcome_surv_time`.
- Exposure column(s).
- Covariate list.

**Locked decisions:**

- Model family per skill (`runmulti_cox` / `runmulti_logit` / `runmulti_lm`,
  etc.).
- Whether to pool across imputations (`fit_mi_models` + `pool_mi_models`)
  when Phase 5 ran.

**Open decisions:**

- Subgroup variables (subgroup-sensitivity).
- Mediator(s) and exposure-mediator interaction (mediation).
- Caliper / replace / target estimand (propensity).

**Expected outputs (per skill):**

| Skill | Output |
|-------|--------|
| `ukbsci-regression` | `04-results/02-cox_main.csv` (and variants) |
| `ukbsci-survival` | `04-results/02b-km_data.csv` |
| `ukbsci-propensity` | `04-results/05-propensity_balance.csv` |
| `ukbsci-mediation` | `04-results/06-mediation.csv` |
| `ukbsci-subgroup-sensitivity` | `04-results/03-cox_subgroup.csv`, `04-results/04-cox_sensitivity.csv` |

All outputs are aggregate (estimates, CIs, p-values) and may be shared with
the local agent.

---

## Phase 7 → `ukbsci-ml` (optional)

**Inputs:**

- `03-cohort/cohort_clean.csv`.
- Feature list and outcome column.

**Locked decisions:**

- Task (`binary`, `multiclass`, `survival`).
- Train/test split seed (record in `06-note/ukbsci-pipeline.md`).

**Expected outputs:**

- `04-results/07-ml_metrics.csv`.
- `04-results/08-ml_shap_summary.csv`.

---

## Phase 8 → `ukbsci-plot`

**Inputs:**

- Every CSV in `04-results/`.
- Style decision (D8): palette family.

**Expected outputs:**

- `05-figs/Fig01-..-FigNN-*.pdf` + matching `.png`.
- `05-figs/data/Fig*.csv` source tables.

After the plotting phase, the workflow skill writes
`06-note/ukbsci-README.md` from the wrap-up template.

---

## Aggregate-only sharing manifest (final step)

The workflow skill produces a manifest enumerating aggregate files and
non-identifying rendered figures that may be shared with the local agent.
The user must review the manifest before sharing and remove any file that
contains row-level values, identifiers, exact dates, data previews, or raw RAP
fields. Generate it with:

```r
exportable <- c(
  list.files("/mnt/project/<area>/04-results", full.names = TRUE),
  list.files("/mnt/project/<area>/05-figs",   full.names = TRUE,
             pattern = "\\.(png|pdf)$"),
  # Include figure data CSVs only after confirming they are aggregate
  # (per-bin, per-stratum, per-feature, or per-cluster), never per-participant.
  list.files("/mnt/project/<area>/05-figs/data", full.names = TRUE,
             pattern = "\\.csv$"),
  list.files("/mnt/project/<area>/06-note",   full.names = TRUE),
  "/mnt/project/<area>/03-cohort/cohort_flow.csv",
  "/mnt/project/<area>/02-extract/manifest.csv"
)
writeLines(exportable, "/mnt/project/<area>/06-note/shareable_aggregate_outputs.txt")
```

The workflow skill then tells the user:

> Everything in `shareable_aggregate_outputs.txt` should be aggregate-level or
> a rendered non-identifying figure before it is shared with the local agent.
> Everything else (extract / cohort
> participant-level tables) must remain in the RAP project and out of the
> agent context.
