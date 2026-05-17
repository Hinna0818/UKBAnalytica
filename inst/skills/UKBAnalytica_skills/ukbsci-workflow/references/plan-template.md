# `ukbsci-pipeline.md` — Phase 0 plan template

Fill this in for every new study. Save the populated file to
`/mnt/project/<study-area>/06-note/ukbsci-pipeline.md`. **Wait for user
approval (`"开始分析"` / `"start"`) before Phase 1.**

---

## 1. Project metadata

- **Title**:
- **PI / analyst**:
- **Date drafted** (UTC):
- **UKB application ID**:
- **RAP project**:
- **`UKBAnalytica` version**: `<packageVersion("UKBAnalytica")>`
- **Censor date**: `as.Date("YYYY-MM-DD")` (must match the user's UKB release)
- **Expected analytic sample size**: ~

## 2. Hypothesis

> One sentence, falsifiable. Example: *"Among UK Biobank participants free
> of aortic aneurysm at baseline, higher baseline systolic blood pressure
> is independently associated with incident AA after adjustment for age, sex,
> BMI, smoking status, and diabetes history."*

## 3. PICO (or its observational analogue)

- **P (Population)**:
- **E (Exposure)**:
- **C (Comparator)**:
- **O (Outcome)**:
- **T (Time horizon)**:

## 4. Variables

### 4.1 Outcome

| Disease key | Definition source | Sources for incidence |
|-------------|-------------------|------------------------|
| `<primary>` | `get_predefined_diseases()` or custom | `c("ICD10","ICD9","Death")` |

### 4.2 Exposure

| Variable | UKB field ID(s) | Type | Transformation |
|----------|------------------|------|-----------------|

### 4.3 Covariates

| Variable | UKB field ID(s) | Type | Rationale |
|----------|------------------|------|-----------|

### 4.4 Stratifiers / subgroups (optional)

| Variable | Levels | Why |
|----------|--------|-----|

### 4.5 Sensitivity flips

| Variant | Spec | Rationale |
|---------|------|-----------|
| Strict hospital | `outcome_sources = c("ICD10","ICD9")` | drop death only |
| Algorithmic | `outcome_sources = "Algorithm"` | UKB-validated dates |
| Exclude early events | drop events < 1y from baseline | reverse-causation |

## 5. Phase tree (8 phases)

```
01-rap_extract.R                [ukbsci-rap-extract]
  ├─ Input:    field IDs + variable sets above
  ├─ Output:   /mnt/project/<area>/02-extract/pheno.csv
  ├─ Mode:     sync | async  (decide from rap_plan_extract()$n_columns)
  └─ Decision: which dictionary version; whether array fields needed

02-cohort_build.R               [ukbsci-cohort]
  ├─ Input:    02-extract/pheno.csv
  ├─ Output:   03-cohort/cohort.csv + 03-cohort/cohort_flow.csv
  └─ Decision: primary_disease; prevalent_sources; outcome_sources

03-preprocess.R                 [ukbsci-preprocess]
  ├─ Input:    03-cohort/cohort.csv
  ├─ Output:   03-cohort/cohort_clean.csv
  └─ Decision: negative-code handling; composite variables

04-baseline_table.R             [ukbsci-baseline]
  ├─ Input:    03-cohort/cohort_clean.csv
  ├─ Output:   04-results/01-baseline_table.csv
  └─ Decision: strata column; factor / continuous lists

05-imputation.R   (optional)    [ukbsci-imputation]
  ├─ Input:    03-cohort/cohort_clean.csv
  ├─ Output:   03-cohort/imp_list.rds
  └─ Decision: m, methods per variable

06-regression.R                 [ukbsci-regression]
  ├─ Input:    cohort_clean.csv  OR  imp_list.rds
  ├─ Output:   04-results/02-cox_main.csv
  └─ Decision: model family; transforms; mi-pool yes/no

06b-survival.R                  [ukbsci-survival]
  ├─ Input:    cohort_clean.csv
  ├─ Output:   04-results/02b-km_data.csv
  └─ Decision: stratifier; competing risks vs cause-specific

07-propensity.R   (optional)    [ukbsci-propensity]
08-mediation.R    (optional)    [ukbsci-mediation]
09-subgroup_sensitivity.R       [ukbsci-subgroup-sensitivity]
10-ml.R           (optional)    [ukbsci-ml]
11-plots.R                      [ukbsci-plot]
  ├─ Input:    04-results/*.csv
  ├─ Output:   05-figs/Fig*.pdf|.png  + 05-figs/data/*.csv
  └─ Decision: palette family (ukbsci_clinical / ukbsci_diverging / ukbsci_sequential)
```

## 6. Decision points

| # | Phase | Question |
|---|-------|----------|
| D0 | Phase 0 | Approve / revise the plan? |
| D1 | Phase 1 | Sync vs async? (decide from `plan$n_columns`) |
| D2 | Phase 2 | `primary_disease`, `prevalent_sources`, `outcome_sources` |
| D3 | Phase 3 | Negative-code list, composite variable formulas |
| D4 | Phase 4 | Strata column; which vars in Table 1 |
| D5 | Phase 5 | Run imputation? `m`, methods |
| D6 | Phase 6 | Which models? Covariates? Sensitivity variants? |
| D7 | Phase 7 | Run ML? Which models and feature set? |
| D8 | Phase 8 | Palette family; figure aspect ratios |

## 7. Deliverables

| Deliverable | Path | Intended use |
|-------------|------|--------------|
| Table 1 | `04-results/01-baseline_table.csv` | Manuscript Table 1 |
| Main Cox table | `04-results/02-cox_main.csv` | Manuscript Table 2 / Fig 2 source |
| Subgroup Cox table | `04-results/03-cox_subgroup.csv` | Manuscript Fig 3 source |
| Sensitivity Cox | `04-results/04-cox_sensitivity.csv` | Supplementary |
| KM curves | `05-figs/Fig03-km.{pdf,png}` | Manuscript Fig 1 |
| Forest plot (main) | `05-figs/Fig02-forest_main.{pdf,png}` | Manuscript Fig 2 |
| Subgroup forest | `05-figs/Fig04-forest_subgroup.{pdf,png}` | Manuscript Fig 3 |
| Volcano (if multi-exposure) | `05-figs/Fig05-volcano.{pdf,png}` | Manuscript Fig 4 |
| ML metrics | `04-results/07-ml_metrics.csv` | Supplementary |
| SHAP summary | `05-figs/Fig06-shap_summary.{pdf,png}` | Supplementary |

## 8. RAP guardrails reminder

The local agent must not read or inspect real participant-level RAP files or
R objects. It should generate scripts for RAP execution and interpret only
aggregate outputs returned by the user. All participant-level files stay under
`/mnt/project/<area>/...`. Only the following aggregate outputs or
non-identifying rendered figures can be shared with the agent:

- `cohort_flow.csv` (aggregate counts)
- `01-baseline_table.csv` (aggregate summary)
- `02-cox_main.csv` and every other `04-results/*.csv` (coefficients + CIs)
- `05-figs/*.{pdf,png}` (rendered figures)
- `05-figs/data/*.csv` (figure source — verify it contains aggregates only,
  not per-participant rows)
- `06-note/*.md` (text)
- `07-log/*.log` only after confirming they contain no row-level values,
  identifiers, exact dates, or data previews

Anything in `02-extract/` or `03-cohort/` that retains per-participant rows
**must not leave RAP**.

---

> **Ready?** Reply `开始分析` to proceed to Phase 1, or `调整计划` and tell me
> which sections to change.
