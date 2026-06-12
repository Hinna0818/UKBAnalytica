---
name: ukbsci-workflow
description: >
  End-to-end orchestrator for a UK Biobank Research Analysis Platform (RAP)
  study using the UKBAnalytica R package. Plans, sequences, and supervises the
  full pipeline — from RAP phenotype extraction, through disease-cohort
  construction with reusable time skeletons, variable preprocessing, baseline
  tables, multiple imputation, regression / survival / propensity / mediation /
  subgroup-sensitivity analyses, machine learning, omics workflows, and final
  publication-ready plots — by routing
  each phase to the right specialist skill (ukbsci-rap-extract, ukbsci-cohort,
  ukbsci-preprocess, ukbsci-baseline, ukbsci-imputation, ukbsci-regression,
  ukbsci-survival, ukbsci-propensity, ukbsci-mediation,
  ukbsci-subgroup-sensitivity, ukbsci-proteomics, ukbsci-ml, ukbsci-plot).
  Use this skill when the user asks for an end-to-end UKB study, an analysis
  plan covering multiple phases, or "RAP-to-publication" guidance. Triggers:
  end-to-end UKB analysis, UKB pipeline, RAP to publication, full study plan,
  端到端 UKB 分析, RAP 到论文, 完整流程, 项目计划, /ukbsci-workflow.
  Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-workflow — UK Biobank end-to-end study orchestrator

## 0. Role of this skill

This skill is the **plan-and-route** layer. It does *not* do extraction,
cohort building, regression, ML, or plotting on its own — those are delegated
to the specialist `ukbsci-*` skills. Its responsibility is:

1. **Diagnose** the user's study scope from their first message.
2. **Produce a plan** (the document at `06-note/ukbsci-pipeline.md`) covering
   every phase, every output file, and every decision point.
3. **Hand off** each phase to the correct specialist skill in the right order.
4. **Gate** progress at known decision points so the user can pause / redirect.
5. **Write a wrap-up note** at the end (`06-note/ukbsci-README.md`) listing
   the actual deliverables.

Think of it as an analyst-PM: it owns the project plan, not the math.

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

---

## 1. When to load this skill

- The user asks for an end-to-end UKB study (`"help me run this analysis from
  RAP to publication"`).
- The user describes a multi-phase project and wants to know the right order.
- The user has already done some phases and wants you to plan the rest.
- The user is starting a new study and isn't sure which `ukbsci-*` skills are
  relevant.

## 2. When NOT to load this skill

- The user has a single, well-scoped question (`"how do I parse ICD-10
  codes?"`). Route directly to the relevant specialist skill.
- The user is debugging a specific function call. Route directly.
- The user wants a plot only. Route to `ukbsci-plot`.

---

## 3. Standard project layout

All outputs live inside RAP project storage. The canonical tree:

```
/mnt/project/<study-area>/
├── 01-script/
│   ├── 01-rap_extract.R
│   ├── 02-cohort_build.R
│   ├── 03-preprocess.R
│   ├── 04-baseline_table.R
│   ├── 05-imputation.R              # optional
│   ├── 06-regression.R
│   ├── 06b-survival.R
│   ├── 07-propensity.R              # optional
│   ├── 08-mediation.R               # optional
│   ├── 09-subgroup_sensitivity.R    # optional
│   ├── 10-ml.R                      # optional
│   └── 11-plots.R
├── 02-extract/                       # outputs from ukbsci-rap-extract
│   ├── env_check.rds
│   ├── manifest.csv
│   └── pheno.csv                    # stays on RAP
├── 03-cohort/                        # outputs from ukbsci-cohort
│   ├── cohort.csv
│   ├── time_skeleton.rds
│   ├── cohort_flow.csv              # aggregate — shareable with agent
│   └── source_compare.csv
├── 04-results/
│   ├── 01-baseline_table.csv
│   ├── 02-cox_main.csv
│   ├── 03-cox_subgroup.csv
│   ├── 04-cox_sensitivity.csv
│   ├── 05-propensity_balance.csv
│   ├── 06-mediation.csv
│   ├── 07-ml_metrics.csv
│   └── 08-ml_shap_summary.csv
├── 05-figs/
│   ├── Fig01-flow.pdf / .png
│   ├── Fig02-forest_main.pdf / .png
│   ├── Fig03-km.pdf / .png
│   ├── Fig04-volcano.pdf / .png
│   ├── Fig05-shap_summary.pdf / .png
│   └── data/                        # per-figure source tables
├── 06-note/
│   ├── ukbsci-pipeline.md           # Phase 0 plan
│   └── ukbsci-README.md             # Phase 8 wrap-up
└── 07-log/
    └── <one .log per script>
```

Naming rules the skill must enforce:

- Scripts are **numbered + descriptive**, one per phase.
- Figures are **dual-format** (`.pdf` + `.png` @ 300 dpi).
- Result tables sit in `04-results/` numbered to match the figure or table.
- Per-figure source data sits in `05-figs/data/` so reviewers can rebuild.
- Every script's `stdout`/`stderr` is captured to `07-log/<script>.log`.

---

## 4. The eight phases

### Phase 0 — Plan generation (this skill, always first)

**Before writing any analysis code**, produce `06-note/ukbsci-pipeline.md`.
The first generated script should also verify the active UKBAnalytica version
and key interfaces. If functions such as `build_survival_dataset()` or
`ukb_ml_flow()` are unavailable, install or load the latest project-approved
UKBAnalytica source before continuing.

For skill-maintenance checks, run the shared smoke test after interface
changes:

```r
Rscript inst/skills/UKBAnalytica_skills/evals/smoke_test_core.R
```

Required sections:

1. **Project metadata**: title, primary outcome, exposure, covariates,
   target population, censor date, expected sample size.
2. **Hypothesis & statistical question** — one sentence, framed as a
   testable claim.
3. **Phase-by-phase tree diagram** (ASCII), naming every script, its inputs,
   outputs, and the user-visible decisions inside it.
4. **Decision points** (numbered) — places where the agent *must* pause and
   ask the user.
5. **Deliverables** — final tables and figures, with intended downstream use
   (manuscript table, supplementary figure, etc.).

Sample tree (use as a template — adapt to the user's study):

```
01-rap_extract.R           [skill: ukbsci-rap-extract]
  ├─ Input:   field IDs from §1.exposure, §1.covariates, outcome sources
  ├─ Output:  /mnt/project/<area>/02-extract/pheno.csv (RAP-resident)
  └─ Decision: sync vs async; whether to include rare fields

02-cohort_build.R          [skill: ukbsci-cohort]
  ├─ Input:   02-extract/pheno.csv
  ├─ Output:  03-cohort/cohort.csv, 03-cohort/cohort_flow.csv
  └─ Decision: primary_disease; prevalent_sources; outcome_sources

03-preprocess.R            [skill: ukbsci-preprocess]
  ├─ Input:   03-cohort/cohort.csv
  ├─ Output:  03-cohort/cohort_clean.csv
  └─ Decision: which negative-code conversions; derived variables

04-baseline_table.R        [skill: ukbsci-baseline]
  ├─ Input:   03-cohort/cohort_clean.csv
  ├─ Output:  04-results/01-baseline_table.csv
  └─ Decision: strata column; which factor / continuous variables

05-imputation.R (optional) [skill: ukbsci-imputation]
  ├─ Input:   03-cohort/cohort_clean.csv
  ├─ Output:  03-cohort/imp_list.rds
  └─ Decision: m; method per variable; mi-pool aggregation

06-regression.R            [skill: ukbsci-regression]
  ├─ Input:   03-cohort/cohort_clean.csv (or imp_list.rds)
  ├─ Output:  04-results/02-cox_main.csv
  └─ Decision: model family; covariates; transformations

06b-survival.R             [skill: ukbsci-survival]
  ├─ Input:   03-cohort/cohort_clean.csv
  ├─ Output:  04-results/02b-km_data.csv, 05-figs/data/Fig03-km.csv
  └─ Decision: stratifier; competing risks vs cause-specific

07-propensity.R (optional) [skill: ukbsci-propensity]
08-mediation.R (optional)  [skill: ukbsci-mediation]
09-subgroup_sensitivity.R  [skill: ukbsci-subgroup-sensitivity]
10-ml.R (optional)         [skill: ukbsci-ml]
11-plots.R                 [skill: ukbsci-plot]
  ├─ Input:   04-results/*.csv
  ├─ Output:  05-figs/Fig01-..-FigNN-*.{pdf,png}
  └─ Decision: palette family (ukbsci_clinical / ukbsci_diverging / ukbsci_sequential)
```

**Always present the plan to the user verbatim** and wait for one of two
responses:

- `"开始分析"` / `"start"` → proceed to Phase 1.
- `"调整计划"` / `"revise"` → ask which sections to change.

Do not start Phase 1 without explicit approval.

### Phase 1 — RAP extraction

Hand off to `ukbsci-rap-extract`. Pass:

- The field IDs / variable sets derived from the plan.
- The chosen mode (`sync` vs `job`) based on `plan$n_columns`.
- The exact output path (`/mnt/project/<area>/02-extract/pheno.csv`).
- A RAP environment check and field manifest:
  `ukb_check_rap_env()` + `ukb_create_extraction_manifest()`.

Verify before continuing:

```r
stopifnot(file.exists("/mnt/project/<area>/02-extract/pheno.csv"))
```

### Phase 2 — Cohort build

Hand off to `ukbsci-cohort`. Confirm the four decision points
(primary disease, prevalent/outcome source sets, censor date) with the user
**before** the script runs. Prefer `ukb_time_skeleton()` first, then
`build_survival_dataset(time_skeleton = ...)`. After it runs, present the
aggregate attrition flow; for custom inclusion/exclusion rules use
`ukb_participant_flow()`.

### Phase 3 — Preprocessing

Hand off to `ukbsci-preprocess`. Standard steps: negative-code sanitation,
unit harmonisation, composite variables (BP, diet score, air pollution).

### Phase 4 — Baseline (Table 1)

Hand off to `ukbsci-baseline`. Confirm strata variable with the user.
Output goes to `04-results/01-baseline_table.csv` as an aggregate table that
may be shared with the local agent.

### Phase 5 — Imputation (optional)

Only run if missingness in covariates is non-trivial. Hand off to
`ukbsci-imputation`. Output is an imputed list, not a single table — every
downstream regression must pool across imputations.

### Phase 6 — Statistical modeling

Hand off to one or more of:

- `ukbsci-regression` (Cox, logistic, linear, batch + extensions)
- `ukbsci-survival` (KM, competing risks, time-dependent ROC)
- `ukbsci-propensity` (PSM, IPTW, balance plots)
- `ukbsci-mediation` (single / multi-mediator, sensitivity)
- `ukbsci-subgroup-sensitivity` (interaction tests, exclusion criteria)

The exact subset depends on the study design. Confirm with the user.

Before running expensive models, generate a small simulated-data smoke test or
use package demo data to validate function signatures. Never expose real
participant rows to the local agent for debugging.

### Phase 7 — Machine learning (optional)

Hand off to `ukbsci-ml`. Common configuration: outcome from cohort, features
from preprocessed table, models = `c("xgboost", "glmnet", "ranger")`. Standard
deliverables: ROC/PR/calibration/DCA + SHAP summary.

### Phase 8 — Plots & wrap-up

Hand off to `ukbsci-plot`. Read source tables from `04-results/`, write
figures to `05-figs/`, source data per figure to `05-figs/data/`.

Then **write the wrap-up note** at `06-note/ukbsci-README.md` containing:

1. Citation header for `UKBAnalytica`.
2. ASCII tree of *actual* files produced (use `list.files()` — do not paste a
   template).
3. Key findings (3–5 bullets) summarising effect sizes, n's, and any
   sensitivity-vs-main discrepancies.
4. Caveats: imputation status, source selection, censoring date, any sources
   omitted.

Present a one-paragraph summary to the user with figure count, log file
status, and the top results.

---

## 5. Decision points the skill must enforce

| # | Phase | Pause and ask |
|---|-------|---------------|
| D0 | Phase 0 | Approve / revise the plan |
| D1 | Phase 1 | Sync vs async extraction |
| D2 | Phase 2 | `primary_disease`, `prevalent_sources`, `outcome_sources` |
| D3 | Phase 3 | Negative-code list, composite variable formulas |
| D4 | Phase 4 | Strata variable, which covariates appear in Table 1 |
| D5 | Phase 5 | Run imputation? `m`, methods |
| D6 | Phase 6 | Which models? Covariate set? Sensitivity variants? |
| D7 | Phase 7 | Run ML? Which models and feature set? |
| D8 | Phase 8 | Palette family; figure aspect ratios |

If the user has already pre-specified a decision, skip the prompt — but echo
back the choice ("Got it — using `primary_disease = AA` per your earlier
message").

---

## 6. Citation header (every emitted script)

```r
###############################################################################
# UKBAnalytica Citation:
# He N. UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for
# UK Biobank RAP Data. R package version 1.0.0.
# https://github.com/Hinna0818/UKBAnalytica
###############################################################################
```

---

## 7. Log capture (every script)

```bash
mkdir -p /mnt/project/<area>/07-log

Rscript 01-script/01-rap_extract.R \
   > /mnt/project/<area>/07-log/01-rap_extract.log 2>&1

Rscript 01-script/02-cohort_build.R \
   > /mnt/project/<area>/07-log/02-cohort_build.log 2>&1

# ...continue for every script
```

The wrap-up note must record `wc -l` for each log so the user can verify
non-empty runs.

---

## 8. Common pitfalls

1. **Skipping Phase 0.** Do not start coding before the plan is approved.
   Even when the user is impatient — the plan unblocks every subsequent
   decision.
2. **Phase-skipping.** Don't run regression before preprocessing produces the
   clean cohort, even if "it's just a quick check".
3. **Inconsistent paths.** The plan fixes `/mnt/project/<area>/...`; every
   sub-skill must honor it. If a sub-skill tries to write to `~/`, intercept.
4. **Forgetting prevalent NAs in regression.** `outcome_status = NA` for
   prevalent — that is *intentional*. Sub-skills must not impute or zero
   these out.
5. **Pooled-imputation mismatch.** When `ukbsci-imputation` is used, every
   downstream regression must use `fit_mi_models()` + `pool_mi_models()` —
   not a single-imputation point estimate.
6. **Source mismatches between plan and run.** If during Phase 6 the user
   changes `outcome_sources`, re-run Phase 2 instead of editing
   `cohort.csv` in place.
7. **Plot styling drift.** All figures should share one palette family from
   `ukbsci-plot`. Mixing within a single manuscript looks unprofessional.

---

## 9. Related skills

| Phase | Skill |
|-------|-------|
| 1 | `ukbsci-rap-extract` |
| 2 | `ukbsci-cohort` |
| 3 | `ukbsci-preprocess` (P5) |
| 4 | `ukbsci-baseline` (P3) |
| 5 | `ukbsci-imputation` (P4) |
| 6 | `ukbsci-regression` (P3), `ukbsci-survival` (P3), `ukbsci-propensity` (P4), `ukbsci-mediation` (P4), `ukbsci-subgroup-sensitivity` (P4) |
| 7 | `ukbsci-ml` (P5), `ukbsci-proteomics` (P5) |
| 8 | `ukbsci-plot` (P6) |

(Roadmap phase tags in parentheses indicate when the specialist skill will
ship in the full pack — see `supp/UKBAnalytica-skill-roadmap.md`.)

---

## 10. References

- [`references/plan-template.md`](references/plan-template.md) — fillable
  `ukbsci-pipeline.md` skeleton
- [`references/wrap-up-template.md`](references/wrap-up-template.md) —
  `ukbsci-README.md` skeleton for Phase 8
- [`references/handoff-checklist.md`](references/handoff-checklist.md) — what
  the workflow skill must hand to each specialist skill
