# UK Biobank Research Scoping Report Template

Use this template when `$ukb-research` returns a full Markdown report.
Adapt section length to the user's needs, but preserve the core logic.

```markdown
# Research Scoping Report: <Topic>

## 1. Executive Summary

Briefly state the proposed research question, why it matters, whether it is
feasible in UK Biobank, and the recommended primary design.

## 2. User Question and Assumptions

- User-provided keywords:
- Proposed exposure or feature domain:
- Proposed outcome:
- Study population:
- Assumptions made during this scoping assessment:

## 3. Literature Background

Summarize the relevant evidence in 2-4 concise paragraphs. Distinguish:

- what has already been studied in UK Biobank or similar cohorts;
- what is biologically or clinically plausible;
- what remains underexplored;
- whether prior papers provide reusable code definitions or analytical designs.

All factual claims about previous studies must be cited with PMID, DOI, or URL.

## 4. Existing Work and Limitations

### 4.1 What has already been done

Summarize the closest prior studies. Focus on design features that shape the
new research idea: cohort, exposure, outcome, data modality, statistical model,
validation strategy, and main conclusion.

### 4.2 Limitations of existing studies

Use citation-supported bullets. Common limitation categories include:

- phenotype ascertainment or endpoint heterogeneity;
- cross-sectional design or weak temporality;
- limited exposure windows or single-time-point measurements;
- insufficient control of confounding or mediation pathways;
- lack of nonlinear, subgroup, or sensitivity analyses;
- limited omics integration or biological interpretation;
- lack of internal/external validation;
- poor reproducibility or missing codelists.

### 4.3 Actionable research gap

State one primary gap that can be addressed with UK Biobank and available
software. Avoid vague claims such as "few studies have examined this" unless
the literature table supports it.

### 4.4 Proposed gap-resolving strategy

Explain how the proposed design addresses the gap. Link each design choice to
a limitation above.

## 5. Evidence Table

| Study | Data source | Exposure/features | Outcome | Design/method | Key finding | Reusable elements |
|---|---|---|---|---|---|---|
| Author/year, PMID/DOI | UKB / other |  |  |  |  | phenotype codes, covariates, model, software |

## 6. Feasibility in UK Biobank

| Domain | Assessment | Rationale |
|---|---|---|
| Phenotype ascertainment | High/Moderate/Low | ICD/self-report/death/cancer/first occurrence availability |
| Exposure availability | High/Moderate/Low | field IDs, repeat measures, assay platform |
| Temporality | High/Moderate/Low | baseline exposure before incident outcome |
| Event count | High/Moderate/Low | expected prevalence/incidence or literature estimate |
| Covariates | High/Moderate/Low | demographics, lifestyle, socioeconomic, clinical variables |
| Missingness risk | High/Moderate/Low | likely missing variables and handling plan |
| Novelty | High/Moderate/Low | comparison with prior studies |
| Overall feasibility | High/Moderate/Low | final judgement |

## 7. Recommended Study Design

### 7.1 Primary objective

State one primary objective in a testable form.

### 7.2 Study population

- Source population:
- Inclusion criteria:
- Exclusion criteria:
- Time zero:
- Follow-up end:

### 7.3 Exposure definition

Specify variables, coding, scaling, transformations, categories, and whether
the exposure is baseline, longitudinal, omics-derived, or externally linked.

### 7.4 Outcome definition

Specify disease definitions, code systems, and sources. Prefer:

- `get_disease_catalog()` to inspect code evidence;
- `get_predefined_diseases()` for built-in endpoints;
- `create_disease_definition()` only when the endpoint is custom.

### 7.5 Covariates

List core and optional covariates with rationale.

### 7.6 Primary analysis

Describe the statistical model, effect measure, covariate adjustment, and
multiple-testing correction if relevant.

### 7.7 Secondary analyses

Examples: stratified analysis, interaction testing, nonlinear modeling,
protein/metabolite enrichment, prediction modeling, mediation, or calibration.

### 7.8 Sensitivity analyses

Examples: alternative endpoint sources, complete-case vs imputation, lagged
event exclusion, additional covariate adjustment, competing-risk model,
negative controls, or alternative exposure windows.

### 7.9 How this design addresses the gap

Map the proposed analyses to the limitations identified in Section 4.

| Prior limitation | Proposed design response | Expected output |
|---|---|---|
| Example: weak temporality | prospective incident-outcome design | HR estimates and lag sensitivity |

## 8. Expected Results and Interpretation

Describe plausible result patterns, including positive, null, and inconsistent
findings. Avoid overstating causality unless the design supports it.

## 9. Implementation Plan

| Step | Task | Recommended software | UKBAnalytica functions | Output |
|---|---|---|---|---|
| 1 | Field planning | R / UKBAnalytica | `rap_plan_extract()`, `ukb_create_extraction_manifest()` | extraction manifest |
| 2 | Baseline preprocessing | R / UKBAnalytica | `preprocess_baseline()`, `get_variable_info()` | clean covariates |
| 3 | Phenotype definition | R / UKBAnalytica | `get_disease_catalog()`, `get_predefined_diseases()`, `build_survival_dataset()` | survival-ready cohort |
| 4 | Descriptive analysis | R / UKBAnalytica | `create_baseline_table()` | Table 1 |
| 5 | Main modeling | R / UKBAnalytica | `runmulti_cox()` / `runmulti_logit()` / `run_rcs()` | effect estimates |
| 6 | Advanced analysis | R/Python/external tools | project-specific | enrichment, ML, MR, GWAS, figures |

## 10. Figures and Tables

Propose manuscript-ready outputs:

- participant flow diagram;
- baseline characteristics table;
- main association forest plot or volcano plot;
- nonlinear RCS curve if needed;
- subgroup/sensitivity summary;
- omics enrichment/PPI plot or ML ROC/SHAP plots when relevant.

## 11. Limitations and Risk Mitigation

List design-specific limitations and how to address them.

## 12. Reference Verification

Report whether references were verified.

```bash
python3 inst/skills/ukb-research/scripts/verify_references.py <report>.md \
  --output reference_audit.json
```

Summary:

- PMIDs checked:
- DOIs checked:
- Unverified references:
- Notes:

## 13. Recommended Next Steps

1. Confirm phenotype and exposure availability.
2. Generate a RAP extraction manifest.
3. Run a small simulated-data script locally or an aggregate-only pilot on RAP.
4. Finalize statistical analysis plan before full execution.

## References

Use PMID, DOI, and URLs where available. Do not include references that cannot
be traced to a real source unless explicitly marked as unverified.
```
