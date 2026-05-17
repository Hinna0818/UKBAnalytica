# UKBAnalytica Skill Standards

## Privacy Boundary

All UKBAnalytica skills must distinguish prohibited raw-data export from
acceptable analytical outputs.

Do not export:

- UK Biobank RAP individual-level raw data.
- Direct identifiers such as `eid`.
- Exact participant dates.
- Raw RAP fields.
- Row-level source tables that can be linked back to participants.

Generally exportable:

- De-identified rendered figures.
- Aggregate baseline tables and participant-flow counts.
- Regression coefficients, hazard ratios, odds ratios, CIs, and P values.
- Model metrics, ROC/PR/calibration/DCA curve data, and feature-level SHAP
  summaries.
- Enrichment, PPI, and feature-level omics results.

## Package-First Examples

Use UKBAnalytica functions in examples whenever possible. Avoid re-implementing
package logic inside skills.

Good examples:

- `rap_plan_extract()` before extraction.
- `preprocess_baseline()` for standardized covariates.
- `get_predefined_diseases()` and `build_survival_dataset()` for endpoints.
- `runmulti_cox()` for protein-outcome Cox screens.
- `ukb_ml_flow()` or `ukb_ml_compare_flows()` for ML workflows.
- `plot_regression_volcano()`, `plot_ml_roc_compare()`, and
  `plot_shap_beeswarm()` for figures.

## Version Reminder

When a workflow uses recently added functions, instruct users to install or
load the latest UKBAnalytica source before running examples.

## Skill Naming

- Use `ukbsci-*` for UKB scientific analysis skills.
- Use direct task names for meta-skills, e.g. `skill-creator`.
- Keep names lowercase and hyphenated.
