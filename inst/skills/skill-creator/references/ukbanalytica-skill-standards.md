# UKBAnalytica Skill Standards

## Privacy Boundary

All UKBAnalytica skills must make clear that a local AI agent is used for
script generation, workflow planning, package guidance, and interpretation of
aggregate outputs only. The agent must not read, inspect, summarize, or process
real UK Biobank participant-level RAP data.

The agent may:

- Generate scripts to run inside RAP.
- Review package code, documentation, function signatures, and simulated-data
  examples.
- Interpret aggregate outputs created inside RAP.

The agent must not receive, read, inspect, summarize, or process:

- UK Biobank RAP individual-level raw data.
- Direct identifiers such as `eid`.
- Exact participant dates.
- Raw RAP fields.
- Row-level source tables that can be linked back to participants.
- Real `.csv`, `.rds`, `.fst`, `.parquet`, `.feather`, or R objects containing
  participant rows, even when identifiers were removed.
- Console output, screenshots, tracebacks, or logs containing row-level values.

Acceptable material:

- De-identified rendered figures.
- Aggregate baseline tables and participant-flow counts.
- Regression coefficients, hazard ratios, odds ratios, CIs, and P values.
- Model metrics, ROC/PR/calibration/DCA curve data, and feature-level SHAP
  summaries; never raw per-row prediction tables or SHAP matrices.
- Enrichment, PPI, and feature-level omics results.

If uncertain, assume the object is participant-level and keep it inside RAP.

## Schema-Only Examples

If a skill needs runnable examples, use schema-only synthetic data. The user
may describe a dataframe by column names and roles, but the agent must not ask
for or inspect real UKB rows. Example code should include a small synthetic-data
smoke test and a separate RAP-executable script for the real analysis.

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
