# UKBAnalytica News

## UKBAnalytica 1.1.0 (2026-08-01)

### RAP phenotyping and genetic analysis

- Added Spark-first primary-care GP querying, clinical and registration parsing, coverage assessment, and coverage-aware participant-level disease integration.
- Added aligned GWAS/PheWAS phenotype preparation and auditable wrappers for REGENIE and PLINK2, including raw official argument-token escape hatches.
- Added genotype conversion planning across BED, PGEN, BGEN, VCF, and BCF formats.
- Added published UK Biobank PRS catalogue/loading interfaces and custom PLINK2 PRS calculation from normalized, build-checked, allele-harmonized weights.
- Added multi-PRS wide scoring, chromosome-specific weight splitting, BGEN/PGEN/BED inputs, controlled parallel execution, and resumable chromosome runs.

### Documentation

- Added complete GWAS/PheWAS, primary-care GP, and PRS workflow chapters to the Quarto website.
- Expanded README examples for RAP extraction, disease phenotyping, genetic analyses, and reproducible PRS calculation.

## UKBAnalytica 1.0.0 (2026-05-16)

### Official release
- Released UKBAnalytica 1.0.0 as an integrated RAP-oriented workflow package for UK Biobank epidemiological analyses.
- Expanded the package description to reflect the current end-to-end scope, including RAP extraction planning, predefined variables and diseases, cohort construction, statistical modeling, machine learning, proteomics analysis, and publication-oriented visualization.

### Workflow coverage
- Consolidated predefined disease definitions and baseline variable mappings for rapid phenotype and covariate setup.
- Standardized the analysis flow from RAP data extraction and preprocessing to multi-source endpoint definition, prevalent/incident classification, survival-ready cohort generation, and downstream analyses.
- Added high-level machine-learning interfaces for single-model workflows and comparison across feature sets or algorithms, with ROC, threshold selection, calibration, SHAP interpretation, and validation-set evaluation helpers.
- Extended proteomics utilities for volcano plots, Gene Ontology enrichment, STRING-based PPI retrieval, topology metrics, fast-greedy community detection, and cluster-level enrichment analysis.
- Refined visualization helpers for forest plots, volcano plots, calibration and decision-curve analysis, ROC comparison, SHAP summaries, enrichment plots, and manuscript-ready figure export.

### Recent updates
- Added optional supplementary nested cross-validation through `ukb_ml_nested_cv()` and `ukb_ml_workflow(nested_cv = TRUE)`, with fold-level performance and feature-selection stability while preserving the default frozen-test workflow.
- Added `ukb_time_skeleton()` to create a reusable follow-up time skeleton with baseline date, death date, loss-to-follow-up date, administrative censoring, follow-up end reason, and valid follow-up indicators.
- Added optional `time_skeleton` support to `build_survival_dataset()` while preserving the default survival workflow.
- Added RAP-protected dictionary utilities, including `ukb_download_rap_dictionary()`, `ukb_query_dictionary()`, and `ukb_validate_columns()`, for official RAP dictionary lookup, Chinese/English field search, and column validation.
- Added the built-in `ukb_dictionary_zh` metadata dataset for Chinese UKB field-path lookup.
- Extended `run_regression()` with `covariate_sets` for nested epidemiological models such as crude, partially adjusted, and fully adjusted analyses.
- Added lightweight publication-style plotting helpers for heatmaps, stacked bars, violin plots, and scatter plots.
- Excluded local `tests/` from package builds and remote tracking.

## UKBAnalytica 0.6.2.2（2026-04-29）

### Disease phenotyping
- Added UK Biobank cancer registry support as a new `CancerRegistry` disease source using fields 40006, 40005, 40011, and 40012.
- Added `cancer_icd10_pattern`, `cancer_histology`, and `cancer_behaviour` to `create_disease_definition()`.
- Added predefined `Lung_Cancer` with cancer registry, ICD-10, and death-registry ascertainment.
- Added UK Biobank First Occurrence support as a new `FirstOccurrence` disease source.
- Added `first_occurrence_fields` and `first_occurrence_source_fields` to `create_disease_definition()`.
- Added internal parsing for singular `p13xxxx` First Occurrence date/source fields, including UKB special date coding 819 handling.
- Added common First Occurrence field IDs to predefined diseases where the disease definition matches 3-character ICD-10 field granularity.

### Workflow helpers
- Added `ukb_clean_missing()` for converting common UKB non-response labels and numeric missing codes into analysis-ready values.
- Added `ukb_snapshot()` to record row/column counts, missingness, complete rows, object size, and deltas across analysis pipeline checkpoints.

### Machine learning
- Added the high-level `ukb_ml_workflow()` API for binary, multiclass, and continuous non-survival ML with a frozen final test set.
- Added standardized ML split, feature selection, tuning, threshold-learning, final-refit, and frozen-test evaluation helpers: `ukb_ml_as_split()`, enhanced `ukb_ml_split_data()`, `ukb_ml_feature_select()`, `ukb_ml_tune()`, `ukb_ml_threshold()`, `ukb_ml_fit_final()`, and `ukb_ml_evaluate_test()`.
- Preserved the legacy `split_ratio` style in `ukb_ml_split_data()` by keeping `$internal_validation` as an alias for the held-out split.
- Extended `ukb_shap()` to support `ukb_ml_workflow` and `ukb_ml_final` objects, defaulting to the frozen test set for workflow objects.
- Added lightweight `rpart` decision tree and `naive_bayes` model backends to `ukb_ml_workflow()`.
- Added optional lazy ML dependency installation via `options(UKBAnalytica.auto_install_ml = TRUE)`; by default, optional model packages are checked only when the selected model needs them and are not installed automatically.
- Marked legacy ML model, CV, prediction, importance, comparison, and evaluation helpers as deprecated with pointers to the new workflow APIs.
- Added `ukb_ml_survival_workflow()` and survival-specific split, feature-selection, tuning, final-refit, and frozen-test evaluation helpers for time-to-event ML.
- Added `model = "cox"` as the lightweight default survival ML backend and aligned survival prediction output with the new workflow object structure.
- Marked legacy `ukb_ml_survival()` as deprecated in favor of `ukb_ml_survival_workflow()`.
- Added simulated-data tests for binary classification, multiclass classification, continuous regression, manual split validation, and threshold behavior.
- Added simulated-data tests for survival ML splitting, manual split validation, Cox workflow fitting, prediction, and frozen-test evaluation.

### RAP phenotype extraction
- Added R-native RAP phenotype extraction helpers that wrap DNAnexus `dx extract_dataset` and RAP `table-exporter`.
- Added `rap_find_dataset()`, `rap_list_fields()`, `rap_plan_extract()`, `rap_extract_pheno()`, and `rap_submit_extract()`.
- Added dry-run extraction planning so users can inspect matched fields and table-exporter command metadata before launching jobs.
- Added support for `variables = ...` using UKBAnalytica predefined baseline mappings, while preserving `field_id = ...` for all instances and arrays of a UKB field.
- Kept Python extraction helper scripts under `inst/python/` as legacy/helper entry points.

### Documentation and tests
- Expanded the RAP extraction chapter with a complete RAP JupyterLab -> Terminal -> R workflow, including output paths such as `/mnt/project`.
- Updated README examples to recommend R-native RAP extraction for new users.
- Split machine learning into a dedicated documentation chapter focused on the simplest standard `ukb_ml_workflow()` path.
- Added offline tests for RAP extraction field parsing, field planning, predefined-variable mapping, and table-exporter dry-run formatting.

## UKBAnalytica 0.6.2.1（2026-04-24）

### Disease definition and phenotyping updates
- Added `OPCS4` operative procedure support for hospital summary operations via `p41272` + `p41282_a*`.
- Added `opcs4_pattern` to `create_disease_definition()` so procedure evidence is opt-in and ignored by default when unspecified.
- Extended case extraction and survival workflows to accept `OPCS4` in `sources`, `prevalent_sources`, and `outcome_sources`.
- Added predefined arrhythmia endpoints combining ICD-10 and OPCS4 where appropriate: `Arrhythmia`, `Ventricular_Arrhythmia`, `AV_Block`, `Intraventricular_Block`, and `SVT`.
- Extended predefined `Atrial_Fibrillation` with OPCS4 support for procedure-augmented atrial arrhythmia ascertainment.

### Documentation
- Updated disease definition chapter and main vignette examples to document `opcs4_pattern` and arrhythmia phenotyping with `ICD10 + OPCS4`.
- Updated `README.md` with an ICD-10 + OPCS4 phenotyping example and clarified the default opt-in behavior for procedure data.


## UKBAnalytica 0.6.2 (2026-04-18)

### Survival workflow and participant flow reporting
- Enhanced `build_survival_dataset()` with `show_flow` to print step-by-step participant attrition in terminal for wide output.
- Added flow summary fields for each step (`n_before`, `n_after`, `excluded`, retention rates from previous/raw cohort).
- Attached the underlying flow table to returned wide-format result via `attr(result, "participant_flow")`.
- Added optional `dt_threads` in `build_survival_dataset()` to let users temporarily configure `data.table` thread count for large runs.

### Robust date parsing and stability improvements
- Added internal `.safe_as_date()` utility (`R/date_utils.R`) to parse mixed date formats safely and convert malformed values to `NA` with warnings instead of stopping execution.
- Replaced direct `as.Date()` calls in key pipelines with `.safe_as_date()` (ICD, death, baseline, incident-time utilities, and case extraction paths).
- Fixed self-report date conversion in `parse_self_reported_illnesses()` to handle malformed year values (`Inf`, `-Inf`, `NaN`, non-numeric strings) without `charToDate` crashes.

### Algorithm source compatibility
- Updated algorithm field lookup to support both `p{field}_i0` and `p{field}` naming conventions for date/source fields.
- Improved diagnostic messages when algorithm columns are unavailable.
- Updated disease definition documentation for algorithm field naming compatibility.

### Diabetes diagnosis method refinement
- Refined diabetes phenotyping workflow for prospective analyses by clarifying multi-source case ascertainment logic (ICD-10, ICD-9, and self-report) under unified earliest-date rules.
- Improved practical support for diabetes endpoint setup across predefined disease definitions (`Diabetes`, `T1DM`, `T2DM`) in cohort construction workflows.

### Machine learning
- Added new exported helper `ukb_ml_split_data()` for train/internal-validation splitting.
- Supports optional stratified sampling by a categorical variable and reproducible splitting with `seed`.
- Added manual page `man/ukb_ml_split_data.Rd` and NAMESPACE export.

## UKBAnalytica 0.6.1 (2026-04-07)
add sensitivity analysis module and refine the docs.
- add `select_incident_by_years()` utility to split incident cases within or after a year cutoff from enrollment.

## UKBAnalytica 0.6.0 (2026-04-02)

### New Module: Machine Learning

#### Core ML Functions (`ml_model.R`)
- `ukb_ml_model()`: Unified interface for training ML models
  - Random Forest (`ranger`)
  - XGBoost (`xgboost`)
  - Elastic Net (`glmnet`)
  - SVM (`e1071`)
  - Neural Network (`nnet`)
  - Logistic/Linear regression
- `ukb_ml_predict()`: Generate predictions
- `ukb_ml_cv()`: K-fold cross-validation with optional repeats
- `ukb_ml_compare()`: Compare multiple models
- `ukb_ml_importance()`: Extract variable importance

#### Model Evaluation (`ml_evaluate.R`)
- `ukb_ml_metrics()`: Compute performance metrics (AUC, accuracy, etc.)
- `ukb_ml_roc()`: ROC curve analysis with CI
- `ukb_ml_calibration()`: Calibration curve with Brier score and ECE
- `ukb_ml_confusion()`: Confusion matrix

#### SHAP Interpretation (`ml_shap.R`)
- `ukb_shap()`: Compute SHAP values for model interpretation
- `ukb_shap_summary()`: Feature importance from SHAP
- `ukb_shap_dependence()`: Single feature SHAP analysis
- `ukb_shap_force()`: Single observation explanation

#### Survival ML (`ml_survival.R`)
- `ukb_ml_survival()`: Survival machine learning models
  - Random Survival Forest (`randomForestSRC`)
  - GBM Survival (`gbm`)
  - CoxNet regularized Cox (`glmnet`)
- `ukb_ml_survival_predict()`: Survival probability prediction
- `ukb_ml_survival_importance()`: Variable importance
- `ukb_ml_survival_shap()`: SHAP for survival models

#### Visualization
- `plot_ml_importance()`: Variable importance bar/dot plot
- `plot_ml_roc()`: ROC curve plot
- `plot_ml_calibration()`: Calibration curve plot
- `plot_ml_confusion()`: Confusion matrix heatmap
- `plot_ml_compare()`: Model comparison plot
- `plot_shap_summary()`: SHAP beeswarm/bar plot
- `plot_shap_dependence()`: SHAP dependence plot
- `plot_shap_force()`: SHAP waterfall plot

### Dependencies (Suggests)
- Added: `ranger`, `xgboost`, `glmnet`, `e1071`, `nnet`, `fastshap`, `pROC`, `randomForestSRC`

### Documentation
- Updated Advanced Analysis chapter with ML module
- Updated README with ML examples

---

## UKBAnalytica 0.5.0 (2026-04-01)

### New Modules: Advanced Statistical Analysis

#### Subgroup Analysis (`subgroup.R`)
- `run_subgroup_analysis()`: Stratified analysis with interaction p-values
- `run_multi_subgroup()`: Batch analysis across multiple subgroup variables
- Supports Cox, logistic, and linear regression models

#### Propensity Score Methods (`propensity.R`)
- `estimate_propensity_score()`: PS estimation via logistic regression or GBM
- `match_propensity()`: 1:k nearest neighbor matching with caliper
- `calculate_weights()`: IPTW weights (ATE, ATT, ATC)
- `assess_balance()`: Covariate balance assessment with SMD
- `run_weighted_analysis()`: Weighted regression analysis

#### Mediation Analysis (`mediation.R`)
- `run_mediation()`: Causal mediation analysis (wrapping regmedint)
- `run_multi_mediator()`: Test multiple mediators
- Supports linear, logistic, and Cox outcome models
- Effects: CDE, PNDE, TNIE, TE, Proportion Mediated

#### Multiple Imputation Pooling (`mi_pool.R`)
- `pool_mi_models()`: Combine regression results using Rubin's Rules
- `fit_mi_models()`: Fit models across imputed datasets
- `create_imputation_list()`: Convert to mitools imputationList
- `pool_custom_estimates()`: Pool custom statistics
- Supports lm, logistic, poisson, cox, negbin models
- Reports FMI (Fraction of Missing Information)

#### Visualization (`visualization.R`)
- `plot_forest()`: Forest plots for subgroup/regression results
- `plot_km_curve()`: Kaplan-Meier survival curves
- `plot_ps_distribution()`: Propensity score distribution (histogram/density)
- `plot_balance()`: Covariate balance before/after matching
- `plot_calibration()`: Calibration plots
- `plot_mediation()`: Mediation effect plots (bar, decomposition, path diagram)
- `plot_mediation_forest()`: Multi-mediator forest plot
- `plot_mi_pooled()`: MI pooled results forest plot
- `plot_mi_diagnostics()`: FMI and variance diagnostics

### Documentation
- New chapter: Advanced Analysis Modules (`docs/08-advanced-analysis.Rmd`)
- Updated technical design document with all module specifications
- Updated README with advanced analysis examples

### Dependencies (Suggests)
- Added: `MatchIt`, `gbm`, `regmedint`, `mitools`, `MASS`, `cobalt`

## UKBAnalytica 0.4.0 (2026-02-05)
Fix bug in `survival.R`: person who has primary disease before initial time will be set `NA` in survival time (in order to distinguish it from person who has primary disease after initial time, with `non-NA` survival time). 

## UKBAnalytica 0.3.0 (2026-01-26)

Add `variable_preprocess.R` module for preprocessing baseline variables.

## UKBAnalytica 0.2.0 (2026-01-25)

### Major changes
- Refactored the primary analysis interface to return a cohort-retaining **wide-format** dataset by default.
- Added a `primary_disease` argument to compute `outcome_status` and `outcome_surv_time` for a single primary endpoint.
- Added `prevalent_sources` and `outcome_sources` argument into `build_survival_dataset` function to manage self-report bias.

### New features
- Multi-source phenotyping support with configurable `sources` (ICD-10, ICD-9, self-report, death).
- Cohort-level follow-up time computation with administrative censoring and death censoring.
- Expanded predefined disease definitions to cover common conditions for rapid prototyping.

### Data acquisition (RAP)
- Added Python utilities under `inst/python/` to extract:
  - Demographic fields (user-specified UKB field IDs; optional ID file input).
  - Metabolomics (all fields; plus non-ratio subset driven by `inst/extdata/metabolites_non_ratio.txt`).
  - Proteomics (batch extract with optional merge).

### Documentation
- Updated README to prioritize data acquisition (RAP) before survival endpoint construction.
- Added a package overview figure in `man/figures/`.

## UKBAnalytica 0.1.0

- Initial release: parsing UKB RAP exports and generating survival-analysis-ready datasets.
