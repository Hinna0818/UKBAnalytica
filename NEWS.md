# UKBAnalytica News

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
- `run_sensitivity_mediation()`: Sensitivity analysis for unmeasured confounding
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
- Added Python utilities under `inst/python/` to download:
  - Demographic fields (user-specified UKB field IDs; optional ID file input).
  - Metabolomics (all fields; plus non-ratio subset driven by `inst/extdata/metabolites_non_ratio.txt`).
  - Proteomics (batch download with optional merge).

### Documentation
- Updated README to prioritize data acquisition (RAP) before survival endpoint construction.
- Added a package overview figure in `man/figures/`.

## UKBAnalytica 0.1.0

- Initial release: parsing UKB RAP exports and generating survival-analysis-ready datasets.
