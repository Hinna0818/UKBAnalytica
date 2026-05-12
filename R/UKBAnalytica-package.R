#' UKBAnalytica: UK Biobank Data Processing and Survival Analysis Toolkit
#'
#' @description
#' A high-performance R package for processing UK Biobank (UKB) Research Analysis
#' Platform (RAP) data exports. Designed for epidemiological studies requiring
#' efficient extraction of diagnosis records and generation of survival analysis datasets.
#'
#' @details
#' \strong{Core Capabilities:}
#' \itemize{
#'   \item Parse ICD-10/ICD-9 diagnosis codes from mixed-format data
#'   \item Parse OPCS4 operative procedure codes from hospital summary operations
#'   \item Process self-reported illness data with fractional year conversion
#'   \item Integrate death registry data as diagnosis sources
#'   \item Generate Cox regression-ready survival datasets
#'   \item Support flexible data source selection for sensitivity analyses
#' }
#'
#' \strong{Key Functions:}
#' \itemize{
#'   \item \code{\link{parse_icd10_diagnoses}}: Extract ICD-10 hospital diagnoses
#'   \item \code{\link{parse_icd9_diagnoses}}: Extract ICD-9 hospital diagnoses
#'   \item \code{\link{parse_opcs4_procedures}}: Extract OPCS4 hospital procedures
#'   \item \code{\link{parse_self_reported_illnesses}}: Extract self-reported conditions
#'   \item \code{\link{parse_death_records}}: Extract death registry data
#'   \item \code{\link{build_survival_dataset}}: Generate survival analysis data
#'   \item \code{\link{extract_cases_by_source}}: Flexible source-specific extraction
#' }
#'
#' \strong{UKB Data Fields:}
#' \itemize{
#'   \item ICD-10: \code{p41270} (codes) + \code{p41280_a*} (dates)
#'   \item ICD-9: \code{p41271} (codes) + \code{p41281_a*} (dates)
#'   \item OPCS4: \code{p41272} (codes) + \code{p41282_a*} (dates)
#'   \item Self-report: \code{p20002_i*_a*} (codes) + \code{p20008_i*_a*} (years)
#'   \item Death: \code{p40001/p40002} (causes) + \code{p40000} (dates)
#' }
#'
#' @references
#' UK Biobank Data Showcase: \url{https://biobank.ndph.ox.ac.uk/showcase/}
#'
#' @docType package
#' @name UKBAnalytica
#' @aliases UKBAnalytica-package
#'
#' @import data.table
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom stringi stri_extract_all_regex stri_trim_both
#' @importFrom stats quantile relevel sd terms model.matrix predict approx
#' @importFrom utils head modifyList
#' @importFrom grid unit
#'
"_PACKAGE"

# Global variables declaration (centralized)
# All NSE variables used in data.table and ggplot2 expressions
utils::globalVariables(c(
  # Core parsing module
  ".", "eid", "idx", "icd10_code", "icd9_code", "opcs4_code", "sr_code", "death_code",
  "diag_date", "diag_year", "death_date", "source", "disease",
  "earliest_date", "col_name", "date_col", "col", "instance", "array_idx",
  "cause_type", "baseline_date", "end_date", "status", "surv_time",
  "prevalent_case", "diagnosis_source", "p41270", "p41271", "p41272",
  "default_surv_time", "control_surv_time", "primary_outcome_prevalent",
  "outcome_status", "outcome_surv_time", "outcome_prevalent",
  # Algorithm and spirometry module
  "algo_date", "algo_source",
  "COPD_prevalent", "COPD_prevalent_icd", "COPD_prevalent_spirometry",
  "COPD_prevalent_source", "COPD_incident", "COPD_spirometry", "COPD_gold_grade",
  "T1DM_history", "T2DM_history", "T2DM_history_enhanced", "Diabetes_history",
  "diabetes_hba1c", "diabetes_subtype", "hba1c_value",
  "default_end", "incident_surv_time",
  "fev1", "fvc", "fev1_fvc_ratio", "fev1_predicted_pct", "obstruction",
  # Survival module

  "i.status", "i.prevalent_case", "i.earliest_date", "i.diagnosis_source", "i.surv_time",
  # Subgroup/propensity module
  "subgroup_var_name", "subgroup_level", "n_total", "n_event",
  "ps", "weight", "match_id", "match_distance", "..covariates",
  "mean_treated", "mean_control", "smd", "balanced",
  "subgroup", "estimate", "lower95", "upper95", "ci_label", "x_pos",
  "ymin", "ymax", "time", "surv", "strata", "n.risk", "n.event",
  "group", "method", "smd_before", "smd_after",
  "observed", "predicted", "bin", "lower", "upper", "count",
  # Correlation module
  "Var1", "Var2", "Correlation", "label",
  # Variable preprocess module
  "fruit_score", "vegetable_score", "fish_score", "meat_score", "diet_score",
  "p1289_i0", "p1299_i0", "p1309_i0", "p1319_i0", "p1329_i0", "p1339_i0",
  # Mediation module
  "effect", "effect_label", "color_group", "est", "p", "p_label", "neg_log10_p",
  "effect_type", "value", "abs_value", "proportion", "mediator", "se",
  "pvalue", "sig", "name", "tnie", "pnde", "te", "pm", "exp_est", "exp_lower", "exp_upper",
  # Multiple imputation module
  "term", "std.error", "statistic", "conf.low", "conf.high", "fmi",
  "fmi_pct", "level",
  # Visualization module
  "after_stat", "x", "y", "xend", "yend",
  # Machine learning module
  "importance", "variable", "sensitivity", "specificity", "model",
  "feature", "mean_abs_shap", "shap", "value_norm", "feature_value",
  "shap_value", "color_value", "direction", "feature_label", "metric",
  "Predicted", "Actual", "Freq", "Proportion", "Label"
))
