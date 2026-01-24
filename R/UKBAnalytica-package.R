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
#' @importFrom stringi stri_extract_all_regex stri_trim_both
#'
"_PACKAGE"

# Global variables declaration to pass R CMD check
utils::globalVariables(c(
  ".", "eid", "idx", "icd10_code", "icd9_code", "sr_code", "death_code",
  "diag_date", "diag_year", "death_date", "source", "disease",
  "earliest_date", "col_name", "date_col", "col", "instance", "array_idx",
  "cause_type", "baseline_date", "end_date", "status", "surv_time",
  "prevalent_case", "diagnosis_source", "p41270", "p41271",
  "default_surv_time", "control_surv_time", "primary_outcome_prevalent",
  "outcome_status", "outcome_surv_time", "outcome_prevalent"
))
