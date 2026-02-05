# Suppress R CMD check NOTEs for data.table variables
if(getRversion() >= "2.15.1") utils::globalVariables(c(
  "i.status", "i.prevalent_case", "i.earliest_date", "i.diagnosis_source", "i.surv_time"
))

#' @title Build Survival Analysis Dataset
#'
#' @description
#' Integrates diagnosis data from multiple sources (ICD-10, ICD-9, self-report, death)
#' to generate a survival dataset. By default, returns a wide table that retains
#' all participants and adds disease history/incident indicators plus follow-up
#' time for a primary disease.
#'
#' @param dt A data.table or data.frame containing complete UKB data.
#' @param disease_definitions Named list of disease definitions (see \code{\link{create_disease_definition}}).
#' @param prevalent_sources Character vector specifying data sources for identifying
#'   prevalent (baseline) cases. Self-report is recommended here since participants
#'   reporting a disease at baseline clearly had it before enrollment. Default includes
#'   all sources: "ICD10", "ICD9", "Self-report", "Death".
#' @param outcome_sources Character vector specifying data sources for defining
#'   incident outcomes. Self-report is typically excluded here because self-reported
#'   diagnosis dates are imprecise (year only) and less reliable for prospective
#'   endpoint ascertainment. Default: "ICD10", "ICD9", "Death".
#' @param censor_date Administrative censoring date (default: "2023-10-31").
#' @param baseline_col Column name for baseline assessment date (default: "p53_i0").
#' @param primary_disease Disease key used to compute follow-up time and event
#'   status (must be in \code{disease_definitions}). If NULL, the first disease
#'   in the list is used.
#' @param output Output format: \code{"wide"} (default) returns the original data
#'   with disease indicator columns; \code{"long"} returns case-level records.
#' @param include_all Logical; when \code{output = "long"}, if TRUE includes the full
#'   cohort with non-cases coded as status = 0.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{<Disease>_history}{1 if prevalent case (from prevalent_sources), 0 otherwise}
#'     \item{<Disease>_incident}{1 if incident case (from outcome_sources), 0 otherwise}
#'     \item{outcome_status}{Event indicator for primary disease (1=event, 0=censored, NA=prevalent case)}
#'     \item{outcome_surv_time}{Follow-up time in years for primary disease (NA for prevalent cases)}
#'   }
#'
#' @details
#' This function supports separate source definitions for prevalent case exclusion
#' and outcome ascertainment. This is important because:
#' \itemize{
#'   \item Self-reported conditions at baseline clearly indicate pre-existing disease
#'     and should be used for prevalent case identification.
#'   \item However, self-reported incident events during follow-up have imprecise dates
#'     (year only) and lower validity, making them unsuitable for outcome definition.
#' }
#'
#' Case classification logic:
#' \itemize{
#'   \item \strong{Prevalent case}: Earliest diagnosis date (from prevalent_sources) <= baseline date.
#'     These participants have \code{outcome_status = NA} and \code{outcome_surv_time = NA}
#'     because they are not at risk for incident disease.
#'   \item \strong{Incident case}: Earliest diagnosis date (from outcome_sources) > baseline date
#'   \item \strong{Censored}: No diagnosis by end of follow-up (status = 0)
#' }
#'
#' Follow-up time calculation (controlled by primary_disease):
#' \itemize{
#'   \item Prevalent case (primary disease): NA (not at risk)
#'   \item Incident case: (diagnosis_date - baseline_date) / 365.25
#'   \item Censored: (min(death_date, censor_date) - baseline_date) / 365.25
#' }
#'
#' @examples
#' \dontrun{
#' ukb_data <- data.table::fread("ukb_data.csv")
#' diseases <- get_predefined_diseases()[c("AA", "Hypertension")]
#'
#' # Use self-report for prevalent exclusion, but not for outcome
#' surv_data <- build_survival_dataset(
#'   ukb_data,
#'   diseases,
#'   prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
#'   outcome_sources = c("ICD10", "Death"),
#'   primary_disease = "AA"
#' )
#' }
#'
#' @import data.table
#' @export
build_survival_dataset <- function(dt,
                                    disease_definitions,
                                    prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
                                    outcome_sources = c("ICD10", "ICD9", "Death"),
                                    censor_date = as.Date("2023-10-31"),
                                    baseline_col = "p53_i0",
                                    primary_disease = NULL,
                                    output = c("wide", "long"),
                                    include_all = TRUE) {

  output <- match.arg(output)

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  if (!baseline_col %in% names(dt)) {
    stop(sprintf("Baseline column not found: %s", baseline_col))
  }

  message("[build_survival_dataset] Extracting diagnosis records...")
  message(sprintf("  Prevalent sources: %s", paste(prevalent_sources, collapse = ", ")))
  message(sprintf("  Outcome sources: %s", paste(outcome_sources, collapse = ", ")))

  # Extract cases for prevalent identification (includes self-report)
  prevalent_cases_dt <- extract_cases_by_source(
    dt = dt,
    disease_definitions = disease_definitions,
    sources = prevalent_sources,
    censor_date = censor_date,
    baseline_col = baseline_col
  )

  # Extract cases for outcome (excludes self-report by default)
  outcome_cases_dt <- extract_cases_by_source(
    dt = dt,
    disease_definitions = disease_definitions,
    sources = outcome_sources,
    censor_date = censor_date,
    baseline_col = baseline_col
  )

  if (output == "long" && !include_all) {
    message("[build_survival_dataset] Complete")
    return(outcome_cases_dt[, .(eid, disease, status, surv_time, prevalent_case, earliest_date, diagnosis_source)])
  }

  message("[build_survival_dataset] Calculating survival times...")

  baseline_dt <- dt[, .(eid, baseline_date = as.Date(get(baseline_col)))]
  death_dates <- get_death_dates(dt)
  all_eids <- data.table::merge.data.table(baseline_dt, death_dates, by = "eid", all.x = TRUE)

  all_eids[, end_date := pmin(death_date, censor_date, na.rm = TRUE)]
  all_eids[is.na(end_date), end_date := censor_date]
  all_eids[, default_surv_time := as.numeric(end_date - baseline_date) / 365.25]

  diseases <- names(disease_definitions)
  if (length(diseases) == 0) {
    stop("No disease definitions provided")
  }

  if (is.null(primary_disease)) {
    primary_disease <- diseases[[1]]
  }
  if (!primary_disease %in% diseases) {
    stop(sprintf("Primary disease '%s' not found in disease definitions", primary_disease))
  }

  if (output == "long") {
    full_list <- lapply(diseases, function(d) {
      # Use prevalent_sources for history, outcome_sources for incident
      d_prevalent <- prevalent_cases_dt[disease == d]
      d_outcome <- outcome_cases_dt[disease == d]

      cohort <- data.table::copy(all_eids)
      cohort[, `:=`(
        disease = d,
        status = 0L,
        prevalent_case = FALSE,
        earliest_date = as.Date(NA),
        diagnosis_source = NA_character_,
        surv_time = default_surv_time
      )]

      # Mark prevalent cases (from prevalent_sources including self-report)
      if (nrow(d_prevalent) > 0) {
        prevalent_eids <- d_prevalent[prevalent_case == TRUE, eid]
        cohort[eid %in% prevalent_eids, prevalent_case := TRUE]
      }

      # Mark incident cases and update survival time (from outcome_sources)
      if (nrow(d_outcome) > 0) {
        cohort[d_outcome,
          `:=`(
            status = i.status,
            earliest_date = i.earliest_date,
            diagnosis_source = i.diagnosis_source,
            surv_time = i.surv_time
          ),
          on = "eid"
        ]
      }

      cohort[, c("baseline_date", "death_date", "end_date", "default_surv_time") := NULL]
      cohort
    })

    full_cohort <- data.table::rbindlist(full_list, use.names = TRUE, fill = TRUE)
    full_cohort <- full_cohort[!is.na(surv_time) & surv_time >= 0]
    data.table::setorder(full_cohort, disease, eid)

    message("[build_survival_dataset] Complete")

    return(full_cohort[, .(eid, disease, status, surv_time, prevalent_case, earliest_date, diagnosis_source)])
  }

  # Generate wide format: history from prevalent_sources, incident from outcome_sources
  wide_dt <- generate_wide_format_dual_source(
    dt,
    disease_definitions,
    prevalent_sources = prevalent_sources,
    outcome_sources = outcome_sources,
    censor_date = censor_date,
    baseline_col = baseline_col
  )

  primary_outcome_cases <- outcome_cases_dt[disease == primary_disease]
  primary_prevalent_cases <- prevalent_cases_dt[disease == primary_disease & prevalent_case == TRUE]
  
  # Get all prevalent eids for primary disease
  prevalent_eids <- character(0)
  if (nrow(primary_prevalent_cases) > 0) {
    prevalent_eids <- primary_prevalent_cases$eid
  }
  
  outcome_dt <- data.table::copy(all_eids)
  outcome_dt[, `:=`(
    outcome_status = 0L,
    outcome_surv_time = default_surv_time
  )]

  # Set outcome status and time from outcome_sources only
  # Only for participants who are NOT prevalent cases
  if (nrow(primary_outcome_cases) > 0) {
    non_prevalent_outcomes <- primary_outcome_cases[!eid %in% prevalent_eids]
    
    if (nrow(non_prevalent_outcomes) > 0) {
      outcome_dt[non_prevalent_outcomes,
        `:=`(
          outcome_status = i.status,
          outcome_surv_time = i.surv_time
        ),
        on = "eid"
      ]
    }
  }
  
  # CRITICAL: Set outcome_status = NA and outcome_surv_time = NA for prevalent cases
  # These participants had the primary disease at or before baseline, 

  # so they are not at risk for incident disease and should be excluded from survival analysis
  if (length(prevalent_eids) > 0) {
    outcome_dt[eid %in% prevalent_eids, `:=`(
      outcome_status = NA_integer_,
      outcome_surv_time = NA_real_
    )]
  }

  outcome_dt[, c("baseline_date", "death_date", "end_date", "default_surv_time") := NULL]

  result_dt <- data.table::merge.data.table(
    data.table::copy(dt),
    wide_dt,
    by = "eid",
    all.x = TRUE
  )
  result_dt <- data.table::merge.data.table(result_dt, outcome_dt, by = "eid", all.x = TRUE)
  data.table::setorder(result_dt, eid)

  message("[build_survival_dataset] Complete")

  return(result_dt)
}


#' @title Build Full Cohort Survival Dataset
#'
#' @description
#' Extends \code{\link{build_survival_dataset}} to include non-cases (controls)
#' for each disease, creating a complete cohort for survival analysis.
#'
#' @inheritParams build_survival_dataset
#' @param exclude_prevalent Logical; if TRUE, excludes prevalent cases from output.
#'
#' @return A data.table with complete cohort survival data.
#'
#' @keywords internal
build_full_cohort <- function(dt,
                               disease_definitions,
                               prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
                               outcome_sources = c("ICD10", "ICD9", "Death"),
                               censor_date = as.Date("2023-10-31"),
                               baseline_col = "p53_i0",
                               primary_disease = NULL,
                               exclude_prevalent = TRUE) {

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  full_cohort <- build_survival_dataset(
    dt,
    disease_definitions,
    prevalent_sources = prevalent_sources,
    outcome_sources = outcome_sources,
    censor_date = censor_date,
    baseline_col = baseline_col,
    primary_disease = primary_disease,
    output = "long",
    include_all = TRUE
  )

  if (exclude_prevalent) {
    full_cohort <- full_cohort[prevalent_case == FALSE | is.na(prevalent_case)]
  }

  full_cohort <- full_cohort[!is.na(surv_time) & surv_time > 0]
  data.table::setorder(full_cohort, disease, eid)

  return(full_cohort)
}
