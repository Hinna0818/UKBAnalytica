#' @title Build Survival Analysis Dataset
#'
#' @description
#' Integrates diagnosis data from multiple sources (ICD-10, ICD-9, self-report, death)
#' to generate a Cox regression-ready survival dataset with proper case classification.
#'
#' @param dt A data.table or data.frame containing complete UKB data.
#' @param disease_definitions Named list of disease definitions (see \code{\link{create_disease_definition}}).
#' @param censor_date Administrative censoring date (default: "2023-10-31").
#' @param baseline_col Column name for baseline assessment date (default: "p53_i0").
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{disease}{Disease name}
#'     \item{status}{Event indicator (1 = incident case, 0 = censored)}
#'     \item{surv_time}{Follow-up time in years}
#'     \item{prevalent_case}{TRUE if diagnosed before baseline}
#'     \item{earliest_date}{Date of first diagnosis}
#'     \item{diagnosis_source}{Data source of earliest diagnosis}
#'   }
#'
#' @details
#' Case classification logic:
#' \itemize{
#'   \item \strong{Prevalent case}: Earliest diagnosis date <= baseline date
#'   \item \strong{Incident case}: Earliest diagnosis date > baseline date (status = 1)
#'   \item \strong{Censored}: No diagnosis by end of follow-up (status = 0)
#' }
#'
#' Follow-up time calculation:
#' \itemize{
#'   \item Incident case: (diagnosis_date - baseline_date) / 365.25
#'   \item Censored: min(death_date, censor_date) - baseline_date) / 365.25
#' }
#'
#' @examples
#' \dontrun{
#' ukb_data <- data.table::fread("ukb_data.csv")
#' diseases <- get_predefined_diseases()[c("AA", "Hypertension")]
#' surv_data <- build_survival_dataset(ukb_data, diseases)
#' }
#'
#' @import data.table
#' @export
build_survival_dataset <- function(dt,
                                    disease_definitions,
                                    censor_date = as.Date("2023-10-31"),
                                    baseline_col = "p53_i0") {

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  if (!baseline_col %in% names(dt)) {
    stop(sprintf("Baseline column not found: %s", baseline_col))
  }

  message("[build_survival_dataset] Extracting diagnosis records...")

  # Extract all diagnosis sources
  icd10_long <- parse_icd10_diagnoses(dt)
  icd9_long <- parse_icd9_diagnoses(dt)
  sr_long <- parse_self_reported_illnesses(dt, baseline_col)
  death_long <- parse_death_records(dt)
  death_dates <- get_death_dates(dt)
  baseline_dt <- dt[, .(eid, baseline_date = as.Date(get(baseline_col)))]

  message("[build_survival_dataset] Processing disease definitions...")

  # Process each disease
  results_list <- lapply(names(disease_definitions), function(disease_key) {
    def <- disease_definitions[[disease_key]]
    disease_name <- if (!is.null(def$name)) def$name else disease_key

    diagnosis_sources <- list()

    # ICD-10
    if (!is.null(def$icd10_pattern) && nrow(icd10_long) > 0) {
      filtered <- filter_icd10_codes(icd10_long, def$icd10_pattern, disease_key)
      if (nrow(filtered) > 0) {
        diagnosis_sources$icd10 <- aggregate_icd10_earliest(filtered)
      }
    }

    # ICD-9
    if (!is.null(def$icd9_pattern) && nrow(icd9_long) > 0) {
      filtered <- filter_icd9_codes(icd9_long, def$icd9_pattern, disease_key)
      if (nrow(filtered) > 0) {
        diagnosis_sources$icd9 <- aggregate_icd9_earliest(filtered)
      }
    }

    # Self-report
    if (!is.null(def$sr_codes) && length(def$sr_codes) > 0 && nrow(sr_long) > 0) {
      filtered <- filter_self_report_codes(sr_long, def$sr_codes, disease_key)
      if (nrow(filtered) > 0) {
        diagnosis_sources$sr <- aggregate_self_report_earliest(filtered)
      }
    }

    # Death
    if (!is.null(def$icd10_pattern) && nrow(death_long) > 0) {
      filtered <- filter_death_codes(death_long, def$icd10_pattern, disease_key)
      if (nrow(filtered) > 0) {
        diagnosis_sources$death <- aggregate_death_as_diagnosis(filtered)
      }
    }

    if (length(diagnosis_sources) == 0) return(NULL)

    all_diagnoses <- data.table::rbindlist(diagnosis_sources, use.names = TRUE, fill = TRUE)

    # Select earliest diagnosis per participant
    earliest_per_person <- all_diagnoses[
      ,
      {
        min_idx <- which.min(earliest_date)
        list(earliest_date = earliest_date[min_idx], diagnosis_source = source[min_idx])
      },
      by = .(eid, disease)
    ]

    return(earliest_per_person)
  })

  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) {
    warning("[build_survival_dataset] No cases found for any disease")
    return(data.table::data.table(
      eid = integer(0), disease = character(0), status = integer(0),
      surv_time = numeric(0), prevalent_case = logical(0),
      earliest_date = as.Date(character(0)), diagnosis_source = character(0)
    ))
  }

  diagnosis_dt <- data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)

  message("[build_survival_dataset] Calculating survival times...")

  # Merge baseline and death dates
  surv_dt <- data.table::merge.data.table(diagnosis_dt, baseline_dt, by = "eid", all.x = TRUE)
  surv_dt <- data.table::merge.data.table(surv_dt, death_dates, by = "eid", all.x = TRUE)

  # Classify cases
  surv_dt[, prevalent_case := !is.na(earliest_date) & earliest_date <= baseline_date]
  surv_dt[, status := as.integer(!is.na(earliest_date) & earliest_date > baseline_date)]

  # Calculate follow-up time
  surv_dt[, end_date := data.table::fifelse(
    status == 1L, earliest_date,
    pmin(death_date, censor_date, na.rm = TRUE)
  )]
  surv_dt[is.na(end_date), end_date := censor_date]
  surv_dt[, surv_time := as.numeric(end_date - baseline_date) / 365.25]

  # Handle invalid values
  surv_dt[surv_time < 0, `:=`(surv_time = NA_real_, status = NA_integer_)]
  surv_dt[, c("baseline_date", "death_date", "end_date") := NULL]
  data.table::setorder(surv_dt, disease, eid)

  message("[build_survival_dataset] Complete")

  return(surv_dt[, .(eid, disease, status, surv_time, prevalent_case, earliest_date, diagnosis_source)])
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
#' @export
build_full_cohort <- function(dt,
                               disease_definitions,
                               censor_date = as.Date("2023-10-31"),
                               baseline_col = "p53_i0",
                               exclude_prevalent = TRUE) {

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  cases_dt <- build_survival_dataset(dt, disease_definitions, censor_date, baseline_col)
  all_eids <- dt[, .(eid, baseline_date = as.Date(get(baseline_col)))]
  death_dates <- get_death_dates(dt)
  all_eids <- data.table::merge.data.table(all_eids, death_dates, by = "eid", all.x = TRUE)

  diseases <- unique(cases_dt$disease)

  full_cohort_list <- lapply(diseases, function(d) {
    d_cases <- cases_dt[disease == d]
    control_eids <- setdiff(all_eids$eid, d_cases$eid)

    controls <- all_eids[eid %in% control_eids]
    controls[, `:=`(
      disease = d, status = 0L, prevalent_case = FALSE,
      earliest_date = as.Date(NA), diagnosis_source = NA_character_
    )]

    controls[, end_date := pmin(death_date, censor_date, na.rm = TRUE)]
    controls[is.na(end_date), end_date := censor_date]
    controls[, surv_time := as.numeric(end_date - baseline_date) / 365.25]
    controls[, c("baseline_date", "death_date", "end_date") := NULL]

    cohort <- data.table::rbindlist(list(d_cases, controls), use.names = TRUE, fill = TRUE)
    return(cohort)
  })

  full_cohort <- data.table::rbindlist(full_cohort_list, use.names = TRUE)

  if (exclude_prevalent) {
    full_cohort <- full_cohort[prevalent_case == FALSE | is.na(prevalent_case)]
  }

  full_cohort <- full_cohort[!is.na(surv_time) & surv_time > 0]
  data.table::setorder(full_cohort, disease, eid)

  return(full_cohort)
}
