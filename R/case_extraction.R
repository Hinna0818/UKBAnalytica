#' @title Extract Cases by Specified Data Sources
#'
#' @description
#' Flexibly extracts disease cases using user-specified data sources.
#' Enables main analysis with strict case definitions (e.g., ICD-10 only)
#' and sensitivity analyses with broader definitions (e.g., all sources).
#'
#' @param dt A data.table or data.frame containing complete UKB data.
#' @param disease_definitions Named list of disease definitions.
#' @param sources Character vector specifying data sources to include.
#'   Valid options: "ICD10", "ICD9", "Self-report", "Death".
#' @param censor_date Administrative censoring date.
#' @param baseline_col Column name for baseline assessment date.
#'
#' @return A data.table with case-level survival data from specified sources.
#'
#' @details
#' This function is designed for epidemiological studies requiring:
#' \itemize{
#'   \item Main analysis with hospital-confirmed diagnoses only
#'   \item Sensitivity analyses including self-reported conditions
#'   \item Source-specific case counts for methods reporting
#' }
#'
#' @examples
#' \dontrun{
#' diseases <- get_predefined_diseases()[c("AA", "Hypertension")]
#'
#' # Main analysis: ICD-10 only
#' main <- extract_cases_by_source(dt, diseases, sources = "ICD10")
#'
#' # Sensitivity: All sources
#' sens <- extract_cases_by_source(dt, diseases,
#'                                  sources = c("ICD10", "ICD9", "Self-report"))
#' }
#'
#' @import data.table
#' @export
extract_cases_by_source <- function(dt,
                                     disease_definitions,
                                     sources = c("ICD10", "ICD9", "Self-report", "Death"),
                                     censor_date = as.Date("2023-10-31"),
                                     baseline_col = "p53_i0") {

  valid_sources <- c("ICD10", "ICD9", "Self-report", "Death")
  sources <- match.arg(sources, valid_sources, several.ok = TRUE)

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  if (!baseline_col %in% names(dt)) {
    stop(sprintf("Baseline column not found: %s", baseline_col))
  }

  message(sprintf("[extract_cases_by_source] Using sources: %s", paste(sources, collapse = ", ")))

  # Extract only requested sources
  icd10_long <- if ("ICD10" %in% sources) parse_icd10_diagnoses(dt) else data.table::data.table()
  icd9_long <- if ("ICD9" %in% sources) parse_icd9_diagnoses(dt) else data.table::data.table()
  sr_long <- if ("Self-report" %in% sources) parse_self_reported_illnesses(dt, baseline_col) else data.table::data.table()
  death_long <- if ("Death" %in% sources) parse_death_records(dt) else data.table::data.table()

  death_dates <- get_death_dates(dt)
  baseline_dt <- dt[, .(eid, baseline_date = as.Date(get(baseline_col)))]

  # Process each disease
  results_list <- lapply(names(disease_definitions), function(disease_key) {
    def <- disease_definitions[[disease_key]]
    diagnosis_sources <- list()

    if ("ICD10" %in% sources && !is.null(def$icd10_pattern) && nrow(icd10_long) > 0) {
      filtered <- filter_icd10_codes(icd10_long, def$icd10_pattern, disease_key)
      if (nrow(filtered) > 0) diagnosis_sources$icd10 <- aggregate_icd10_earliest(filtered)
    }

    if ("ICD9" %in% sources && !is.null(def$icd9_pattern) && nrow(icd9_long) > 0) {
      filtered <- filter_icd9_codes(icd9_long, def$icd9_pattern, disease_key)
      if (nrow(filtered) > 0) diagnosis_sources$icd9 <- aggregate_icd9_earliest(filtered)
    }

    if ("Self-report" %in% sources && !is.null(def$sr_codes) && length(def$sr_codes) > 0 && nrow(sr_long) > 0) {
      filtered <- filter_self_report_codes(sr_long, def$sr_codes, disease_key)
      if (nrow(filtered) > 0) diagnosis_sources$sr <- aggregate_self_report_earliest(filtered)
    }

    if ("Death" %in% sources && !is.null(def$icd10_pattern) && nrow(death_long) > 0) {
      filtered <- filter_death_codes(death_long, def$icd10_pattern, disease_key)
      if (nrow(filtered) > 0) diagnosis_sources$death <- aggregate_death_as_diagnosis(filtered)
    }

    if (length(diagnosis_sources) == 0) return(NULL)

    all_diagnoses <- data.table::rbindlist(diagnosis_sources, use.names = TRUE, fill = TRUE)

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
    warning("[extract_cases_by_source] No cases found")
    return(data.table::data.table(
      eid = integer(0), disease = character(0), earliest_date = as.Date(character(0)),
      diagnosis_source = character(0), prevalent_case = logical(0),
      status = integer(0), surv_time = numeric(0)
    ))
  }

  diagnosis_dt <- data.table::rbindlist(results_list, use.names = TRUE, fill = TRUE)

  # Calculate survival metrics
  surv_dt <- data.table::merge.data.table(diagnosis_dt, baseline_dt, by = "eid", all.x = TRUE)
  surv_dt <- data.table::merge.data.table(surv_dt, death_dates, by = "eid", all.x = TRUE)

  surv_dt[, prevalent_case := !is.na(earliest_date) & earliest_date <= baseline_date]
  surv_dt[, status := as.integer(!is.na(earliest_date) & earliest_date > baseline_date)]

  surv_dt[, end_date := data.table::fifelse(
    status == 1L, earliest_date,
    pmin(death_date, censor_date, na.rm = TRUE)
  )]
  surv_dt[is.na(end_date), end_date := censor_date]
  surv_dt[, surv_time := as.numeric(end_date - baseline_date) / 365.25]
  surv_dt[surv_time < 0, `:=`(surv_time = NA_real_, status = NA_integer_)]

  surv_dt[, c("baseline_date", "death_date", "end_date") := NULL]
  data.table::setorder(surv_dt, disease, eid)

  return(surv_dt)
}


#' @title Generate Wide-Format Disease Status Table
#'
#' @description
#' Transforms case-level data into a wide-format table with one row per participant.
#' Each disease generates two columns: \code{_history} (prevalent) and \code{_incident}.
#'
#' @inheritParams extract_cases_by_source
#' @param include_dates Logical; if TRUE, includes diagnosis date columns.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{{Disease}_history}{1 if prevalent case, 0 otherwise (covariate use)}
#'     \item{{Disease}_incident}{1 if incident case, 0 otherwise (outcome use)}
#'     \item{{Disease}_date}{(Optional) Earliest diagnosis date}
#'   }
#'
#' @examples
#' \dontrun{
#' diseases <- get_predefined_diseases()[c("AA", "Hypertension", "Diabetes")]
#'
#' # For Cox regression:
#' # - Use _history columns as covariates
#' # - Use _incident column of primary outcome as event indicator
#' wide_dt <- generate_wide_format(dt, diseases, sources = "ICD10")
#' }
#'
#' @export
generate_wide_format <- function(dt,
                                  disease_definitions,
                                  sources = c("ICD10", "ICD9", "Self-report", "Death"),
                                  censor_date = as.Date("2023-10-31"),
                                  baseline_col = "p53_i0",
                                  include_dates = FALSE) {

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  surv_long <- extract_cases_by_source(dt, disease_definitions, sources, censor_date, baseline_col)
  all_eids <- dt[, .(eid)]
  wide_dt <- data.table::copy(all_eids)

  diseases <- unique(surv_long$disease)

  for (d in diseases) {
    d_data <- surv_long[disease == d]

    d_wide <- data.table::merge.data.table(
      all_eids,
      d_data[, .(eid, prevalent_case, status, earliest_date)],
      by = "eid", all.x = TRUE
    )

    d_wide[, (paste0(d, "_history")) := as.integer(prevalent_case == TRUE & !is.na(prevalent_case))]
    d_wide[, (paste0(d, "_incident")) := as.integer(status == 1L & !is.na(status))]

    if (include_dates) {
      d_wide[, (paste0(d, "_date")) := earliest_date]
    }

    cols_to_keep <- c("eid", paste0(d, "_history"), paste0(d, "_incident"))
    if (include_dates) cols_to_keep <- c(cols_to_keep, paste0(d, "_date"))

    d_wide <- d_wide[, cols_to_keep, with = FALSE]
    wide_dt <- data.table::merge.data.table(wide_dt, d_wide, by = "eid", all.x = TRUE)
  }

  # Replace NA with 0 for binary columns
  for (col in grep("_(history|incident)$", names(wide_dt), value = TRUE)) {
    data.table::set(wide_dt, which(is.na(wide_dt[[col]])), col, 0L)
  }

  data.table::setorder(wide_dt, eid)
  return(wide_dt)
}


#' @title Compare Case Counts Across Data Sources
#'
#' @description
#' Generates a summary table comparing case counts from different data sources.
#' Useful for methods sections and sensitivity analysis planning.
#'
#' @param dt A data.table containing complete UKB data.
#' @param disease_definitions Named list of disease definitions.
#' @param baseline_col Column name for baseline date.
#'
#' @return A data.table with case counts by source and combination.
#'
#' @examples
#' \dontrun{
#' diseases <- get_predefined_diseases()[c("AA", "Hypertension")]
#' comparison <- compare_data_sources(dt, diseases)
#' print(comparison)
#' }
#'
#' @export
compare_data_sources <- function(dt,
                                  disease_definitions,
                                  baseline_col = "p53_i0") {

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  diseases <- names(disease_definitions)

  comparison_list <- lapply(diseases, function(d) {
    single_def <- disease_definitions[d]

    icd10_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = "ICD10", baseline_col = baseline_col),
      error = function(e) data.table::data.table()
    )

    icd9_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = "ICD9", baseline_col = baseline_col),
      error = function(e) data.table::data.table()
    )

    sr_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = "Self-report", baseline_col = baseline_col),
      error = function(e) data.table::data.table()
    )

    hospital_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = c("ICD10", "ICD9"), baseline_col = baseline_col),
      error = function(e) data.table::data.table()
    )

    all_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = c("ICD10", "ICD9", "Self-report"), baseline_col = baseline_col),
      error = function(e) data.table::data.table()
    )

    data.table::data.table(
      disease = d,
      ICD10_total = nrow(icd10_cases),
      ICD10_incident = sum(icd10_cases$status == 1, na.rm = TRUE),
      ICD10_prevalent = sum(icd10_cases$prevalent_case == TRUE, na.rm = TRUE),
      ICD9_total = nrow(icd9_cases),
      Self_report_total = nrow(sr_cases),
      Hospital_total = nrow(hospital_cases),
      All_sources_total = nrow(all_cases),
      All_sources_incident = sum(all_cases$status == 1, na.rm = TRUE),
      All_sources_prevalent = sum(all_cases$prevalent_case == TRUE, na.rm = TRUE)
    )
  })

  result <- data.table::rbindlist(comparison_list)
  return(result)
}


#' @title Prepare Analysis-Ready Dataset with Primary Outcome
#'
#' @description
#' Generates a complete analysis-ready dataset with all diseases as covariates
#' and specified primary outcome with survival time and status.
#'
#' @param dt A data.table containing complete UKB data.
#' @param disease_definitions Named list of disease definitions.
#' @param primary_outcome Name of the primary outcome disease.
#' @param sources Character vector of data sources to use.
#' @param censor_date Administrative censoring date.
#' @param baseline_col Column name for baseline date.
#' @param exclude_prevalent_outcome Logical; if TRUE, excludes prevalent cases of primary outcome.
#'
#' @return A data.table ready for Cox regression with:
#'   \describe{
#'     \item{{Disease}_history}{Covariate columns for all diseases}
#'     \item{{Disease}_incident}{Incident case indicators}
#'     \item{outcome_status}{Primary outcome event indicator}
#'     \item{outcome_surv_time}{Follow-up time for primary outcome}
#'   }
#'
#' @examples
#' \dontrun{
#' diseases <- get_predefined_diseases()[c("AA", "Hypertension", "Diabetes")]
#'
#' # AA as primary outcome, adjusting for hypertension and diabetes history
#' analysis_dt <- prepare_analysis_dataset(
#'   dt, diseases, primary_outcome = "AA", sources = "ICD10"
#' )
#'
#' # Cox regression
#' library(survival)
#' coxph(Surv(outcome_surv_time, outcome_status) ~
#'       Hypertension_history + Diabetes_history, data = analysis_dt)
#' }
#'
#' @export
prepare_analysis_dataset <- function(dt,
                                      disease_definitions,
                                      primary_outcome,
                                      sources = c("ICD10", "ICD9", "Self-report", "Death"),
                                      censor_date = as.Date("2023-10-31"),
                                      baseline_col = "p53_i0",
                                      exclude_prevalent_outcome = TRUE) {

  if (!primary_outcome %in% names(disease_definitions)) {
    stop(sprintf("Primary outcome '%s' not found in disease definitions", primary_outcome))
  }

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  # Generate wide format with all diseases
  wide_dt <- generate_wide_format(dt, disease_definitions, sources, censor_date, baseline_col)

  # Get survival data for primary outcome
  surv_long <- extract_cases_by_source(dt, disease_definitions, sources, censor_date, baseline_col)
  primary_data <- surv_long[disease == primary_outcome, .(eid, status, surv_time, prevalent_case)]

  # Get default survival time for non-cases
  death_dates <- get_death_dates(dt)
  baseline_dt <- dt[, .(eid, baseline_date = as.Date(get(baseline_col)))]
  control_info <- data.table::merge.data.table(baseline_dt, death_dates, by = "eid", all.x = TRUE)
  control_info[, end_date := pmin(death_date, censor_date, na.rm = TRUE)]
  control_info[is.na(end_date), end_date := censor_date]
  control_info[, control_surv_time := as.numeric(end_date - baseline_date) / 365.25]

  # Merge primary outcome data
  wide_dt <- data.table::merge.data.table(wide_dt, primary_data, by = "eid", all.x = TRUE)
  wide_dt <- data.table::merge.data.table(
    wide_dt, control_info[, .(eid, control_surv_time)], by = "eid", all.x = TRUE
  )

  # Set outcome columns
  wide_dt[, outcome_status := data.table::fifelse(is.na(status), 0L, status)]
  wide_dt[, outcome_surv_time := data.table::fifelse(is.na(surv_time), control_surv_time, surv_time)]
  wide_dt[, outcome_prevalent := data.table::fifelse(is.na(prevalent_case), FALSE, prevalent_case)]

  # Clean up
  wide_dt[, c("status", "surv_time", "prevalent_case", "control_surv_time") := NULL]

  # Exclude prevalent cases of primary outcome
  if (exclude_prevalent_outcome) {
    n_before <- nrow(wide_dt)
    wide_dt <- wide_dt[outcome_prevalent == FALSE]
    n_excluded <- n_before - nrow(wide_dt)
    if (n_excluded > 0) {
      message(sprintf("[prepare_analysis_dataset] Excluded %d prevalent cases", n_excluded))
    }
  }

  wide_dt[, outcome_prevalent := NULL]
  wide_dt <- wide_dt[!is.na(outcome_surv_time) & outcome_surv_time > 0]
  data.table::setorder(wide_dt, eid)

  return(wide_dt)
}


#' @title Extract Disease History (Prevalent Cases) for Covariates
#'
#' @description
#' Extracts prevalent case status (diagnosed before baseline) for specified diseases.
#' Designed for use as covariates in sensitivity analyses or covariate adjustment.
#' Returns a wide-format table with one binary column per disease.
#'
#' @param dt A data.table or data.frame containing complete UKB data.
#' @param diseases Character vector of disease names to extract.
#'   Must match keys in \code{disease_definitions} or predefined disease names.
#' @param disease_definitions Named list of disease definitions. If NULL,
#'   uses \code{\link{get_predefined_diseases}}.
#' @param sources Character vector specifying data sources.
#'   Default: "ICD10". Options: "ICD10", "ICD9", "Self-report", "Death".
#' @param baseline_col Column name for baseline assessment date.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{{Disease}_history}{1 if prevalent case, 0 otherwise (one column per disease)}
#'   }
#'
#' @details
#' This function is specifically designed for extracting covariate data in
#' epidemiological studies. Common use cases:
#' \itemize{
#'   \item Adjusting for baseline comorbidities in Cox regression
#'   \item Sensitivity analyses with different case definitions
#'   \item Creating propensity score matching variables
#' }
#'
#' The function only returns history (prevalent) columns, not incident columns,
#' to clearly separate covariate extraction from outcome definition.
#'
#' @examples
#' \dontrun{
#' # Extract hypertension and diabetes history using ICD-10 only
#' history_icd10 <- extract_disease_history(
#'   dt = ukb_data,
#'   diseases = c("Hypertension", "Diabetes"),
#'   sources = "ICD10"
#' )
#'
#' # Sensitivity: Include self-reported conditions
#' history_all <- extract_disease_history(
#'   dt = ukb_data,
#'   diseases = c("Hypertension", "Diabetes"),
#'   sources = c("ICD10", "ICD9", "Self-report")
#' )
#'
#' # Merge with main analysis dataset
#' analysis_dt <- merge(outcome_data, history_icd10, by = "eid")
#' }
#'
#' @export
extract_disease_history <- function(dt,
                                     diseases,
                                     disease_definitions = NULL,
                                     sources = "ICD10",
                                     baseline_col = "p53_i0") {

  # Validate inputs
  if (!is.character(diseases) || length(diseases) == 0) {
    stop("'diseases' must be a non-empty character vector")
  }

  valid_sources <- c("ICD10", "ICD9", "Self-report", "Death")
  sources <- match.arg(sources, valid_sources, several.ok = TRUE)

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  # Use predefined diseases if not provided
 if (is.null(disease_definitions)) {
    disease_definitions <- get_predefined_diseases()
  }

  # Validate disease names
  missing_diseases <- setdiff(diseases, names(disease_definitions))
  if (length(missing_diseases) > 0) {
    stop(sprintf("Disease(s) not found in definitions: %s",
                 paste(missing_diseases, collapse = ", ")))
  }

  # Subset to requested diseases
  selected_defs <- disease_definitions[diseases]

  message(sprintf("[extract_disease_history] Extracting %d disease(s) from: %s",
                  length(diseases), paste(sources, collapse = ", ")))

  # Extract cases
  surv_long <- extract_cases_by_source(
    dt = dt,
    disease_definitions = selected_defs,
    sources = sources,
    baseline_col = baseline_col
  )

  # Get all eids
  all_eids <- dt[, .(eid)]
  result_dt <- data.table::copy(all_eids)

  # Create history column for each disease
  for (d in diseases) {
    d_data <- surv_long[disease == d & prevalent_case == TRUE, .(eid)]
    d_data[, (paste0(d, "_history")) := 1L]

    result_dt <- data.table::merge.data.table(
      result_dt,
      d_data,
      by = "eid",
      all.x = TRUE
    )

    # Fill NA with 0
    hist_col <- paste0(d, "_history")
    data.table::set(result_dt, which(is.na(result_dt[[hist_col]])), hist_col, 0L)
  }

  data.table::setorder(result_dt, eid)

  # Summary message
  for (d in diseases) {
    n_cases <- sum(result_dt[[paste0(d, "_history")]])
    message(sprintf("  %s_history: %d prevalent cases", d, n_cases))
  }

  return(result_dt)
}


#' @title Extract Disease History with Multiple Source Comparisons
#'
#' @description
#' Extracts prevalent case status from multiple data source combinations
#' simultaneously for sensitivity analysis comparison. Returns a table
#' with separate columns for each source definition.
#'
#' @param dt A data.table or data.frame containing complete UKB data.
#' @param diseases Character vector of disease names to extract.
#' @param disease_definitions Named list of disease definitions.
#' @param baseline_col Column name for baseline date.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{{Disease}_history_ICD10}{Prevalent case using ICD-10 only}
#'     \item{{Disease}_history_hospital}{Prevalent case using ICD-10 + ICD-9}
#'     \item{{Disease}_history_all}{Prevalent case using all sources}
#'   }
#'
#' @examples
#' \dontrun{
#' # Get all sensitivity variants at once
#' history_comparison <- extract_disease_history_sensitivity(
#'   dt = ukb_data,
#'   diseases = c("Hypertension", "Diabetes")
#' )
#'
#' # Compare: Hypertension_history_ICD10 vs Hypertension_history_all
#' }
#'
#' @export
extract_disease_history_sensitivity <- function(dt,
                                                 diseases,
                                                 disease_definitions = NULL,
                                                 baseline_col = "p53_i0") {

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  if (is.null(disease_definitions)) {
    disease_definitions <- get_predefined_diseases()
  }

  message("[extract_disease_history_sensitivity] Generating sensitivity variants...")

  # ICD-10 only
  hist_icd10 <- extract_disease_history(
    dt, diseases, disease_definitions,
    sources = "ICD10", baseline_col = baseline_col
  )
  old_names <- paste0(diseases, "_history")
  new_names <- paste0(diseases, "_history_ICD10")
  data.table::setnames(hist_icd10, old_names, new_names)

  # Hospital (ICD-10 + ICD-9)
  hist_hospital <- extract_disease_history(
    dt, diseases, disease_definitions,
    sources = c("ICD10", "ICD9"), baseline_col = baseline_col
  )
  new_names <- paste0(diseases, "_history_hospital")
  data.table::setnames(hist_hospital, old_names, new_names)

  # All sources
  hist_all <- extract_disease_history(
    dt, diseases, disease_definitions,
    sources = c("ICD10", "ICD9", "Self-report"), baseline_col = baseline_col
  )
  new_names <- paste0(diseases, "_history_all")
  data.table::setnames(hist_all, old_names, new_names)

  # Merge all variants
  result_dt <- data.table::merge.data.table(hist_icd10, hist_hospital, by = "eid")
  result_dt <- data.table::merge.data.table(result_dt, hist_all, by = "eid")

  # Summary
  message("\n[Summary] Prevalent case counts by source:")
  for (d in diseases) {
    n_icd10 <- sum(result_dt[[paste0(d, "_history_ICD10")]])
    n_hosp <- sum(result_dt[[paste0(d, "_history_hospital")]])
    n_all <- sum(result_dt[[paste0(d, "_history_all")]])
    message(sprintf("  %s: ICD10=%d, Hospital=%d, All=%d", d, n_icd10, n_hosp, n_all))
  }

  return(result_dt)
}

