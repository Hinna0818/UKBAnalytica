#' @title Build Survival Analysis Dataset
#'
#' @description
#' Integrates diagnosis data from multiple sources (ICD-10, ICD-9, self-report,
#' death, OPCS4 procedures, cancer registry records, First Occurrence fields,
#' algorithm) to generate a survival dataset. By default, returns a wide
#' table that retains all participants and adds disease history/incident indicators
#' plus follow-up time for a primary disease.
#'
#' @param dt A data.table or data.frame containing complete UKB data.
#' @param disease_definitions Named list of disease definitions (see \code{\link{create_disease_definition}}).
#' @param prevalent_sources Character vector specifying data sources for identifying
#'   prevalent (baseline) cases. Self-report is recommended here since participants
#'   reporting a disease at baseline clearly had it before enrollment. Default includes
#'   all core sources: "ICD10", "ICD9", "Self-report", "Death".
#'   "OPCS4" can be added for surgical phenotypes when \code{opcs4_pattern}
#'   is supplied in the disease definition.
#'   Also supports "CancerRegistry" for UKB cancer registry outcomes,
#'   "FirstOccurrence" for UKB First Occurrence fields, and "Algorithm" for
#'   UK Biobank algorithmically-defined outcomes.
#' @param outcome_sources Character vector specifying data sources for defining
#'   incident outcomes. Self-report is typically excluded here because self-reported
#'   diagnosis dates are imprecise (year only) and less reliable for prospective
#'   endpoint ascertainment. Default: "ICD10", "ICD9", "Death".
#'   "CancerRegistry" can be added for cancer outcomes; "FirstOccurrence" can
#'   be added when the extracted dataset includes UKB First Occurrence fields
#'   for the disease definition. "OPCS4" can be included when the event of
#'   interest is a surgery or procedure-based phenotype.
#' @param censor_date Administrative censoring date (default: "2023-10-31").
#' @param baseline_col Column name for baseline assessment date (default: "p53_i0").
#' @param time_skeleton Optional output from \code{\link{ukb_time_skeleton}}.
#'   When supplied, its \code{baseline_date} is used for prevalent/incident
#'   classification and its participant-specific \code{followup_end_date} is
#'   used to calculate default follow-up time for non-cases. The
#'   \code{censor_date} argument remains the common administrative censoring
#'   date used by endpoint extraction.
#' @param primary_disease One or more disease keys used to compute follow-up
#'   time and event status (all must be in \code{disease_definitions}). If
#'   \code{NULL}, the first disease in the list is used. For wide output, a
#'   single key preserves the legacy \code{outcome_status} and
#'   \code{outcome_surv_time} columns. Multiple keys produce independent
#'   \code{<Disease>_status} and \code{<Disease>_surv_time} columns in one call.
#' @param output Output format: \code{"wide"} (default) returns the original data
#'   with disease indicator columns; \code{"long"} returns case-level records.
#' @param include_all Logical; when \code{output = "long"}, if TRUE includes the full
#'   cohort with non-cases coded as status = 0.
#' @param show_flow Logical; if \code{TRUE} and \code{output = "wide"}, prints
#'   a step-by-step participant attrition table in the terminal, including
#'   counts before/after each filter and retention rates.
#' @param dt_threads Optional integer. If provided, temporarily sets
#'   \code{data.table} thread count via \code{data.table::setDTthreads()} for this
#'   function call, and restores the previous thread setting on exit.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{\code{<Disease>_history}}{1 if prevalent case (from prevalent_sources), 0 otherwise}
#'     \item{\code{<Disease>_incident}}{1 if incident case (from outcome_sources), 0 otherwise}
#'     \item{outcome_status}{For one primary disease, event indicator
#'       (1=event, 0=censored, NA=prevalent case).}
#'     \item{outcome_surv_time}{For one primary disease, follow-up time in
#'       years (NA for prevalent cases).}
#'     \item{\code{<Disease>_status}}{For multiple primary diseases, an
#'       independent event indicator for each disease.}
#'     \item{\code{<Disease>_surv_time}}{For multiple primary diseases,
#'       independent follow-up time in years for each disease.}
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
#'   \item OPCS4 procedure dates are often useful for procedure-defined endpoints or
#'     surgical history, but may occur later than the true biological disease onset.
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
#' Follow-up time calculation (independently for every selected
#' \code{primary_disease}):
#' \itemize{
#'   \item Prevalent case (primary disease): NA (not at risk)
#'   \item Incident case: (diagnosis_date - baseline_date) / 365.25
#'   \item Censored: (min(death_date, censor_date) - baseline_date) / 365.25
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
                                    time_skeleton = NULL,
                                    primary_disease = NULL,
                                    output = c("wide", "long"),
                                    include_all = TRUE,
                                    show_flow = TRUE,
                                    dt_threads = NULL) {

  output <- match.arg(output)

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  if (!is.null(dt_threads)) {
    if (!is.numeric(dt_threads) || length(dt_threads) != 1 || is.na(dt_threads) || dt_threads < 1) {
      stop("`dt_threads` must be NULL or a single positive integer", call. = FALSE)
    }
    dt_threads <- as.integer(dt_threads)
    old_threads <- data.table::getDTthreads()
    data.table::setDTthreads(threads = dt_threads)
    on.exit(data.table::setDTthreads(threads = old_threads), add = TRUE)
    message(sprintf(
      "[build_survival_dataset] data.table threads: %d (will restore to %d on exit)",
      dt_threads, old_threads
    ))
  }

  if (is.null(time_skeleton) && !baseline_col %in% names(dt)) {
    stop(sprintf("Baseline column not found: %s", baseline_col))
  }

  skeleton_info <- .survival_prepare_time_skeleton(
    time_skeleton = time_skeleton,
    dt = dt,
    baseline_col = baseline_col
  )
  if (!is.null(skeleton_info)) {
    dt <- skeleton_info$dt
    baseline_col <- skeleton_info$baseline_col
  }

  diseases <- names(disease_definitions)
  if (length(diseases) == 0L) {
    stop("No disease definitions provided", call. = FALSE)
  }

  if (is.null(primary_disease)) {
    primary_disease <- diseases[[1L]]
  }
  if (
    !is.character(primary_disease) ||
      length(primary_disease) == 0L ||
      anyNA(primary_disease) ||
      any(!nzchar(primary_disease))
  ) {
    stop("`primary_disease` must be NULL or a non-empty character vector.", call. = FALSE)
  }
  primary_disease <- unique(primary_disease)
  missing_primary <- setdiff(primary_disease, diseases)
  if (length(missing_primary) > 0L) {
    stop(
      sprintf(
        "Primary disease(s) not found in disease definitions: %s",
        paste(missing_primary, collapse = ", ")
      ),
      call. = FALSE
    )
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

  message("[build_survival_dataset] Calculating survival times...")

  if (!is.null(skeleton_info)) {
    all_eids <- data.table::copy(skeleton_info$time_skeleton)
    all_eids <- all_eids[, .(
      eid,
      baseline_date,
      death_date,
      end_date = followup_end_date,
      default_surv_time = followup_time_years,
      time_followup_end_reason = followup_end_reason,
      time_valid_followup = valid_followup
    )]
  } else {
    all_eids <- .ukb_followup_window(
      dt = dt,
      baseline_col = baseline_col,
      censor_date = censor_date
    )
    all_eids[, end_date := followup_end_date]
    all_eids[, default_surv_time := as.numeric(end_date - baseline_date) / 365.25]
    all_eids[, `:=`(
      time_followup_end_reason = data.table::fcase(
        !is.na(death_date) & end_date == death_date,
        "death",
        !is.na(lost_to_followup_date) & end_date == lost_to_followup_date,
        "lost_to_followup",
        default = "administrative_censoring"
      ),
      time_valid_followup = !is.na(baseline_date) & !is.na(end_date) & default_surv_time >= 0
    )]
    all_eids[, c("lost_to_followup_date", "followup_end_date") := NULL]
  }

  if (output == "long") {
    if (isTRUE(show_flow)) {
      message("[build_survival_dataset] `show_flow` is currently available for output='wide'; skipping flow print for long output.")
    }

    full_list <- lapply(diseases, function(d) {
      # Use prevalent_sources for history, outcome_sources for incident
      d_prevalent <- prevalent_cases_dt[disease == d]
      d_outcome <- outcome_cases_dt[disease == d]

      cohort <- data.table::copy(all_eids)
      cohort[, `:=`(
        disease = d,
        status = data.table::fifelse(
          get("time_valid_followup"),
          0L,
          NA_integer_
        ),
        prevalent_case = FALSE,
        earliest_date = as.Date(NA),
        diagnosis_source = NA_character_,
        surv_time = data.table::fifelse(
          get("time_valid_followup"),
          get("default_surv_time"),
          NA_real_
        )
      )]

      # Mark prevalent cases (from prevalent_sources including self-report)
      if (nrow(d_prevalent) > 0) {
        prevalent_records <- d_prevalent[prevalent_case == TRUE]
        prevalent_eids <- prevalent_records$eid
        if (nrow(prevalent_records) > 0L) {
          cohort[prevalent_records, `:=`(
            prevalent_case = TRUE,
            status = NA_integer_,
            surv_time = NA_real_,
            earliest_date = i.earliest_date,
            diagnosis_source = i.diagnosis_source
          ), on = "eid"]
        }
      } else {
        prevalent_eids <- cohort$eid[0]
      }

      # Mark incident cases and update survival time (from outcome_sources)
      if (nrow(d_outcome) > 0) {
        non_prevalent_outcome <- d_outcome[!eid %in% prevalent_eids]
        cohort[non_prevalent_outcome,
          `:=`(
            status = i.status,
            earliest_date = i.earliest_date,
            diagnosis_source = i.diagnosis_source,
            surv_time = i.surv_time
          ),
          on = "eid"
        ]
      }

      cohort[, c("baseline_date", "death_date", "end_date", "default_surv_time",
                  "time_followup_end_reason", "time_valid_followup") := NULL]
      cohort
    })

    full_cohort <- data.table::rbindlist(full_list, use.names = TRUE, fill = TRUE)
    data.table::setorder(full_cohort, disease, eid)

    message("[build_survival_dataset] Complete")
    output_columns <- c(
      "eid", "disease", "status", "surv_time", "prevalent_case",
      "earliest_date", "diagnosis_source"
    )
    if (!isTRUE(include_all)) {
      full_cohort <- full_cohort[status == 1L & prevalent_case == FALSE]
    }
    return(full_cohort[, output_columns, with = FALSE])
  }

  # Generate wide format: history from prevalent_sources, incident from outcome_sources
  wide_dt <- generate_wide_format_dual_source(
    dt,
    disease_definitions,
    prevalent_sources = prevalent_sources,
    outcome_sources = outcome_sources,
    censor_date = censor_date,
    baseline_col = baseline_col,
    prevalent_long = prevalent_cases_dt,
    outcome_long = outcome_cases_dt
  )

  endpoint_columns <- .survival_endpoint_columns(primary_disease)
  outcome_dt <- .survival_build_wide_endpoints(
    all_eids = all_eids,
    prevalent_cases_dt = prevalent_cases_dt,
    outcome_cases_dt = outcome_cases_dt,
    primary_diseases = primary_disease,
    endpoint_columns = endpoint_columns
  )

  result_dt <- data.table::merge.data.table(
    data.table::copy(dt),
    wide_dt,
    by = "eid",
    all.x = TRUE
  )
  result_dt <- data.table::merge.data.table(result_dt, outcome_dt, by = "eid", all.x = TRUE)
  internal_time_cols <- grep("^\\.ukba_time_skeleton_", names(result_dt), value = TRUE)
  if (length(internal_time_cols) > 0L) {
    result_dt[, (internal_time_cols) := NULL]
  }
  data.table::setorder(result_dt, eid)

  if (isTRUE(show_flow)) {
    flow_list <- lapply(seq_along(primary_disease), function(index) {
      disease <- primary_disease[[index]]
      columns <- endpoint_columns[[index]]
      flow_dt <- .survival_participant_flow(
        result_dt = result_dt,
        raw_n = nrow(dt),
        disease = disease,
        status_col = columns[["status"]],
        time_col = columns[["time"]]
      )

      flow_print <- data.table::copy(flow_dt)
      flow_print$retained_from_prev <- ifelse(
        is.na(flow_print$retained_from_prev),
        NA_character_,
        sprintf("%.2f%%", flow_print$retained_from_prev * 100)
      )
      flow_print$retained_from_raw <- ifelse(
        is.na(flow_print$retained_from_raw),
        NA_character_,
        sprintf("%.2f%%", flow_print$retained_from_raw * 100)
      )

      message(sprintf(
        "[build_survival_dataset] Participant flow (primary disease: %s)",
        disease
      ))
      print(flow_print)
      flow_dt
    })

    if (length(flow_list) == 1L) {
      attr(result_dt, "participant_flow") <- flow_list[[1L]]
    } else {
      named_flow <- lapply(seq_along(flow_list), function(index) {
        flow <- data.table::copy(flow_list[[index]])
        flow[, disease := primary_disease[[index]]]
        data.table::setcolorder(flow, c("disease", setdiff(names(flow), "disease")))
        flow
      })
      attr(result_dt, "participant_flow") <- data.table::rbindlist(
        named_flow,
        use.names = TRUE
      )
    }
  }

  message("[build_survival_dataset] Complete")

  return(result_dt)
}

.survival_endpoint_columns <- function(primary_diseases) {
  multiple <- length(primary_diseases) > 1L
  lapply(primary_diseases, function(disease) {
    if (multiple) {
      c(
        status = paste0(disease, "_status"),
        time = paste0(disease, "_surv_time")
      )
    } else {
      c(status = "outcome_status", time = "outcome_surv_time")
    }
  })
}

.survival_build_wide_endpoints <- function(all_eids,
                                           prevalent_cases_dt,
                                           outcome_cases_dt,
                                           primary_diseases,
                                           endpoint_columns) {
  outcome_dt <- data.table::copy(all_eids)

  for (index in seq_along(primary_diseases)) {
    disease_key <- primary_diseases[[index]]
    status_col <- endpoint_columns[[index]][["status"]]
    time_col <- endpoint_columns[[index]][["time"]]

    data.table::set(outcome_dt, j = status_col, value = 0L)
    data.table::set(
      outcome_dt,
      j = time_col,
      value = as.numeric(outcome_dt[["default_surv_time"]])
    )

    prevalent_records <- prevalent_cases_dt[
      disease == disease_key & prevalent_case == TRUE
    ]
    prevalent_eids <- prevalent_records[["eid"]]

    disease_outcomes <- outcome_cases_dt[disease == disease_key]
    if (length(prevalent_eids) > 0L && nrow(disease_outcomes) > 0L) {
      disease_outcomes <- disease_outcomes[!eid %in% prevalent_eids]
    }

    if (nrow(disease_outcomes) > 0L) {
      outcome_index <- match(disease_outcomes[["eid"]], outcome_dt[["eid"]])
      keep <- !is.na(outcome_index)
      if (any(keep)) {
        data.table::set(
          outcome_dt,
          i = outcome_index[keep],
          j = status_col,
          value = as.integer(disease_outcomes[["status"]][keep])
        )
        data.table::set(
          outcome_dt,
          i = outcome_index[keep],
          j = time_col,
          value = as.numeric(disease_outcomes[["surv_time"]][keep])
        )
      }
    }

    if (length(prevalent_eids) > 0L) {
      prevalent_index <- match(prevalent_eids, outcome_dt[["eid"]])
      prevalent_index <- prevalent_index[!is.na(prevalent_index)]
      if (length(prevalent_index) > 0L) {
        data.table::set(
          outcome_dt,
          i = prevalent_index,
          j = status_col,
          value = NA_integer_
        )
        data.table::set(
          outcome_dt,
          i = prevalent_index,
          j = time_col,
          value = NA_real_
        )
      }
    }
  }

  outcome_dt[, c("baseline_date", "death_date", "end_date", "default_surv_time") := NULL]
  outcome_dt
}

.survival_participant_flow <- function(result_dt,
                                       raw_n,
                                       disease,
                                       status_col,
                                       time_col) {
  history_col <- paste0(disease, "_history")
  status_vec <- result_dt[[status_col]]
  time_vec <- result_dt[[time_col]]
  idx_non_missing_status <- !is.na(status_vec)

  history_rule <- sprintf("Keep %s == 0", history_col)
  idx_non_prevalent <- idx_non_missing_status
  if (history_col %in% names(result_dt)) {
    history_vec <- result_dt[[history_col]]
    idx_non_prevalent <- idx_non_missing_status &
      !is.na(history_vec) &
      history_vec == 0L
  } else {
    history_rule <- sprintf("Column %s not found; skip this filter", history_col)
  }

  idx_valid_time <- idx_non_prevalent &
    !is.na(time_vec) &
    time_vec >= 0

  result_n <- nrow(result_dt)
  non_missing_status_n <- sum(idx_non_missing_status)
  non_prevalent_n <- sum(idx_non_prevalent)
  valid_time_n <- sum(idx_valid_time)

  flow_dt <- data.table::data.table(
    step = c(
      "Raw cohort",
      "After build_survival_dataset",
      sprintf("Keep non-missing %s", status_col),
      sprintf("Exclude baseline prevalent %s", disease),
      sprintf("Keep valid %s", time_col)
    ),
    n_before = as.integer(c(
      raw_n,
      raw_n,
      result_n,
      non_missing_status_n,
      non_prevalent_n
    )),
    n_after = as.integer(c(
      raw_n,
      result_n,
      non_missing_status_n,
      non_prevalent_n,
      valid_time_n
    )),
    exclusion_rule = c(
      "None",
      "Function output row count check",
      sprintf("Exclude %s is NA", status_col),
      history_rule,
      sprintf("Keep !is.na(%s) & %s >= 0", time_col, time_col)
    )
  )

  data.table::set(
    flow_dt,
    j = "excluded",
    value = flow_dt[["n_before"]] - flow_dt[["n_after"]]
  )
  data.table::set(
    flow_dt,
    j = "retained_from_prev",
    value = data.table::fifelse(
      flow_dt[["n_before"]] > 0,
      flow_dt[["n_after"]] / flow_dt[["n_before"]],
      NA_real_
    )
  )
  data.table::set(
    flow_dt,
    j = "retained_from_raw",
    value = if (raw_n > 0L) {
      flow_dt[["n_after"]] / raw_n
    } else {
      rep(NA_real_, nrow(flow_dt))
    }
  )
  flow_dt
}

.survival_prepare_time_skeleton <- function(time_skeleton,
                                            dt,
                                            baseline_col = "p53_i0") {
  if (is.null(time_skeleton)) {
    return(NULL)
  }
  if (!data.table::is.data.table(time_skeleton)) {
    time_skeleton <- data.table::as.data.table(time_skeleton)
  } else {
    time_skeleton <- data.table::copy(time_skeleton)
  }

  required <- c("eid", "baseline_date", "followup_end_date", "followup_time_years")
  missing <- setdiff(required, names(time_skeleton))
  if (length(missing) > 0L) {
    stop(
      "`time_skeleton` must contain columns: ",
      paste(required, collapse = ", "),
      ". Missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(time_skeleton$eid) > 0L) {
    stop("`time_skeleton$eid` must be unique.", call. = FALSE)
  }
  if (!"death_date" %in% names(time_skeleton)) {
    time_skeleton[, death_date := as.Date(NA)]
  }
  if (!"followup_end_reason" %in% names(time_skeleton)) {
    time_skeleton[, followup_end_reason := NA_character_]
  }
  if (!"valid_followup" %in% names(time_skeleton)) {
    time_skeleton[, valid_followup := !is.na(baseline_date) &
                    !is.na(followup_end_date) &
                    !is.na(followup_time_years) &
                    followup_time_years >= 0]
  }

  time_skeleton[, baseline_date := .safe_as_date(baseline_date, col_name = "time_skeleton$baseline_date")]
  time_skeleton[, followup_end_date := .safe_as_date(followup_end_date, col_name = "time_skeleton$followup_end_date")]
  time_skeleton[, death_date := .safe_as_date(death_date, col_name = "time_skeleton$death_date")]
  time_skeleton[, followup_time_years := suppressWarnings(as.numeric(followup_time_years))]

  missing_eids <- setdiff(dt$eid, time_skeleton$eid)
  if (length(missing_eids) > 0L) {
    stop(
      "`time_skeleton` is missing ",
      length(missing_eids),
      " participant(s) from `dt`.",
      call. = FALSE
    )
  }

  skeleton_dt <- time_skeleton[eid %in% dt$eid]
  temp_baseline_col <- ".ukba_time_skeleton_baseline_date"
  while (temp_baseline_col %in% names(dt)) {
    temp_baseline_col <- paste0(temp_baseline_col, "_")
  }

  temp_followup_col <- ".ukba_time_skeleton_followup_end_date"
  while (temp_followup_col %in% names(dt)) {
    temp_followup_col <- paste0(temp_followup_col, "_")
  }

  merged_dt <- data.table::merge.data.table(
    data.table::copy(dt),
    skeleton_dt[, .(
      eid,
      .ukba_tmp_baseline = baseline_date,
      .ukba_tmp_followup_end = followup_end_date
    )],
    by = "eid",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setnames(merged_dt, ".ukba_tmp_baseline", temp_baseline_col)
  data.table::setnames(merged_dt, ".ukba_tmp_followup_end", temp_followup_col)
  data.table::setorder(merged_dt, eid)
  data.table::setorder(skeleton_dt, eid)

  list(
    dt = merged_dt,
    baseline_col = temp_baseline_col,
    time_skeleton = skeleton_dt
  )
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
                               exclude_prevalent = TRUE,
                               dt_threads = NULL) {

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
    include_all = TRUE,
    dt_threads = dt_threads
  )

  if (exclude_prevalent) {
    full_cohort <- full_cohort[prevalent_case == FALSE | is.na(prevalent_case)]
  }

  full_cohort <- full_cohort[!is.na(surv_time) & surv_time > 0]
  data.table::setorder(full_cohort, disease, eid)

  return(full_cohort)
}
