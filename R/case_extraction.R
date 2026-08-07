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
#'   Valid options: "ICD10", "ICD9", "Self-report", "Death", "OPCS4",
#'   "CancerRegistry", "FirstOccurrence", "Algorithm".
#'   "OPCS4" uses hospital inpatient summary operations
#'   (\code{p41272} + \code{p41282_a*}) and requires
#'   \code{opcs4_pattern} in the disease definition.
#'   "Algorithm" uses UK Biobank algorithmically-defined outcomes (Category 42)
#'   which combine multiple data sources with high positive predictive value.
#'   Requires \code{algo_date_field} in the disease definition.
#'   If \code{algo_source_field} is also provided, output
#'   \code{diagnosis_source} is refined as \code{"Algorithm_<source_code>"}.
#'   "FirstOccurrence" uses UK Biobank First Occurrence date fields
#'   (Category 1712, \code{p13xxxx}) and requires
#'   \code{first_occurrence_fields} in the disease definition.
#'   "CancerRegistry" uses UK Biobank cancer register records
#'   (\code{p40006_i*} + \code{p40005_i*}) and requires
#'   \code{cancer_icd10_pattern} in the disease definition.
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
#'   \item Procedure-augmented definitions for surgical phenotypes using OPCS4
#'   \item Cancer registry ascertainment for malignant neoplasm endpoints
#'   \item First Occurrence fields for UKB's pre-mapped 3-character ICD-10 outcomes
#'   \item Source-specific case counts for methods reporting
#'   \item UK Biobank algorithmically-defined outcomes for validated case ascertainment
#' }
#'
#' The "Algorithm" source reads date fields from UK Biobank Category 42
#' (Algorithmically-defined outcomes). These are pre-computed by the UK Biobank
#' outcome adjudication group, combining self-report, hospital admissions,
#' and death records with high positive predictive value.
#' Records with date \code{1900-01-01} are excluded (date unknown).
#' If a source field is available in the definition, it is propagated into
#' \code{diagnosis_source} as \code{"Algorithm_<source_code>"}.
#'
#' The "FirstOccurrence" source reads singular UKB fields such as
#' \code{p131298_i0} or \code{p131298} for I21 first reported. Values with
#' UKB special date coding 819 (\code{1900-01-01}, \code{1901-01-01},
#' \code{1902-02-02}, \code{1903-03-03}, \code{1909-09-09}, and
#' \code{2037-07-07}) are excluded.
#'
#' @export
extract_cases_by_source <- function(dt,
                                     disease_definitions,
                                     sources = c("ICD10", "ICD9", "Self-report", "Death"),
                                     censor_date = as.Date("2023-10-31"),
                                     baseline_col = "p53_i0") {

  valid_sources <- c("ICD10", "ICD9", "Self-report", "Death", "OPCS4", "CancerRegistry", "FirstOccurrence", "Algorithm")
  sources <- match.arg(sources, valid_sources, several.ok = TRUE)

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  if (!baseline_col %in% names(dt)) {
    stop(sprintf("Baseline column not found: %s", baseline_col))
  }

  message(sprintf("[extract_cases_by_source] Using sources: %s", paste(sources, collapse = ", ")))

  # Extract only requested sources
  icd10_long <- if ("ICD10" %in% sources) parse_icd10_diagnoses(dt) else data.table()
  icd9_long <- if ("ICD9" %in% sources) parse_icd9_diagnoses(dt) else data.table()
  sr_long <- if ("Self-report" %in% sources) parse_self_reported_illnesses(dt, baseline_col) else data.table()
  death_long <- if ("Death" %in% sources) parse_death_records(dt) else data.table()
  opcs4_long <- if ("OPCS4" %in% sources) parse_opcs4_procedures(dt) else data.table()
  cancer_long <- if ("CancerRegistry" %in% sources) parse_cancer_registry(dt) else data.table()

  followup_window <- .ukb_followup_window(
    dt = dt,
    baseline_col = baseline_col,
    censor_date = censor_date
  )

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

    death_pattern <- if (!is.null(def$death_icd10)) def$death_icd10 else def$icd10_pattern
    if ("Death" %in% sources && !is.null(death_pattern) && nrow(death_long) > 0) {
      filtered <- filter_death_codes(death_long, death_pattern, disease_key)
      if (nrow(filtered) > 0) diagnosis_sources$death <- aggregate_death_as_diagnosis(filtered)
    }

    if ("OPCS4" %in% sources && !is.null(def$opcs4_pattern) && nrow(opcs4_long) > 0) {
      filtered <- filter_opcs4_codes(opcs4_long, def$opcs4_pattern, disease_key)
      if (nrow(filtered) > 0) diagnosis_sources$opcs4 <- aggregate_opcs4_earliest(filtered)
    }

    if ("CancerRegistry" %in% sources && !is.null(def$cancer_icd10_pattern) && nrow(cancer_long) > 0) {
      filtered <- filter_cancer_registry(
        cancer_long,
        pattern = def$cancer_icd10_pattern,
        disease_label = disease_key,
        histology = def$cancer_histology,
        behaviour = def$cancer_behaviour
      )
      if (nrow(filtered) > 0) {
        diagnosis_sources$cancer <- aggregate_cancer_registry_earliest(filtered)
      }
    }

    if ("FirstOccurrence" %in% sources && !is.null(def$first_occurrence_fields)) {
      fo_long <- .parse_first_occurrence_records(
        dt = dt,
        fields = def$first_occurrence_fields,
        source_fields = def$first_occurrence_source_fields
      )
      if (nrow(fo_long) > 0) {
        diagnosis_sources$first_occurrence <- .aggregate_first_occurrence_earliest(fo_long, disease_key)
      }
    }

    # Algorithm source: UK Biobank algorithmically-defined outcomes (Category 42)
    if ("Algorithm" %in% sources && !is.null(def$algo_date_field)) {
      algo_date_fields <- as.integer(def$algo_date_field)
      algo_source_fields <- def$algo_source_field
      if (is.null(algo_source_fields)) {
        algo_source_fields <- rep(NA_integer_, length(algo_date_fields))
      }
      if (length(algo_source_fields) != length(algo_date_fields)) {
        stop(
          "Algorithm date/source field lengths differ for disease '",
          disease_key,
          "'. Recreate the definition with create_disease_definition().",
          call. = FALSE
        )
      }

      algorithm_records <- vector("list", length(algo_date_fields))
      for (algorithm_index in seq_along(algo_date_fields)) {
        date_field <- algo_date_fields[[algorithm_index]]
        source_field <- algo_source_fields[[algorithm_index]]
        algo_col_candidates <- c(
          paste0("p", date_field, "_i0"),
          paste0("p", date_field)
        )
        algo_col <- algo_col_candidates[algo_col_candidates %in% names(dt)][1]
        if (is.na(algo_col)) {
          message(sprintf(
            "  [Algorithm] No date column found for %s field %s, skipping that field",
            disease_key,
            date_field
          ))
          next
        }

        algo_source_col <- NA_character_
        if (!is.na(source_field)) {
          algo_source_candidates <- c(
            paste0("p", source_field, "_i0"),
            paste0("p", source_field)
          )
          algo_source_col <- algo_source_candidates[
            algo_source_candidates %in% names(dt)
          ][1]
        }
        has_algo_source <- !is.na(algo_source_col)

        if (has_algo_source) {
          algo_dt <- dt[, .(
            eid,
            algo_date = .safe_as_date(get(algo_col), col_name = algo_col),
            algo_source = as.character(get(algo_source_col))
          )]
        } else {
          algo_dt <- dt[, .(
            eid,
            algo_date = .safe_as_date(get(algo_col), col_name = algo_col),
            algo_source = NA_character_
          )]
        }

        algo_dt <- algo_dt[
          !is.na(algo_date) & algo_date != as.Date("1900-01-01")
        ]
        if (nrow(algo_dt) > 0L) {
          algo_dt[, source := fifelse(
            !is.na(algo_source) & trimws(algo_source) != "",
            paste0("Algorithm_", trimws(algo_source)),
            "Algorithm"
          )]
          algo_dt[, `:=`(
            disease = disease_key,
            earliest_date = algo_date,
            diagnosis_date_quality = "recorded"
          )]
          algo_dt[, c("algo_date", "algo_source") := NULL]
          algorithm_records[[algorithm_index]] <- algo_dt
        }
      }
      algorithm_records <- algorithm_records[
        !vapply(algorithm_records, is.null, logical(1))
      ]
      if (length(algorithm_records) > 0L) {
        diagnosis_sources$algo <- rbindlist(
          algorithm_records,
          use.names = TRUE,
          fill = TRUE
        )
      }
    }

    if (length(diagnosis_sources) == 0) return(NULL)

    all_diagnoses <- rbindlist(diagnosis_sources, use.names = TRUE, fill = TRUE)

    if (!"diagnosis_date_quality" %in% names(all_diagnoses)) {
      all_diagnoses[, "diagnosis_date_quality" := NA_character_]
    }
    earliest_per_person <- all_diagnoses[
      ,
      {
        dated <- which(!is.na(earliest_date))
        selected <- if (length(dated) > 0L) {
          dated[[which.min(earliest_date[dated])]]
        } else {
          1L
        }
        list(
          earliest_date = earliest_date[[selected]],
          diagnosis_source = source[[selected]],
          diagnosis_date_quality = .SD[["diagnosis_date_quality"]][[selected]]
        )
      },
      by = .(eid, disease)
    ]

    return(earliest_per_person)
  })

  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) {
    warning("[extract_cases_by_source] No cases found")
    return(data.table(
      eid = integer(0), disease = character(0), earliest_date = as.Date(character(0)),
      diagnosis_source = character(0), prevalent_case = logical(0),
      status = integer(0), surv_time = numeric(0)
    ))
  }

  diagnosis_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

  # Calculate survival metrics
  surv_dt <- merge.data.table(
    diagnosis_dt,
    followup_window,
    by = "eid",
    all.x = TRUE,
    sort = FALSE
  )

  surv_dt[, "undated_diagnosis" :=
    is.na(get("earliest_date")) & !is.na(get("diagnosis_source"))]
  surv_dt[, prevalent_case := fifelse(
    get("undated_diagnosis"),
    NA,
    !is.na(earliest_date) & earliest_date <= baseline_date
  )]
  surv_dt[, status := as.integer(
    !is.na(earliest_date) &
      !is.na(baseline_date) &
      !is.na(followup_end_date) &
      earliest_date > baseline_date &
      earliest_date <= followup_end_date
  )]
  surv_dt[is.na(baseline_date) | is.na(followup_end_date), status := NA_integer_]
  surv_dt[get("undated_diagnosis") %in% TRUE, status := NA_integer_]

  surv_dt[, end_date := fifelse(
    status == 1L, earliest_date,
    followup_end_date
  )]
  surv_dt[is.na(end_date), end_date := followup_end_date]
  surv_dt[, surv_time := as.numeric(end_date - baseline_date) / 365.25]
  surv_dt[surv_time < 0, `:=`(surv_time = NA_real_, status = NA_integer_)]

  surv_dt[, c(
    "baseline_date", "death_date", "lost_to_followup_date",
    "followup_end_date", "end_date", "undated_diagnosis"
  ) := NULL]
  setorder(surv_dt, disease, eid)

  return(surv_dt)
}


#' Extract participant-level disease diagnosis status
#'
#' @description
#' Defines whether each participant has a selected disease using one or more
#' UK Biobank evidence sources. This is the recommended public helper when the
#' goal is disease ascertainment rather than construction of a full survival
#' cohort. For survival-ready endpoints, use \code{\link{build_survival_dataset}}.
#'
#' @param dt A data.table or data.frame containing UKB data.
#' @param disease Character vector of disease keys or disease names.
#' @param disease_definitions Optional named list of disease definitions. If
#'   \code{NULL}, \code{\link{get_predefined_diseases}} is used.
#' @param sources Character vector of evidence sources. Valid options are
#'   \code{"ICD10"}, \code{"ICD9"}, \code{"Self-report"}, \code{"Death"},
#'   \code{"OPCS4"}, \code{"CancerRegistry"}, \code{"FirstOccurrence"}, and
#'   \code{"Algorithm"}.
#' @param censor_date Administrative censoring date.
#' @param baseline_col Column name for baseline assessment date.
#' @param include_all Logical. If \code{TRUE}, return one row per participant
#'   per disease, including non-cases. If \code{FALSE}, return diagnosed
#'   participants only.
#'
#' @return A data.table with participant-level diagnosis status, first diagnosis
#'   date, diagnosis source, prevalent and incident indicators, and survival
#'   fields returned by \code{\link{extract_cases_by_source}} where available.
#' @export
#'
extract_disease_diagnosis <- function(dt,
                                      disease,
                                      disease_definitions = NULL,
                                      sources = c("ICD10", "ICD9", "Self-report", "Death"),
                                      censor_date = as.Date("2023-10-31"),
                                      baseline_col = "p53_i0",
                                      include_all = TRUE) {
  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }
  if (!"eid" %in% names(dt)) {
    stop("Column 'eid' was not found in `dt`.", call. = FALSE)
  }
  if (!is.character(disease) || length(disease) == 0L || anyNA(disease)) {
    stop("`disease` must be a non-empty character vector.", call. = FALSE)
  }
  if (is.null(disease_definitions)) {
    disease_definitions <- get_predefined_diseases()
  }
  if (!is.list(disease_definitions) || is.null(names(disease_definitions))) {
    stop("`disease_definitions` must be a named list.", call. = FALSE)
  }

  disease_keys <- .resolve_disease_keys(disease, disease_definitions)
  selected_defs <- disease_definitions[disease_keys]

  cases <- extract_cases_by_source(
    dt = dt,
    disease_definitions = selected_defs,
    sources = sources,
    censor_date = censor_date,
    baseline_col = baseline_col
  )
  if (!is.data.table(cases)) {
    cases <- as.data.table(cases)
  }

  all_pairs <- CJ(
    eid = dt[["eid"]],
    disease = disease_keys,
    unique = TRUE
  )

  if (nrow(cases) > 0L) {
    out <- merge.data.table(
      all_pairs,
      cases,
      by = c("eid", "disease"),
      all.x = TRUE
    )
  } else {
    out <- all_pairs
    out[, `:=`(
      earliest_date = as.Date(NA),
      diagnosis_source = NA_character_,
      prevalent_case = NA,
      status = NA_integer_,
      surv_time = NA_real_
    )]
  }

  out[, diagnosed := !is.na(earliest_date) |
        diagnosis_source %in% "Self-report_unknown_date"]
  out[, prevalent_case := fifelse(is.na(prevalent_case), FALSE, prevalent_case)]
  out[, incident_case := fifelse(is.na(status), 0L, as.integer(status == 1L))]
  out[, status := fifelse(is.na(status), 0L, status)]

  setcolorder(
    out,
    c(
      "eid", "disease", "diagnosed", "prevalent_case", "incident_case",
      "earliest_date", "diagnosis_source", "status", "surv_time"
    )
  )
  setorder(out, disease, eid)
  if (!isTRUE(include_all)) {
    out <- out[diagnosed == TRUE]
  }
  out[]
}

.resolve_disease_keys <- function(disease, disease_definitions) {
  keys <- names(disease_definitions)
  out <- character(length(disease))
  def_names <- vapply(disease_definitions, function(x) {
    if (is.list(x) && !is.null(x$name)) as.character(x$name)[1] else NA_character_
  }, character(1))

  for (i in seq_along(disease)) {
    x <- disease[[i]]
    if (x %in% keys) {
      out[[i]] <- x
      next
    }
    hit <- keys[tolower(keys) == tolower(x)]
    if (length(hit) == 1L) {
      out[[i]] <- hit
      next
    }
    hit <- keys[tolower(def_names) == tolower(x)]
    if (length(hit) == 1L) {
      out[[i]] <- hit
      next
    }
    stop(
      "Disease not found in `disease_definitions`: ", x,
      ". Use names(get_predefined_diseases()) to inspect supported keys.",
      call. = FALSE
    )
  }
  unique(out)
}


#' @title Generate Wide-Format with Dual Source Definition
#'
#' @description
#' Internal function that generates wide-format disease status using separate
#' sources for prevalent (history) and incident cases. This supports the common
#' epidemiological practice of using self-report for baseline exclusion but not
#' for outcome ascertainment.
#'
#' @param dt A data.table containing UKB data.
#' @param disease_definitions Named list of disease definitions.
#' @param prevalent_sources Sources for identifying prevalent cases.
#' @param outcome_sources Sources for identifying incident cases.
#' @param censor_date Administrative censoring date.
#' @param baseline_col Column name for baseline date.
#'
#' @return A data.table with _history and _incident columns per disease.
#'
#' @keywords internal
generate_wide_format_dual_source <- function(dt,
                                              disease_definitions,
                                              prevalent_sources,
                                              outcome_sources,
                                              censor_date,
                                              baseline_col,
                                              prevalent_long = NULL,
                                              outcome_long = NULL) {

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  # Extract cases only if precomputed inputs are not provided.
  # This avoids duplicate expensive parsing in high-volume pipelines.
  if (is.null(prevalent_long)) {
    prevalent_long <- extract_cases_by_source(
      dt, disease_definitions, prevalent_sources, censor_date, baseline_col
    )
  }
  if (!is.data.table(prevalent_long)) {
    prevalent_long <- as.data.table(prevalent_long)
  }

  if (is.null(outcome_long)) {
    outcome_long <- extract_cases_by_source(
      dt, disease_definitions, outcome_sources, censor_date, baseline_col
    )
  }
  if (!is.data.table(outcome_long)) {
    outcome_long <- as.data.table(outcome_long)
  }

  all_eids <- dt[, .(eid)]
  wide_dt <- copy(all_eids)

  diseases <- names(disease_definitions)

  for (d in diseases) {
    # History from prevalent_sources
    d_prevalent <- prevalent_long[disease == d]
    # Incident from outcome_sources
    d_outcome <- outcome_long[disease == d]

    d_wide <- copy(all_eids)

    # Mark history (prevalent) from prevalent_sources
    prevalent_eids <- d_prevalent[prevalent_case == TRUE, eid]
    if (nrow(d_prevalent) > 0) {
      d_wide[, (paste0(d, "_history")) := as.integer(eid %in% prevalent_eids)]
    } else {
      d_wide[, (paste0(d, "_history")) := 0L]
    }

    # Mark incident from outcome_sources
    if (nrow(d_outcome) > 0) {
      incident_eids <- d_outcome[
        status == 1L & !eid %in% prevalent_eids,
        eid
      ]
      d_wide[, (paste0(d, "_incident")) := as.integer(eid %in% incident_eids)]
    } else {
      d_wide[, (paste0(d, "_incident")) := 0L]
    }

    # Replace NA with 0
    hist_col <- paste0(d, "_history")
    inc_col <- paste0(d, "_incident")
    set(d_wide, which(is.na(d_wide[[hist_col]])), hist_col, 0L)
    set(d_wide, which(is.na(d_wide[[inc_col]])), inc_col, 0L)

    d_wide <- d_wide[, c("eid", hist_col, inc_col), with = FALSE]
    wide_dt <- merge.data.table(wide_dt, d_wide, by = "eid", all.x = TRUE)
  }

  setorder(wide_dt, eid)
  return(wide_dt)
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
#' @noRd
generate_wide_format <- function(dt,
                                  disease_definitions,
                                  sources = c("ICD10", "ICD9", "Self-report", "Death"),
                                  censor_date = as.Date("2023-10-31"),
                                  baseline_col = "p53_i0",
                                  include_dates = FALSE) {

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  surv_long <- extract_cases_by_source(dt, disease_definitions, sources, censor_date, baseline_col)
  all_eids <- dt[, .(eid)]
  wide_dt <- copy(all_eids)

  diseases <- unique(surv_long$disease)

  for (d in diseases) {
    d_data <- surv_long[disease == d]

    d_wide <- merge.data.table(
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
    wide_dt <- merge.data.table(wide_dt, d_wide, by = "eid", all.x = TRUE)
  }

  # Replace NA with 0 for binary columns
  for (col in grep("_(history|incident)$", names(wide_dt), value = TRUE)) {
    set(wide_dt, which(is.na(wide_dt[[col]])), col, 0L)
  }

  setorder(wide_dt, eid)
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
#' @export
compare_data_sources <- function(dt,
                                  disease_definitions,
                                  baseline_col = "p53_i0") {

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  diseases <- names(disease_definitions)

  comparison_list <- lapply(diseases, function(d) {
    single_def <- disease_definitions[d]

    icd10_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = "ICD10", baseline_col = baseline_col),
      error = function(e) data.table()
    )

    icd9_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = "ICD9", baseline_col = baseline_col),
      error = function(e) data.table()
    )

    sr_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = "Self-report", baseline_col = baseline_col),
      error = function(e) data.table()
    )

    hospital_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = c("ICD10", "ICD9"), baseline_col = baseline_col),
      error = function(e) data.table()
    )

    all_cases <- tryCatch(
      extract_cases_by_source(dt, single_def, sources = c("ICD10", "ICD9", "Self-report"), baseline_col = baseline_col),
      error = function(e) data.table()
    )

    data.table(
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

  result <- rbindlist(comparison_list)
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
#' @noRd
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

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  # Generate wide format with all diseases
  wide_dt <- generate_wide_format(dt, disease_definitions, sources, censor_date, baseline_col)

  # Get survival data for primary outcome
  surv_long <- extract_cases_by_source(dt, disease_definitions, sources, censor_date, baseline_col)
  primary_data <- surv_long[disease == primary_outcome, .(eid, status, surv_time, prevalent_case)]

  # Get default survival time for non-cases
  death_dates <- get_death_dates(dt)
  baseline_dt <- dt[, .(eid, baseline_date = .safe_as_date(get(baseline_col), col_name = baseline_col))]
  control_info <- merge.data.table(baseline_dt, death_dates, by = "eid", all.x = TRUE)
  control_info[, end_date := pmin(death_date, censor_date, na.rm = TRUE)]
  control_info[is.na(end_date), end_date := censor_date]
  control_info[, control_surv_time := as.numeric(end_date - baseline_date) / 365.25]

  # Merge primary outcome data
  wide_dt <- merge.data.table(wide_dt, primary_data, by = "eid", all.x = TRUE)
  wide_dt <- merge.data.table(
    wide_dt, control_info[, .(eid, control_surv_time)], by = "eid", all.x = TRUE
  )

  # Set outcome columns
  wide_dt[, outcome_status := fifelse(is.na(status), 0L, status)]
  wide_dt[, outcome_surv_time := fifelse(is.na(surv_time), control_surv_time, surv_time)]
  wide_dt[, outcome_prevalent := fifelse(is.na(prevalent_case), FALSE, prevalent_case)]

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
  setorder(wide_dt, eid)

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
#'   Default: "ICD10". Options: "ICD10", "ICD9", "Self-report", "Death",
#'   "OPCS4", "CancerRegistry", "FirstOccurrence", "Algorithm".
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

  valid_sources <- c("ICD10", "ICD9", "Self-report", "Death", "OPCS4", "CancerRegistry", "FirstOccurrence", "Algorithm")
  sources <- match.arg(sources, valid_sources, several.ok = TRUE)

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
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
  result_dt <- copy(all_eids)

  # Create history column for each disease
  for (d in diseases) {
    d_data <- surv_long[disease == d & prevalent_case == TRUE, .(eid)]
    d_data[, (paste0(d, "_history")) := 1L]

    result_dt <- merge.data.table(
      result_dt,
      d_data,
      by = "eid",
      all.x = TRUE
    )

    # Fill NA with 0
    hist_col <- paste0(d, "_history")
    set(result_dt, which(is.na(result_dt[[hist_col]])), hist_col, 0L)
  }

  setorder(result_dt, eid)

  # Summary message
  for (d in diseases) {
    n_cases <- sum(result_dt[[paste0(d, "_history")]])
    message(sprintf("  %s_history: %d prevalent cases", d, n_cases))
  }

  return(result_dt)
}


#' @title Extract Baseline Diabetes Subtypes (T1DM/T2DM) with HbA1c Support
#'
#' @description
#' Extracts baseline prevalent Type 1 and Type 2 diabetes using existing
#' source-based disease history logic, and optionally augments Type 2
#' classification using baseline HbA1c.
#'
#' @param dt A data.table or data.frame containing UKB data.
#' @param disease_definitions Named list of disease definitions. If NULL,
#'   uses \code{\link{get_predefined_diseases}}.
#' @param sources Character vector specifying sources for baseline history.
#'   Options: "ICD10", "ICD9", "Self-report", "Death", "CancerRegistry",
#'   "FirstOccurrence", "Algorithm".
#' @param baseline_col Column name for baseline date. Default: \code{"p53_i0"}.
#' @param hba1c_col Column name for baseline HbA1c (mmol/mol).
#'   Default: \code{"p30750_i0"}.
#' @param hba1c_threshold Numeric threshold for diabetes by HbA1c.
#'   Default: \code{48} mmol/mol (equivalent to 6.5 percent).
#' @param include_hba1c Logical. If TRUE (default), HbA1c is used to
#'   augment T2DM classification.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{T1DM_history}{Baseline prevalent T1DM from selected sources (0/1)}
#'     \item{T2DM_history}{Baseline prevalent T2DM from selected sources (0/1)}
#'     \item{diabetes_hba1c}{Baseline HbA1c diabetes flag (0/1/NA)}
#'     \item{T2DM_history_enhanced}{T2DM from source history OR HbA1c criterion (0/1)}
#'     \item{Diabetes_history}{Any baseline diabetes (T1DM or enhanced T2DM) (0/1)}
#'     \item{diabetes_subtype}{"Type1", "Type2", or "No_diabetes"}
#'   }
#'
#' @details
#' This is a baseline classification helper and does not redefine incident
#' event logic. Type 1 has priority when both T1 and T2 evidence are present.
#'
#' @export
extract_diabetes_subtype_baseline <- function(dt,
                                              disease_definitions = NULL,
                                              sources = c("ICD10", "ICD9", "Self-report"),
                                              baseline_col = "p53_i0",
                                              hba1c_col = "p30750_i0",
                                              hba1c_threshold = 48,
                                              include_hba1c = TRUE) {

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  if (is.null(disease_definitions)) {
    disease_definitions <- get_predefined_diseases()
  }

  required_defs <- c("T1DM", "T2DM")
  missing_defs <- setdiff(required_defs, names(disease_definitions))
  if (length(missing_defs) > 0) {
    stop(sprintf(
      "Missing disease definition(s): %s. Please provide T1DM and T2DM in disease_definitions.",
      paste(missing_defs, collapse = ", ")
    ))
  }

  dm_hist <- extract_disease_history(
    dt = dt,
    diseases = required_defs,
    disease_definitions = disease_definitions,
    sources = sources,
    baseline_col = baseline_col
  )

  result_dt <- copy(dm_hist)

  if (include_hba1c) {
    if (!hba1c_col %in% names(dt)) {
      warning(sprintf("HbA1c column not found: %s. Falling back to source-only classification.", hba1c_col))
      result_dt[, diabetes_hba1c := NA_integer_]
    } else {
      hba1c_dt <- unique(dt[, .(eid, hba1c_value = suppressWarnings(as.numeric(get(hba1c_col))))], by = "eid")
      result_dt <- merge.data.table(result_dt, hba1c_dt, by = "eid", all.x = TRUE)
      result_dt[, diabetes_hba1c := fifelse(
        is.na(hba1c_value),
        NA_integer_,
        as.integer(hba1c_value >= hba1c_threshold)
      )]
      result_dt[, hba1c_value := NULL]
    }
  } else {
    result_dt[, diabetes_hba1c := NA_integer_]
  }

  result_dt[, T2DM_history_enhanced := {
    t2_code <- fifelse(is.na(T2DM_history), 0L, T2DM_history)
    hba1c_dm <- fifelse(is.na(diabetes_hba1c), 0L, diabetes_hba1c)
    as.integer((t2_code == 1L) | (hba1c_dm == 1L))
  }]

  result_dt[, Diabetes_history := {
    t1_code <- fifelse(is.na(T1DM_history), 0L, T1DM_history)
    as.integer((t1_code == 1L) | (T2DM_history_enhanced == 1L))
  }]

  result_dt[, diabetes_subtype := "No_diabetes"]
  result_dt[T2DM_history_enhanced == 1L, diabetes_subtype := "Type2"]
  result_dt[T1DM_history == 1L, diabetes_subtype := "Type1"]

  setorder(result_dt, eid)
  return(result_dt[])
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
#' @export
extract_disease_history_sensitivity <- function(dt,
                                                 diseases,
                                                 disease_definitions = NULL,
                                                 baseline_col = "p53_i0") {

  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
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
  setnames(hist_icd10, old_names, new_names)

  # Hospital (ICD-10 + ICD-9)
  hist_hospital <- extract_disease_history(
    dt, diseases, disease_definitions,
    sources = c("ICD10", "ICD9"), baseline_col = baseline_col
  )
  new_names <- paste0(diseases, "_history_hospital")
  setnames(hist_hospital, old_names, new_names)

  # All sources
  hist_all <- extract_disease_history(
    dt, diseases, disease_definitions,
    sources = c("ICD10", "ICD9", "Self-report"), baseline_col = baseline_col
  )
  new_names <- paste0(diseases, "_history_all")
  setnames(hist_all, old_names, new_names)

  # Merge all variants
  result_dt <- merge.data.table(hist_icd10, hist_hospital, by = "eid")
  result_dt <- merge.data.table(result_dt, hist_all, by = "eid")

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
