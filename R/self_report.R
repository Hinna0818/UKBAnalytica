#' @title Parse Self-Reported Illness Records
#'
#' @description
#' Extracts self-reported illness data from UK Biobank touchscreen questionnaire.
#' Converts coded illness data (p20002_i*_a*) and interpolated year of diagnosis
#' (p20008_i*_a*) into a standardized long-format data.table.
#'
#' @param dt A data.table or data.frame containing UKB data with columns:
#'   \code{eid}, \code{p20002_i*_a*}, and \code{p20008_i*_a*} columns.
#' @param baseline_col Column name for baseline date (default: "p53_i0").
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{sr_code}{Self-report illness code}
#'     \item{diag_date}{Approximate date of diagnosis}
#'     \item{source}{Data source identifier ("Self-report")}
#'     \item{instance}{Assessment instance (0, 1, 2, 3)}
#'     \item{array_idx}{Array index within instance}
#'     \item{baseline_report}{Whether the illness was reported at baseline}
#'     \item{diagnosis_date_quality}{How the diagnosis date was obtained}
#'   }
#'
#' @details
#' Year-to-date conversion logic:
#' \itemize{
#'   \item p20008 stores fractional years (e.g., 1983.5 = mid-1983)
#'   \item Fractional part * 12 = approximate month
#'   \item Special values (-1, -3) indicate "don't know" or "prefer not to answer"
#' }
#'
#' @import data.table
#' @export
parse_self_reported_illnesses <- function(dt, baseline_col = "p53_i0") {
  dt <- data.table::as.data.table(data.table::copy(dt))
  if (!"eid" %in% names(dt)) {
    stop("Column 'eid' was not found in 'dt'.", call. = FALSE)
  }
  if (anyNA(dt[["eid"]]) || anyDuplicated(dt[["eid"]])) {
    stop("Column 'eid' must be non-missing and unique.", call. = FALSE)
  }

  # Step 1: Melt p20002_i*_a* code columns
  code_cols <- grep("^p20002_i[0-9]+_a[0-9]+$", names(dt), value = TRUE)
  if (length(code_cols) == 0) {
    message("[parse_self_reported_illnesses] No p20002_i*_a* columns found")
    return(data.table::data.table(
      eid = integer(0), sr_code = integer(0),
      diag_date = as.Date(character(0)), source = character(0),
      instance = integer(0), array_idx = integer(0),
      baseline_report = logical(0), diagnosis_date_quality = character(0)
    ))
  }

  dt[, (code_cols) := lapply(.SD, as.integer), .SDcols = code_cols]

  codes_long <- data.table::melt(
    dt[, c("eid", code_cols), with = FALSE],
    id.vars = "eid", measure.vars = code_cols,
    variable.name = "col_name", value.name = "sr_code", na.rm = FALSE
  )

  # Extract instance and array index from column name
  codes_long[, c("instance", "array_idx") := {
    parts <- regmatches(col_name, regexec("p20002_i([0-9]+)_a([0-9]+)", col_name))
    list(
      as.integer(sapply(parts, function(x) if (length(x) == 3) x[2] else NA)),
      as.integer(sapply(parts, function(x) if (length(x) == 3) x[3] else NA))
    )
  }]
  codes_long[, col_name := NULL]

  # Remove invalid codes (NA or special values)
  codes_long <- codes_long[!is.na(sr_code) & sr_code > 0]

  if (nrow(codes_long) == 0) {
    message("[parse_self_reported_illnesses] No valid self-report codes found")
    return(data.table::data.table(
      eid = integer(0), sr_code = integer(0),
      diag_date = as.Date(character(0)), source = character(0),
      instance = integer(0), array_idx = integer(0),
      baseline_report = logical(0), diagnosis_date_quality = character(0)
    ))
  }

  # Step 2: Melt p20008_i*_a* year columns
  year_cols <- grep("^p20008_i[0-9]+_a[0-9]+$", names(dt), value = TRUE)
  if (length(year_cols) == 0) {
    warning("[parse_self_reported_illnesses] No p20008_i*_a* year columns found")
    result <- data.table::copy(codes_long)
    result[, diag_year := NA_real_]
  } else {
    eids_with_codes <- unique(codes_long$eid)
    dt_sub <- dt[eid %in% eids_with_codes, c("eid", year_cols), with = FALSE]
    dt_sub[, (year_cols) := lapply(.SD, as.numeric), .SDcols = year_cols]

    years_long <- data.table::melt(
      dt_sub, id.vars = "eid", measure.vars = year_cols,
      variable.name = "col_name", value.name = "diag_year", na.rm = FALSE
    )

    years_long[, c("instance", "array_idx") := {
      parts <- regmatches(col_name, regexec("p20008_i([0-9]+)_a([0-9]+)", col_name))
      list(
        as.integer(sapply(parts, function(x) if (length(x) == 3) x[2] else NA)),
        as.integer(sapply(parts, function(x) if (length(x) == 3) x[3] else NA))
      )
    }]
    years_long[, col_name := NULL]

    result <- data.table::merge.data.table(
      codes_long, years_long,
      by = c("eid", "instance", "array_idx"), all.x = TRUE
    )
  }

  # Step 4: Convert fractional year to date (robust to malformed values)
  result[, diag_year := suppressWarnings(as.numeric(diag_year))]
  result[!is.finite(diag_year) | diag_year < 0, diag_year := NA_real_]
  result[, diag_date := {
    year_int <- floor(diag_year)
    month_frac <- diag_year - year_int
    month_int <- pmax(1L, pmin(12L, as.integer(round(month_frac * 12)) + 1L))

    valid_idx <- !is.na(year_int) & is.finite(year_int) & year_int >= 1
    date_str <- rep(NA_character_, .N)
    if (any(valid_idx)) {
      date_str[valid_idx] <- sprintf(
        "%04d-%02d-01",
        as.integer(year_int[valid_idx]),
        month_int[valid_idx]
      )
    }

    .safe_as_date(date_str, col_name = "Self_report_diag_date")
  }]

  result[, diag_year := NULL]
  baseline_dt <- data.table::data.table(
    eid = dt[["eid"]],
    baseline_date = if (baseline_col %in% names(dt)) {
      .safe_as_date(dt[[baseline_col]], col_name = baseline_col)
    } else {
      as.Date(NA)
    }
  )
  result <- data.table::merge.data.table(
    result,
    baseline_dt,
    by = "eid",
    all.x = TRUE,
    sort = FALSE
  )
  result[, "baseline_report" := get("instance") == 0L]
  result[, "diagnosis_date_quality" := data.table::fcase(
    get("baseline_report") & is.na(get("diag_date")) &
      !is.na(get("baseline_date")),
    "baseline_date_imputed",
    get("baseline_report") & !is.na(get("diag_date")) &
      !is.na(get("baseline_date")) & get("diag_date") > get("baseline_date"),
    "baseline_date_clamped",
    get("baseline_report") & !is.na(get("diag_date")),
    "reported_before_or_at_baseline",
    !get("baseline_report") & !is.na(get("diag_date")),
    "estimated",
    default = "unknown"
  )]
  result[
    get("baseline_report") & !is.na(get("baseline_date")) &
      (is.na(get("diag_date")) | get("diag_date") > get("baseline_date")),
    "diag_date" := get("baseline_date")
  ]
  result[, `:=`(baseline_date = NULL, source = "Self-report")]
  data.table::setorder(result, eid, diag_date, na.last = TRUE)

  output_columns <- c(
    "eid", "sr_code", "diag_date", "source", "instance", "array_idx",
    "baseline_report", "diagnosis_date_quality"
  )
  return(result[, output_columns, with = FALSE])
}


#' @title Filter Self-Reported Illness Records by Code
#'
#' @description
#' Filters self-reported illness records by specific UKB illness codes.
#'
#' @param sr_long A data.table from \code{\link{parse_self_reported_illnesses}}.
#' @param codes Integer vector of UKB self-report illness codes.
#' @param disease_label Disease name label to assign to matched records.
#'
#' @return A data.table with filtered records and added \code{disease} column.
#'
#' @details
#' Common UKB self-report codes:
#' \itemize{
#'   \item 1065: High blood pressure
#'   \item 1066: Heart attack
#'   \item 1067: Angina
#'   \item 1068: Stroke
#'   \item 1220: Diabetes
#'   \item 1076: Heart failure
#' }
#'
#' @keywords internal
filter_self_report_codes <- function(sr_long, codes, disease_label) {
  if (!data.table::is.data.table(sr_long)) {
    sr_long <- data.table::as.data.table(sr_long)
  }
  result <- sr_long[sr_code %in% codes]
  result[, disease := disease_label]
  return(result)
}


#' @title Aggregate Earliest Self-Report Date Per Participant
#'
#' @description
#' Computes the earliest self-reported diagnosis date for each participant-disease combination.
#'
#' @param sr_filtered A data.table from \code{\link{filter_self_report_codes}}.
#'
#' @return A data.table with columns: \code{eid}, \code{disease},
#'   \code{earliest_date}, \code{source}.
#'
#' @keywords internal
aggregate_self_report_earliest <- function(sr_filtered) {
  if (!data.table::is.data.table(sr_filtered)) {
    sr_filtered <- data.table::as.data.table(sr_filtered)
  }
  if (!"diagnosis_date_quality" %in% names(sr_filtered)) {
    sr_filtered[, "diagnosis_date_quality" := data.table::fifelse(
      is.na(get("diag_date")),
      "unknown",
      "estimated"
    )]
  }
  result <- sr_filtered[
    ,
    {
      dated <- which(!is.na(diag_date))
      selected <- if (length(dated) > 0L) {
        dated[[which.min(diag_date[dated])]]
      } else {
        1L
      }
      list(
        earliest_date = diag_date[[selected]],
        source = if (is.na(diag_date[[selected]])) {
          "Self-report_unknown_date"
        } else {
          "Self-report"
        },
        diagnosis_date_quality = .SD[["diagnosis_date_quality"]][[selected]]
      )
    },
    by = .(eid, disease)
  ]
  return(result)
}
