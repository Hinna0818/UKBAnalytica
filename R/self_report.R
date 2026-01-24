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
#' @examples
#' \dontrun{
#' ukb_data <- data.table::fread("ukb_data.csv")
#' sr_long <- parse_self_reported_illnesses(ukb_data)
#' }
#'
#' @import data.table
#' @export
parse_self_reported_illnesses <- function(dt, baseline_col = "p53_i0") {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  # Step 1: Melt p20002_i*_a* code columns
  code_cols <- grep("^p20002_i[0-9]+_a[0-9]+$", names(dt), value = TRUE)
  if (length(code_cols) == 0) {
    message("[parse_self_reported_illnesses] No p20002_i*_a* columns found")
    return(data.table::data.table(
      eid = integer(0), sr_code = integer(0),
      diag_date = as.Date(character(0)), source = character(0),
      instance = integer(0), array_idx = integer(0)
    ))
  }

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
      instance = integer(0), array_idx = integer(0)
    ))
  }

  # Step 2: Melt p20008_i*_a* year columns
  year_cols <- grep("^p20008_i[0-9]+_a[0-9]+$", names(dt), value = TRUE)
  if (length(year_cols) == 0) {
    warning("[parse_self_reported_illnesses] No p20008_i*_a* year columns found")
    codes_long[, `:=`(diag_date = as.Date(NA), source = "Self-report")]
    return(codes_long)
  }

  eids_with_codes <- unique(codes_long$eid)
  dt_sub <- dt[eid %in% eids_with_codes, c("eid", year_cols), with = FALSE]

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

  # Step 3: Join codes and years
  result <- data.table::merge.data.table(
    codes_long, years_long,
    by = c("eid", "instance", "array_idx"), all.x = TRUE
  )

  # Step 4: Convert fractional year to date
  result[diag_year < 0, diag_year := NA_real_]
  result[, diag_date := {
    year_int <- floor(diag_year)
    month_frac <- diag_year - year_int
    month_int <- pmax(1L, pmin(12L, as.integer(round(month_frac * 12)) + 1L))
    date_str <- sprintf("%04d-%02d-01", year_int, month_int)
    as.Date(date_str)
  }]

  result[, diag_year := NULL]
  result[, source := "Self-report"]
  data.table::setorder(result, eid, diag_date, na.last = TRUE)

  return(result[, .(eid, sr_code, diag_date, source, instance, array_idx)])
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
#' @export
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
#' @export
aggregate_self_report_earliest <- function(sr_filtered) {
  if (!data.table::is.data.table(sr_filtered)) {
    sr_filtered <- data.table::as.data.table(sr_filtered)
  }
  result <- sr_filtered[
    !is.na(diag_date),
    .(earliest_date = min(diag_date), source = "Self-report"),
    by = .(eid, disease)
  ]
  return(result)
}
