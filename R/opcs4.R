#' @title Parse OPCS4 Hospital Procedure Records
#'
#' @description
#' Extracts OPCS4 operative procedure codes from UK Biobank hospital inpatient
#' summary operations data. Supports the common export shape where
#' \code{p41272} stores a list-string of codes and \code{p41282_a*} stores the
#' corresponding dates, while also tolerating expanded \code{p41272_a*} columns.
#'
#' @param dt A data.table or data.frame containing UKB data with columns:
#'   \code{eid}, and either \code{p41272} or \code{p41272_a*}, plus
#'   \code{p41282_a*} date columns.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{opcs4_code}{OPCS4 procedure code}
#'     \item{diag_date}{Date of first recorded procedure for that code/index}
#'     \item{source}{Data source identifier ("OPCS4")}
#'   }
#'
#' @details
#' The function implements the same index-matching logic used for UKB summary
#' diagnosis fields: the k-th procedure code in \code{p41272} corresponds to the
#' date stored in \code{p41282_a(k-1)} (0-indexed).
#'
#' @examples
#' \dontrun{
#' ukb_data <- data.table::fread("ukb_data.csv")
#' opcs4_long <- parse_opcs4_procedures(ukb_data)
#' }
#'
#' @export
parse_opcs4_procedures <- function(dt) {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  code_cols <- grep("^p41272_a[0-9]+$", names(dt), value = TRUE)

  if ("p41272" %in% names(dt)) {
    p41272_char <- as.character(dt$p41272)
    valid_rows <- which(!is.na(p41272_char) & nchar(p41272_char) > 0 & p41272_char != "NA")

    if (length(valid_rows) == 0) {
      message("[parse_opcs4_procedures] No valid OPCS4 data found")
      return(data.table::data.table(
        eid = integer(0), opcs4_code = character(0),
        diag_date = as.Date(character(0)), source = character(0)
      ))
    }

    codes_long <- dt[valid_rows, ][
      ,
      {
        codes <- stringi::stri_extract_all_regex(as.character(p41272), "[A-Z][0-9]{3}")[[1]]
        if (length(codes) == 0 || all(is.na(codes))) {
          list(idx = integer(0), opcs4_code = character(0))
        } else {
          list(idx = seq_along(codes) - 1L, opcs4_code = codes)
        }
      },
      by = eid
    ]
  } else if (length(code_cols) > 0) {
    dt_sub <- dt[, c("eid", code_cols), with = FALSE]
    dt_sub[, (code_cols) := lapply(.SD, as.character), .SDcols = code_cols]

    codes_long <- data.table::melt(
      dt_sub, id.vars = "eid", measure.vars = code_cols,
      variable.name = "code_col", value.name = "opcs4_code", na.rm = TRUE
    )
    codes_long <- codes_long[!is.na(opcs4_code) & nzchar(opcs4_code) & opcs4_code != "NA"]
    codes_long[, idx := as.integer(sub("^p41272_a", "", code_col))]
    codes_long[, code_col := NULL]
    codes_long[, opcs4_code := stringi::stri_extract_first_regex(opcs4_code, "[A-Z][0-9]{3}")]
    codes_long <- codes_long[!is.na(opcs4_code)]
  } else {
    message("[parse_opcs4_procedures] No p41272 or p41272_a* columns found")
    return(data.table::data.table(
      eid = integer(0), opcs4_code = character(0),
      diag_date = as.Date(character(0)), source = character(0)
    ))
  }

  if (nrow(codes_long) == 0) {
    message("[parse_opcs4_procedures] No OPCS4 codes found")
    return(data.table::data.table(
      eid = integer(0), opcs4_code = character(0),
      diag_date = as.Date(character(0)), source = character(0)
    ))
  }

  date_cols <- grep("^p41282_a[0-9]+$", names(dt), value = TRUE)
  if (length(date_cols) == 0) {
    warning("[parse_opcs4_procedures] No p41282_a* date columns found")
    codes_long[, `:=`(diag_date = as.Date(NA), source = "OPCS4")]
    return(codes_long[, .(eid, opcs4_code, diag_date, source)])
  }

  eids_with_codes <- unique(codes_long$eid)
  dt_sub <- dt[eid %in% eids_with_codes, c("eid", date_cols), with = FALSE]
  dt_sub[, (date_cols) := lapply(.SD, as.character), .SDcols = date_cols]

  dates_long <- data.table::melt(
    dt_sub, id.vars = "eid", measure.vars = date_cols,
    variable.name = "date_col", value.name = "diag_date", na.rm = TRUE
  )
  dates_long[, idx := as.integer(sub("^p41282_a", "", date_col))]
  dates_long[, date_col := NULL]
  dates_long[, diag_date := .safe_as_date(diag_date, col_name = "OPCS4_diag_date")]

  result <- data.table::merge.data.table(
    codes_long, dates_long, by = c("eid", "idx"), all.x = TRUE
  )
  result[, idx := NULL]
  result[, source := "OPCS4"]
  data.table::setorder(result, eid, diag_date, na.last = TRUE)

  result[, .(eid, opcs4_code, diag_date, source)]
}


#' @title Filter OPCS4 Procedure Records by Code Pattern
#'
#' @description
#' Filters OPCS4 procedure records using regular expression pattern matching.
#'
#' @param opcs4_long A data.table from \code{\link{parse_opcs4_procedures}}.
#' @param pattern Regular expression pattern for OPCS4 procedure codes.
#' @param disease_label Disease name label to assign to matched records.
#'
#' @return A data.table with filtered records and added \code{disease} column.
#'
#' @keywords internal
filter_opcs4_codes <- function(opcs4_long, pattern, disease_label) {
  if (!data.table::is.data.table(opcs4_long)) {
    opcs4_long <- data.table::as.data.table(opcs4_long)
  }
  result <- opcs4_long[grepl(pattern, opcs4_code, perl = TRUE)]
  result[, disease := disease_label]
  result
}


#' @title Aggregate Earliest OPCS4 Procedure Date Per Participant
#'
#' @description
#' Computes the earliest procedure date for each participant-disease combination.
#'
#' @param opcs4_filtered A data.table from \code{\link{filter_opcs4_codes}}.
#'
#' @return A data.table with columns: \code{eid}, \code{disease},
#'   \code{earliest_date}, \code{source}.
#'
#' @keywords internal
aggregate_opcs4_earliest <- function(opcs4_filtered) {
  if (!data.table::is.data.table(opcs4_filtered)) {
    opcs4_filtered <- data.table::as.data.table(opcs4_filtered)
  }
  result <- opcs4_filtered[
    !is.na(diag_date),
    .(earliest_date = min(diag_date), source = "OPCS4"),
    by = .(eid, disease)
  ]
  result
}