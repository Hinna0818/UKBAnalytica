#' @title Parse ICD-9 Hospital Diagnosis Records
#'
#' @description
#' Extracts ICD-9 diagnosis codes from UK Biobank hospital inpatient data.
#' Converts the mixed-format storage (Python list string in p41271 + date array
#' in p41281_a*) into a standardized long-format data.table.
#'
#' @param dt A data.table or data.frame containing UKB data with columns:
#'   \code{eid}, \code{p41271}, and \code{p41281_a*} date columns.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{icd9_code}{ICD-9 diagnosis code}
#'     \item{diag_date}{Date of diagnosis}
#'     \item{source}{Data source identifier ("ICD9")}
#'   }
#'
#' @details
#' ICD-9 codes in UKB follow the format: 3-5 digits, optionally prefixed with V or E.
#' The function handles logical NA columns that may occur when all values are missing.
#'
#' @examples
#' \dontrun{
#' ukb_data <- data.table::fread("ukb_data.csv")
#' icd9_long <- parse_icd9_diagnoses(ukb_data)
#' }
#'
#' @import data.table
#' @importFrom stringi stri_extract_all_regex
#' @export
parse_icd9_diagnoses <- function(dt) {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  if (!"p41271" %in% names(dt)) {
    message("[parse_icd9_diagnoses] Column p41271 not found")
    return(data.table::data.table(
      eid = integer(0), icd9_code = character(0),
      diag_date = as.Date(character(0)), source = character(0)
    ))
  }

  # Handle logical NA columns by converting to character
  p41271_char <- as.character(dt$p41271)
  valid_rows <- which(!is.na(p41271_char) & nchar(p41271_char) > 0 & p41271_char != "NA")

  if (length(valid_rows) == 0) {
    message("[parse_icd9_diagnoses] No valid ICD-9 data found")
    return(data.table::data.table(
      eid = integer(0), icd9_code = character(0),
      diag_date = as.Date(character(0)), source = character(0)
    ))
  }

  # Step 1: Parse p41271 string list
  codes_long <- dt[valid_rows, ][
    ,
    {
      codes <- stringi::stri_extract_all_regex(as.character(p41271), "[VE]?[0-9]{3,5}")[[1]]
      if (length(codes) == 0 || all(is.na(codes))) {
        list(idx = integer(0), icd9_code = character(0))
      } else {
        list(idx = seq_along(codes) - 1L, icd9_code = codes)
      }
    },
    by = eid
  ]

  if (nrow(codes_long) == 0) {
    message("[parse_icd9_diagnoses] No ICD-9 codes found")
    return(data.table::data.table(
      eid = integer(0), icd9_code = character(0),
      diag_date = as.Date(character(0)), source = character(0)
    ))
  }

  # Step 2: Melt p41281_a* date columns
  date_cols <- grep("^p41281_a[0-9]+$", names(dt), value = TRUE)
  if (length(date_cols) == 0) {
    warning("[parse_icd9_diagnoses] No p41281_a* date columns found")
    codes_long[, `:=`(diag_date = as.Date(NA), source = "ICD9")]
    return(codes_long[, .(eid, icd9_code, diag_date, source)])
  }

  eids_with_codes <- unique(codes_long$eid)
  dt_sub <- dt[eid %in% eids_with_codes, c("eid", date_cols), with = FALSE]
  dt_sub[, (date_cols) := lapply(.SD, as.character), .SDcols = date_cols]

  dates_long <- data.table::melt(
    dt_sub, id.vars = "eid", measure.vars = date_cols,
    variable.name = "date_col", value.name = "diag_date", na.rm = FALSE
  )
  dates_long[, idx := as.integer(sub("^p41281_a", "", date_col))]
  dates_long[, date_col := NULL]
  dates_long[, diag_date := as.Date(diag_date)]

  # Step 3: Join by eid and index
  result <- data.table::merge.data.table(
    codes_long, dates_long, by = c("eid", "idx"), all.x = TRUE
  )
  result[, idx := NULL]
  result[, source := "ICD9"]
  data.table::setorder(result, eid, diag_date, na.last = TRUE)

  return(result[, .(eid, icd9_code, diag_date, source)])
}


#' @title Filter ICD-9 Records by Code Pattern
#'
#' @description
#' Filters ICD-9 diagnosis records using regular expression pattern matching.
#'
#' @param icd9_long A data.table from \code{\link{parse_icd9_diagnoses}}.
#' @param pattern Regular expression pattern for ICD-9 codes.
#' @param disease_label Disease name label to assign to matched records.
#'
#' @return A data.table with filtered records and added \code{disease} column.
#'
#' @keywords internal
filter_icd9_codes <- function(icd9_long, pattern, disease_label) {
  if (!data.table::is.data.table(icd9_long)) {
    icd9_long <- data.table::as.data.table(icd9_long)
  }
  result <- icd9_long[grepl(pattern, icd9_code, perl = TRUE)]
  result[, disease := disease_label]
  return(result)
}


#' @title Aggregate Earliest ICD-9 Diagnosis Date Per Participant
#'
#' @description
#' Computes the earliest diagnosis date for each participant-disease combination.
#'
#' @param icd9_filtered A data.table from \code{\link{filter_icd9_codes}}.
#'
#' @return A data.table with columns: \code{eid}, \code{disease},
#'   \code{earliest_date}, \code{source}.
#'
#' @keywords internal
aggregate_icd9_earliest <- function(icd9_filtered) {
  if (!data.table::is.data.table(icd9_filtered)) {
    icd9_filtered <- data.table::as.data.table(icd9_filtered)
  }
  result <- icd9_filtered[
    !is.na(diag_date),
    .(earliest_date = min(diag_date), source = "ICD9"),
    by = .(eid, disease)
  ]
  return(result)
}



#' @title Parse ICD-10 Hospital Diagnosis Records
#'
#' @description
#' Extracts ICD-10 diagnosis codes from UK Biobank hospital inpatient data.
#' Converts the mixed-format storage (Python list string in p41270 + date array
#' in p41280_a*) into a standardized long-format data.table.
#'
#' @param dt A data.table or data.frame containing UKB data with columns:
#'   \code{eid}, \code{p41270}, and \code{p41280_a*} date columns.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{icd10_code}{ICD-10 diagnosis code}
#'     \item{diag_date}{Date of diagnosis}
#'     \item{source}{Data source identifier ("ICD10")}
#'   }
#'
#' @details
#' The function implements the Index-Match logic specified in the UKB data dictionary:
#' the k-th element in p41270 corresponds to date column p41280_a(k-1) (0-indexed).
#'
#' Processing pipeline:
#' \enumerate{
#'   \item Parse Python list string format in p41270
#'   \item Melt p41280_a* date columns to long format
#'   \item Join codes and dates by eid and positional index
#' }
#'
#' @examples
#' \dontrun{
#' ukb_data <- data.table::fread("ukb_data.csv")
#' icd10_long <- parse_icd10_diagnoses(ukb_data)
#' }
#'
#' @import data.table
#' @importFrom stringi stri_extract_all_regex
#' @export
parse_icd10_diagnoses <- function(dt) {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  # Step 1: Parse p41270 string list to long format
  codes_long <- dt[
    !is.na(p41270) & p41270 != "",
    {
      codes <- stringi::stri_extract_all_regex(p41270, "[A-Z][0-9]{2,3}")[[1]]
      if (length(codes) == 0 || all(is.na(codes))) {
        list(idx = integer(0), icd10_code = character(0))
      } else {
        list(idx = seq_along(codes) - 1L, icd10_code = codes)
      }
    },
    by = eid
  ]

  if (nrow(codes_long) == 0) {
    message("[parse_icd10_diagnoses] No ICD-10 codes found")
    return(data.table::data.table(
      eid = integer(0), icd10_code = character(0),
      diag_date = as.Date(character(0)), source = character(0)
    ))
  }

  # Step 2: Melt p41280_a* date columns
  date_cols <- grep("^p41280_a[0-9]+$", names(dt), value = TRUE)
  if (length(date_cols) == 0) {
    warning("[parse_icd10_diagnoses] No p41280_a* date columns found")
    codes_long[, `:=`(diag_date = as.Date(NA), source = "ICD10")]
    return(codes_long[, .(eid, icd10_code, diag_date, source)])
  }

  eids_with_codes <- unique(codes_long$eid)
  dt_sub <- dt[eid %in% eids_with_codes, c("eid", date_cols), with = FALSE]
  dt_sub[, (date_cols) := lapply(.SD, as.character), .SDcols = date_cols]

  dates_long <- data.table::melt(
    dt_sub, id.vars = "eid", measure.vars = date_cols,
    variable.name = "date_col", value.name = "diag_date", na.rm = FALSE
  )
  dates_long[, idx := as.integer(sub("^p41280_a", "", date_col))]
  dates_long[, date_col := NULL]
  dates_long[, diag_date := as.Date(diag_date)]

  # Step 3: Join by eid and index
  result <- data.table::merge.data.table(
    codes_long, dates_long, by = c("eid", "idx"), all.x = TRUE
  )
  result[, idx := NULL]
  result[, source := "ICD10"]
  data.table::setorder(result, eid, diag_date, na.last = TRUE)

  return(result[, .(eid, icd10_code, diag_date, source)])
}


#' @title Filter ICD-10 Records by Code Pattern
#'
#' @description
#' Filters ICD-10 diagnosis records using regular expression pattern matching.
#'
#' @param icd10_long A data.table from \code{\link{parse_icd10_diagnoses}}.
#' @param pattern Regular expression pattern for ICD-10 codes (e.g., "^I71" for aortic aneurysm).
#' @param disease_label Disease name label to assign to matched records.
#'
#' @return A data.table with filtered records and added \code{disease} column.
#'
#' @examples
#' \dontrun{
#' icd10_long <- parse_icd10_diagnoses(ukb_data)
#' aa_cases <- filter_icd10_codes(icd10_long, "^I71", "Aortic_Aneurysm")
#' }
#'
#' @keywords internal
filter_icd10_codes <- function(icd10_long, pattern, disease_label) {
  if (!data.table::is.data.table(icd10_long)) {
    icd10_long <- data.table::as.data.table(icd10_long)
  }
  result <- icd10_long[grepl(pattern, icd10_code, perl = TRUE)]
  result[, disease := disease_label]
  return(result)
}


#' @title Aggregate Earliest ICD-10 Diagnosis Date Per Participant
#'
#' @description
#' Computes the earliest diagnosis date for each participant-disease combination.
#' Essential for determining incident vs prevalent cases in survival analysis.
#'
#' @param icd10_filtered A data.table from \code{\link{filter_icd10_codes}}.
#'
#' @return A data.table with columns: \code{eid}, \code{disease},
#'   \code{earliest_date}, \code{source}.
#'
#' @keywords internal
aggregate_icd10_earliest <- function(icd10_filtered) {
  if (!data.table::is.data.table(icd10_filtered)) {
    icd10_filtered <- data.table::as.data.table(icd10_filtered)
  }
  result <- icd10_filtered[
    !is.na(diag_date),
    .(earliest_date = min(diag_date), source = "ICD10"),
    by = .(eid, disease)
  ]
  return(result)
}
