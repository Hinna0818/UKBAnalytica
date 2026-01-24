#' @title Parse Death Registry Records
#'
#' @description
#' Extracts death registry data from UK Biobank linked mortality records.
#' Parses both primary (p40001) and contributing (p40002) causes of death
#' along with death dates (p40000).
#'
#' @param dt A data.table or data.frame containing UKB data with columns:
#'   \code{eid}, \code{p40001_i*}, \code{p40002_i*_a*}, and \code{p40000_i*}.
#'
#' @return A data.table with columns:
#'   \describe{
#'     \item{eid}{Participant identifier}
#'     \item{death_code}{ICD-10 cause of death code}
#'     \item{death_date}{Date of death}
#'     \item{source}{Data source identifier ("Death")}
#'     \item{cause_type}{"primary" or "secondary"}
#'   }
#'
#' @details
#' Death causes serve as definitive diagnosis confirmation. If a participant
#' died from a specific disease, the death date becomes the diagnosis date
#' for that condition (if not previously diagnosed).
#'
#' @examples
#' \dontrun{
#' ukb_data <- data.table::fread("ukb_data.csv")
#' death_long <- parse_death_records(ukb_data)
#' }
#'
#' @import data.table
#' @export
parse_death_records <- function(dt) {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  # Step 1: Extract death dates
  death_date_cols <- grep("^p40000_i[0-9]+$", names(dt), value = TRUE)
  if (length(death_date_cols) == 0) {
    message("[parse_death_records] No p40000 death date columns found")
    return(data.table::data.table(
      eid = integer(0), death_code = character(0),
      death_date = as.Date(character(0)), source = character(0),
      cause_type = character(0)
    ))
  }

  death_dates <- dt[, c("eid", death_date_cols), with = FALSE]
  death_dates[, (death_date_cols) := lapply(.SD, as.character), .SDcols = death_date_cols]
  death_dates_long <- data.table::melt(
    death_dates, id.vars = "eid", measure.vars = death_date_cols,
    variable.name = "col", value.name = "death_date", na.rm = TRUE
  )
  death_dates_long[, death_date := as.Date(death_date)]
  death_dates_long[, col := NULL]

  # Aggregate to earliest death date per participant
  death_dates_agg <- death_dates_long[
    !is.na(death_date),
    .(death_date = min(death_date)),
    by = eid
  ]

  if (nrow(death_dates_agg) == 0) {
    message("[parse_death_records] No death records found")
    return(data.table::data.table(
      eid = integer(0), death_code = character(0),
      death_date = as.Date(character(0)), source = character(0),
      cause_type = character(0)
    ))
  }

  # Step 2: Parse primary causes (p40001_i*)
  primary_cols <- grep("^p40001_i[0-9]+$", names(dt), value = TRUE)
  primary_long <- NULL

  if (length(primary_cols) > 0) {
    dt[, (primary_cols) := lapply(.SD, as.character), .SDcols = primary_cols]
    primary_long <- data.table::melt(
      dt[, c("eid", primary_cols), with = FALSE],
      id.vars = "eid", measure.vars = primary_cols,
      variable.name = "col", value.name = "death_code", na.rm = TRUE
    )
    primary_long <- primary_long[!is.na(death_code) & death_code != ""]
    primary_long[, col := NULL]
    primary_long[, cause_type := "primary"]
  }

  # Step 3: Parse secondary causes (p40002_i*_a*)
  secondary_cols <- grep("^p40002_i[0-9]+_a[0-9]+$", names(dt), value = TRUE)
  secondary_long <- NULL

  if (length(secondary_cols) > 0) {
    dt[, (secondary_cols) := lapply(.SD, as.character), .SDcols = secondary_cols]
    secondary_long <- data.table::melt(
      dt[, c("eid", secondary_cols), with = FALSE],
      id.vars = "eid", measure.vars = secondary_cols,
      variable.name = "col", value.name = "death_code", na.rm = TRUE
    )
    secondary_long <- secondary_long[!is.na(death_code) & death_code != ""]
    secondary_long[, col := NULL]
    secondary_long[, cause_type := "secondary"]
  }

  # Step 4: Combine all causes
  all_causes <- data.table::rbindlist(
    list(primary_long, secondary_long),
    use.names = TRUE, fill = TRUE
  )

  if (nrow(all_causes) == 0) {
    message("[parse_death_records] No death cause codes found")
    return(data.table::data.table(
      eid = integer(0), death_code = character(0),
      death_date = as.Date(character(0)), source = character(0),
      cause_type = character(0)
    ))
  }

  all_causes <- unique(all_causes, by = c("eid", "death_code"))

  # Step 5: Join with death dates
  result <- data.table::merge.data.table(
    all_causes, death_dates_agg, by = "eid", all.x = TRUE
  )
  result[, source := "Death"]
  data.table::setorder(result, eid, cause_type, death_code)

  return(result[, .(eid, death_code, death_date, source, cause_type)])
}


#' @title Filter Death Records by ICD-10 Code Pattern
#'
#' @description
#' Filters death cause records using regular expression pattern matching.
#'
#' @param death_long A data.table from \code{\link{parse_death_records}}.
#' @param pattern Regular expression pattern for ICD-10 death codes.
#' @param disease_label Disease name label to assign to matched records.
#'
#' @return A data.table with filtered records and added \code{disease} column.
#'
#' @export
filter_death_codes <- function(death_long, pattern, disease_label) {
  if (!data.table::is.data.table(death_long)) {
    death_long <- data.table::as.data.table(death_long)
  }
  result <- death_long[grepl(pattern, death_code, perl = TRUE)]
  result[, disease := disease_label]
  return(result)
}


#' @title Aggregate Death as Diagnosis Source
#'
#' @description
#' Uses death date as diagnosis date for participants who died from the target condition.
#'
#' @param death_filtered A data.table from \code{\link{filter_death_codes}}.
#'
#' @return A data.table with columns: \code{eid}, \code{disease},
#'   \code{earliest_date}, \code{source}.
#'
#' @export
aggregate_death_as_diagnosis <- function(death_filtered) {
  if (!data.table::is.data.table(death_filtered)) {
    death_filtered <- data.table::as.data.table(death_filtered)
  }
  result <- death_filtered[
    !is.na(death_date),
    .(earliest_date = min(death_date), source = "Death"),
    by = .(eid, disease)
  ]
  return(result)
}


#' @title Extract Death Dates for All Deceased Participants
#'
#' @description
#' Returns death dates for all deceased participants, used for censoring
#' in survival analysis.
#'
#' @param dt A data.table or data.frame containing UKB data.
#'
#' @return A data.table with columns: \code{eid}, \code{death_date}.
#'
#' @export
get_death_dates <- function(dt) {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }

  death_date_cols <- grep("^p40000_i[0-9]+$", names(dt), value = TRUE)
  if (length(death_date_cols) == 0) {
    return(data.table::data.table(
      eid = integer(0), death_date = as.Date(character(0))
    ))
  }

  dt[, (death_date_cols) := lapply(.SD, as.character), .SDcols = death_date_cols]

  death_dates <- data.table::melt(
    dt[, c("eid", death_date_cols), with = FALSE],
    id.vars = "eid", measure.vars = death_date_cols,
    variable.name = "col", value.name = "death_date", na.rm = TRUE
  )
  death_dates[, death_date := as.Date(death_date)]
  death_dates[, col := NULL]

  result <- death_dates[
    !is.na(death_date),
    .(death_date = min(death_date)),
    by = eid
  ]
  return(result)
}
