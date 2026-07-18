#' Build a UK Biobank follow-up time skeleton
#'
#' @description
#' Creates a reusable participant-level time skeleton for prospective UK
#' Biobank analyses. The function standardizes baseline date, approximate birth
#' date, age at baseline, death date, loss-to-follow-up date, administrative
#' censoring date, follow-up end date, and follow-up time. It does not define
#' disease outcomes; instead, it provides a common time basis that can be reused
#' by endpoint-specific functions such as \code{\link{build_survival_dataset}}.
#'
#' @param data A data.frame or data.table containing UK Biobank columns.
#' @param id_col Participant identifier column. Default \code{"eid"}.
#' @param baseline_col Baseline assessment date column. Default \code{"p53_i0"}.
#' @param birth_year_col Year-of-birth column. Default \code{"p34"}.
#' @param birth_month_col Month-of-birth column. Default \code{"p52"}.
#' @param age_col Age-at-baseline column. Default \code{"p21022"}. If missing,
#'   age is approximated from baseline date and birth year/month when available.
#' @param death_date_cols Death date columns or a regular expression used to
#'   identify them. Default \code{"^(participant\\.)?p40000_i[0-9]+$"}.
#' @param lost_to_followup_col Optional date lost to follow-up column. Default
#'   \code{"p191"}.
#' @param admin_censor_date Administrative censoring date.
#' @param keep_source_dates Logical. If \code{FALSE}, source dates used to define
#'   censoring are removed from the output.
#'
#' @return A data.table with one row per participant and standardized follow-up
#'   time fields.
#' @export
#'
#' @examples
#' demo <- data.frame(
#'   eid = 1:3,
#'   p53_i0 = as.Date(c("2010-01-01", "2011-01-01", "2012-01-01")),
#'   p21022 = c(50, 60, 70),
#'   p40000_i0 = as.Date(c(NA, "2015-01-01", NA))
#' )
#'
#' ukb_time_skeleton(demo, admin_censor_date = as.Date("2020-12-31"))
#'
#' @import data.table
#' @export
ukb_time_skeleton <- function(data,
                              id_col = "eid",
                              baseline_col = "p53_i0",
                              birth_year_col = "p34",
                              birth_month_col = "p52",
                              age_col = "p21022",
                              death_date_cols = "^(participant\\.)?p40000_i[0-9]+$",
                              lost_to_followup_col = "p191",
                              admin_censor_date = as.Date("2023-10-31"),
                              keep_source_dates = TRUE) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("'id_col' must be a single non-empty column name.", call. = FALSE)
  }
  if (!is.character(baseline_col) || length(baseline_col) != 1L || !nzchar(baseline_col)) {
    stop("'baseline_col' must be a single non-empty column name.", call. = FALSE)
  }
  resolved_id_col <- .ukb_time_resolve_single_col(data, id_col)
  resolved_baseline_col <- .ukb_time_resolve_single_col(data, baseline_col)

  if (is.na(resolved_id_col)) {
    stop("ID column not found: ", id_col, call. = FALSE)
  }
  if (is.na(resolved_baseline_col)) {
    stop("Baseline date column not found: ", baseline_col, call. = FALSE)
  }

  admin_censor_date <- .safe_as_date(admin_censor_date, col_name = "admin_censor_date")
  if (length(admin_censor_date) != 1L || is.na(admin_censor_date)) {
    stop("'admin_censor_date' must be coercible to a single non-missing Date.", call. = FALSE)
  }

  dt <- data.table::as.data.table(data.table::copy(data))
  .ukb_time_validate_ids(dt[[resolved_id_col]], resolved_id_col)
  out <- data.table::data.table(
    eid = dt[[resolved_id_col]],
    baseline_date = .safe_as_date(dt[[resolved_baseline_col]], col_name = resolved_baseline_col)
  )

  out[, `:=`(
    birth_year = .ukb_time_integer_col(dt, birth_year_col),
    birth_month = .ukb_time_integer_col(dt, birth_month_col)
  )]

  out[, birth_date_approx := as.Date(NA)]
  valid_birth_month <- !is.na(out$birth_year) &
    !is.na(out$birth_month) &
    out$birth_month >= 1L &
    out$birth_month <= 12L
  if (any(valid_birth_month)) {
    out[valid_birth_month, birth_date_approx := as.Date(sprintf(
      "%04d-%02d-15",
      birth_year,
      birth_month
    ))]
  }

  resolved_age_col <- .ukb_time_resolve_single_col(dt, age_col)
  age_from_col <- .ukb_time_numeric_col(dt, age_col)
  age_from_birth <- as.numeric(out$baseline_date - out$birth_date_approx) / 365.25
  out[, age_at_baseline := ifelse(!is.na(age_from_col), age_from_col, age_from_birth)]
  out[, age_at_baseline_source := ifelse(
    !is.na(age_from_col),
    ifelse(!is.na(resolved_age_col), resolved_age_col, "none"),
    ifelse(!is.na(age_from_birth), "birth_year_month", "missing")
  )]

  death_dates <- .ukb_time_extract_death_dates(dt, id_col = resolved_id_col, death_date_cols = death_date_cols)
  out <- data.table::merge.data.table(out, death_dates, by = "eid", all.x = TRUE)

  lost_col <- .ukb_time_resolve_single_col(dt, lost_to_followup_col)
  if (!is.na(lost_col)) {
    lost_dt <- data.table::data.table(
      eid = dt[[resolved_id_col]],
      lost_to_followup_date = .safe_as_date(dt[[lost_col]], col_name = lost_col)
    )
  } else {
    lost_dt <- data.table::data.table(
      eid = dt[[resolved_id_col]],
      lost_to_followup_date = as.Date(NA)
    )
  }
  out <- data.table::merge.data.table(out, lost_dt, by = "eid", all.x = TRUE, sort = FALSE)

  out[, admin_censor_date := admin_censor_date]
  out[, followup_end_date := pmin(death_date, lost_to_followup_date, admin_censor_date, na.rm = TRUE)]
  out[is.na(followup_end_date), followup_end_date := admin_censor_date]

  out[, followup_end_reason := data.table::fcase(
    !is.na(death_date) & followup_end_date == death_date, "death",
    !is.na(lost_to_followup_date) & followup_end_date == lost_to_followup_date, "lost_to_followup",
    followup_end_date == admin_censor_date, "administrative_censoring",
    default = "unknown"
  )]

  out[, followup_time_days := as.numeric(followup_end_date - baseline_date)]
  out[, followup_time_years := followup_time_days / 365.25]
  out[, `:=`(
    baseline_missing = is.na(baseline_date),
    death_before_baseline = !is.na(death_date) & !is.na(baseline_date) & death_date < baseline_date,
    lost_to_followup_before_baseline = !is.na(lost_to_followup_date) & !is.na(baseline_date) & lost_to_followup_date < baseline_date,
    valid_followup = !is.na(baseline_date) & !is.na(followup_end_date) & followup_time_days >= 0
  )]

  if (!isTRUE(keep_source_dates)) {
    out[, c("death_date", "lost_to_followup_date", "admin_censor_date") := NULL]
  }

  data.table::setorder(out, eid)
  attr(out, "time_skeleton_summary") <- .ukb_time_skeleton_summary(out)
  class(out) <- c("ukb_time_skeleton", class(out))
  out[]
}

#' Create the participant-level observation window used by endpoint functions
#' @keywords internal
#' @noRd
.ukb_followup_window <- function(dt,
                                 baseline_col = "p53_i0",
                                 censor_date = as.Date("2023-10-31")) {
  dt <- data.table::as.data.table(dt)
  if (!"eid" %in% names(dt)) {
    stop("Column 'eid' was not found in the input data.", call. = FALSE)
  }
  if (!baseline_col %in% names(dt)) {
    stop("Baseline column not found: ", baseline_col, call. = FALSE)
  }
  .ukb_time_validate_ids(dt[["eid"]], "eid")
  censor_date <- .safe_as_date(censor_date, col_name = "censor_date")
  if (length(censor_date) != 1L || is.na(censor_date)) {
    stop("`censor_date` must be a single non-missing Date.", call. = FALSE)
  }

  baseline_dt <- data.table::data.table(
    eid = dt[["eid"]],
    baseline_date = .safe_as_date(dt[[baseline_col]], col_name = baseline_col)
  )
  death_dates <- get_death_dates(dt)

  lost_candidates <- c(
    "p191", "p191_i0", "participant.p191", "participant.p191_i0"
  )
  lost_col <- intersect(lost_candidates, names(dt))
  if (length(lost_col) > 0L) {
    lost_dt <- data.table::data.table(
      eid = dt[["eid"]],
      lost_to_followup_date = .safe_as_date(dt[[lost_col[[1]]]], col_name = lost_col[[1]])
    )
  } else {
    lost_dt <- data.table::data.table(
      eid = dt[["eid"]],
      lost_to_followup_date = as.Date(NA)
    )
  }

  out <- data.table::merge.data.table(
    baseline_dt,
    death_dates,
    by = "eid",
    all.x = TRUE,
    sort = FALSE
  )
  out <- data.table::merge.data.table(
    out,
    lost_dt,
    by = "eid",
    all.x = TRUE,
    sort = FALSE
  )

  skeleton_cols <- grep(
    "^\\.ukba_time_skeleton_followup_end_date",
    names(dt),
    value = TRUE
  )
  if (length(skeleton_cols) > 1L) {
    stop("Multiple internal follow-up end columns were found.", call. = FALSE)
  }
  if (length(skeleton_cols) == 1L) {
    skeleton_dt <- data.table::data.table(
      eid = dt[["eid"]],
      skeleton_followup_end_date = .safe_as_date(
        dt[[skeleton_cols[[1]]]],
        col_name = skeleton_cols[[1]]
      )
    )
    out <- data.table::merge.data.table(
      out,
      skeleton_dt,
      by = "eid",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    data.table::set(out, j = "skeleton_followup_end_date", value = as.Date(NA))
  }

  data.table::set(out, j = "calculated_followup_end_date", value = pmin(
    out[["death_date"]],
    out[["lost_to_followup_date"]],
    censor_date,
    na.rm = TRUE
  ))
  data.table::set(out, j = "followup_end_date", value = data.table::fifelse(
    !is.na(out[["skeleton_followup_end_date"]]),
    out[["skeleton_followup_end_date"]],
    out[["calculated_followup_end_date"]]
  ))
  out[is.na(followup_end_date), followup_end_date := censor_date]
  out[, c("skeleton_followup_end_date", "calculated_followup_end_date") := NULL]
  out[]
}

.ukb_time_validate_ids <- function(x, column_name = "eid") {
  if (anyNA(x)) {
    stop("Participant ID column '", column_name, "' contains missing values.", call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop(
      "Participant ID column '", column_name,
      "' must contain one unique row per participant.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.ukb_time_col_exists <- function(dt, col) {
  !is.na(.ukb_time_resolve_single_col(dt, col))
}

.ukb_time_integer_col <- function(dt, col) {
  resolved_col <- .ukb_time_resolve_single_col(dt, col)
  if (!is.na(resolved_col)) {
    return(suppressWarnings(as.integer(dt[[resolved_col]])))
  }
  rep(NA_integer_, nrow(dt))
}

.ukb_time_numeric_col <- function(dt, col) {
  resolved_col <- .ukb_time_resolve_single_col(dt, col)
  if (!is.na(resolved_col)) {
    return(suppressWarnings(as.numeric(dt[[resolved_col]])))
  }
  rep(NA_real_, nrow(dt))
}

.ukb_time_resolve_single_col <- function(dt, col) {
  if (is.null(col) || length(col) != 1L || is.na(col) || !nzchar(col)) {
    return(NA_character_)
  }
  if (col %in% names(dt)) {
    return(col)
  }
  hits <- grep(col, names(dt), value = TRUE)
  if (length(hits) == 1L) {
    return(hits[[1]])
  }
  NA_character_
}

.ukb_time_resolve_date_cols <- function(dt, cols) {
  if (is.null(cols) || length(cols) == 0L) {
    return(character(0))
  }
  exact <- intersect(cols, names(dt))
  if (length(exact) > 0L) {
    return(exact)
  }
  if (length(cols) == 1L && !is.na(cols) && nzchar(cols)) {
    return(grep(cols, names(dt), value = TRUE))
  }
  character(0)
}

.ukb_time_extract_death_dates <- function(dt, id_col, death_date_cols) {
  cols <- .ukb_time_resolve_date_cols(dt, death_date_cols)
  if (length(cols) == 0L) {
    return(data.table::data.table(eid = dt[[id_col]], death_date = as.Date(NA)))
  }

  death_dt <- dt[, c(id_col, cols), with = FALSE]
  data.table::setnames(death_dt, id_col, "eid")
  death_long <- data.table::melt(
    death_dt,
    id.vars = "eid",
    measure.vars = cols,
    variable.name = "death_date_col",
    value.name = "death_date",
    na.rm = TRUE
  )
  death_long[, death_date := .safe_as_date(death_date, col_name = "death_date")]
  if (nrow(death_long) == 0L || all(is.na(death_long[["death_date"]]))) {
    return(data.table::data.table(
      eid = dt[[id_col]],
      death_date = as.Date(NA)
    ))
  }
  death_agg <- death_long[
    !is.na(death_date),
    .(death_date = min(death_date)),
    by = eid
  ]
  out <- data.table::data.table(eid = dt[[id_col]])
  out <- data.table::merge.data.table(out, death_agg, by = "eid", all.x = TRUE)
  out
}

.ukb_time_skeleton_summary <- function(x) {
  data.table::data.table(
    metric = c(
      "n",
      "valid_followup",
      "missing_baseline_date",
      "death_before_baseline",
      "lost_to_followup_before_baseline",
      "administrative_censoring",
      "death",
      "lost_to_followup"
    ),
    value = c(
      nrow(x),
      sum(x$valid_followup, na.rm = TRUE),
      sum(x$baseline_missing, na.rm = TRUE),
      sum(x$death_before_baseline, na.rm = TRUE),
      sum(x$lost_to_followup_before_baseline, na.rm = TRUE),
      sum(x$followup_end_reason == "administrative_censoring", na.rm = TRUE),
      sum(x$followup_end_reason == "death", na.rm = TRUE),
      sum(x$followup_end_reason == "lost_to_followup", na.rm = TRUE)
    )
  )
}
