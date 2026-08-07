# UK Biobank First Occurrence fields are one date field per 3-character ICD-10
# code, with the corresponding source field stored as the next field ID.
.first_occurrence_special_dates <- as.Date(c(
  "1900-01-01", # Code has no event date
  "1901-01-01", # Event date before date of birth
  "1902-02-02", # Event date matches date of birth
  "1903-03-03", # Event date in the same calendar year as birth
  "1909-09-09", # Future placeholder/default
  "2037-07-07"  # Future placeholder/default
))

.first_occurrence_col_candidates <- function(field_id) {
  field_id <- as.character(field_id)
  c(
    paste0("p", field_id),
    paste0("p", field_id, "_i0"),
    paste0("p", field_id, "_i0_a0")
  )
}

.first_occurrence_empty <- function() {
  data.table(
    eid = integer(0),
    fo_field = integer(0),
    fo_source_field = integer(0),
    fo_date = as.Date(character(0)),
    fo_source = character(0)
  )
}

.normalize_first_occurrence_fields <- function(fields, field_name = "first_occurrence_fields") {
  if (is.null(fields)) {
    return(NULL)
  }

  if (is.list(fields) && !is.data.frame(fields)) {
    fields <- unlist(fields, use.names = TRUE)
  }

  if (!is.numeric(fields) && !is.character(fields)) {
    stop(sprintf("'%s' must be an integer or character vector of UKB date field IDs", field_name), call. = FALSE)
  }

  field_ids <- suppressWarnings(as.integer(fields))
  if (any(is.na(field_ids))) {
    stop(sprintf("'%s' contains non-numeric UKB field IDs", field_name), call. = FALSE)
  }
  if (any(field_ids <= 0L)) {
    stop(sprintf("'%s' must contain positive UKB field IDs", field_name), call. = FALSE)
  }

  field_ids <- unique(field_ids)
  names(field_ids) <- NULL
  field_ids
}

.parse_first_occurrence_records <- function(dt, fields, source_fields = NULL) {
  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  fields <- .normalize_first_occurrence_fields(fields)
  if (is.null(fields) || length(fields) == 0) {
    return(.first_occurrence_empty())
  }

  if (!"eid" %in% names(dt)) {
    stop("Column 'eid' is required for FirstOccurrence parsing", call. = FALSE)
  }

  if (is.null(source_fields)) {
    source_fields <- fields + 1L
  } else {
    source_fields <- .normalize_first_occurrence_fields(source_fields, "first_occurrence_source_fields")
    if (length(source_fields) != length(fields)) {
      stop("'first_occurrence_source_fields' must have the same length as 'first_occurrence_fields'", call. = FALSE)
    }
  }

  results <- vector("list", length(fields))
  for (i in seq_along(fields)) {
    field_id <- fields[[i]]
    date_col <- .first_occurrence_col_candidates(field_id)
    date_col <- date_col[date_col %in% names(dt)][1]

    if (is.na(date_col)) {
      next
    }

    source_field <- source_fields[[i]]
    source_col <- .first_occurrence_col_candidates(source_field)
    source_col <- source_col[source_col %in% names(dt)][1]

    if (is.na(source_col)) {
      one <- dt[, .(
        eid,
        fo_field = field_id,
        fo_source_field = source_field,
        fo_date = .safe_as_date(get(date_col), col_name = date_col),
        fo_source = NA_character_
      )]
    } else {
      one <- dt[, .(
        eid,
        fo_field = field_id,
        fo_source_field = source_field,
        fo_date = .safe_as_date(get(date_col), col_name = date_col),
        fo_source = as.character(get(source_col))
      )]
    }

    valid_dates <- !is.na(one[["fo_date"]]) &
      !one[["fo_date"]] %in% .first_occurrence_special_dates
    one <- one[valid_dates]
    if (nrow(one) > 0) {
      results[[i]] <- one
    }
  }

  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) {
    return(.first_occurrence_empty())
  }

  rbindlist(results, use.names = TRUE, fill = TRUE)
}

.aggregate_first_occurrence_earliest <- function(fo_long, disease_label) {
  if (!is.data.table(fo_long)) {
    fo_long <- as.data.table(fo_long)
  }

  if (nrow(fo_long) == 0) {
    return(data.table(
      eid = integer(0),
      disease = character(0),
      earliest_date = as.Date(character(0)),
      source = character(0)
    ))
  }

  fo_long[, disease := disease_label]
  fo_long[
    ,
    {
      dates <- .SD[["fo_date"]]
      sources <- .SD[["fo_source"]]
      min_idx <- which.min(dates)
      raw_source <- sources[min_idx]
      source_label <- if (!is.na(raw_source) && nzchar(trimws(raw_source))) {
        paste0("FirstOccurrence_", trimws(raw_source))
      } else {
        "FirstOccurrence"
      }

      list(
        earliest_date = dates[min_idx],
        source = source_label
      )
    },
    by = .(eid, disease)
  ]
}
