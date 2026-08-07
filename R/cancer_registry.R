.cancer_registry_empty <- function() {
  data.table(
    eid = integer(0),
    cancer_icd10_code = character(0),
    diag_date = as.Date(character(0)),
    cancer_histology = integer(0),
    cancer_behaviour = integer(0),
    source = character(0)
  )
}

.cancer_registry_col <- function(field_id, instance, names_dt) {
  candidates <- c(
    paste0("p", field_id, "_i", instance),
    if (identical(as.integer(instance), 0L)) paste0("p", field_id) else character(0)
  )
  candidates[candidates %in% names_dt][1]
}

.extract_icd10_code <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL", "NAN")] <- NA_character_

  matches <- regexpr("[A-Z][0-9]{2}(\\.?[0-9A-Z]{0,3})?", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  has_match <- !is.na(x) & matches > 0L
  out[has_match] <- regmatches(x, matches)
  out <- gsub("\\.", "", out)
  out
}

.extract_integer_code <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "N/A", "NULL", "NAN")] <- NA_character_
  suppressWarnings(as.integer(sub("^([0-9]+).*$", "\\1", x_chr)))
}

#' @title Parse Cancer Registry Records
#'
#' @description
#' Extracts cancer registry records from UK Biobank fields 40006, 40005, 40011,
#' and 40012 into a standardized long-format table. Field 40006 stores cancer
#' ICD-10 type, 40005 stores diagnosis date, 40011 stores tumour histology, and
#' 40012 stores tumour behaviour.
#'
#' @param dt A data.table or data.frame containing UKB cancer registry columns.
#'
#' @return A data.table with columns: \code{eid}, \code{cancer_icd10_code},
#'   \code{diag_date}, \code{cancer_histology}, \code{cancer_behaviour}, and
#'   \code{source}.
#'
#' @export
parse_cancer_registry <- function(dt) {
  if (!is.data.table(dt)) {
    dt <- as.data.table(dt)
  }

  if (!"eid" %in% names(dt)) {
    stop("Column 'eid' is required for CancerRegistry parsing", call. = FALSE)
  }

  code_cols <- grep("^p40006_i[0-9]+$", names(dt), value = TRUE)
  if ("p40006" %in% names(dt)) {
    code_cols <- unique(c(code_cols, "p40006"))
  }

  if (length(code_cols) == 0) {
    message("[parse_cancer_registry] No p40006 cancer ICD-10 columns found")
    return(.cancer_registry_empty())
  }

  results <- vector("list", length(code_cols))
  names_dt <- names(dt)

  for (i in seq_along(code_cols)) {
    code_col <- code_cols[[i]]
    instance <- if (grepl("^p40006_i[0-9]+$", code_col)) {
      as.integer(sub("^p40006_i", "", code_col))
    } else {
      0L
    }

    date_col <- .cancer_registry_col(40005, instance, names_dt)
    hist_col <- .cancer_registry_col(40011, instance, names_dt)
    beh_col <- .cancer_registry_col(40012, instance, names_dt)

    code_vec <- .extract_icd10_code(dt[[code_col]])
    date_vec <- if (!is.na(date_col)) {
      .safe_as_date(dt[[date_col]], col_name = date_col)
    } else {
      as.Date(rep(NA_character_, nrow(dt)))
    }
    hist_vec <- if (!is.na(hist_col)) .extract_integer_code(dt[[hist_col]]) else rep(NA_integer_, nrow(dt))
    beh_vec <- if (!is.na(beh_col)) .extract_integer_code(dt[[beh_col]]) else rep(NA_integer_, nrow(dt))

    one <- data.table(
      eid = dt[["eid"]],
      cancer_icd10_code = code_vec,
      diag_date = date_vec,
      cancer_histology = hist_vec,
      cancer_behaviour = beh_vec,
      source = "CancerRegistry"
    )
    valid_codes <- !is.na(one[["cancer_icd10_code"]]) & nzchar(one[["cancer_icd10_code"]])
    one <- one[valid_codes]
    if (nrow(one) > 0) {
      results[[i]] <- one
    }
  }

  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) {
    return(.cancer_registry_empty())
  }

  rbindlist(results, use.names = TRUE, fill = TRUE)
}

#' @title Filter Cancer Registry Records by ICD-10 and Tumour Metadata
#'
#' @description
#' Filters cancer registry records using ICD-10 pattern matching and optional
#' tumour histology / behaviour constraints.
#'
#' @param cancer_long A data.table from \code{\link{parse_cancer_registry}}.
#' @param pattern Regular expression pattern for cancer ICD-10 codes.
#' @param disease_label Disease name label to assign to matched records.
#' @param histology Optional integer vector of ICD-O histology codes.
#' @param behaviour Optional integer vector of ICD-O behaviour codes. Use
#'   \code{3L} for malignant tumours.
#'
#' @return A filtered data.table with an added \code{disease} column.
#'
#' @keywords internal
filter_cancer_registry <- function(cancer_long,
                                   pattern,
                                   disease_label,
                                   histology = NULL,
                                   behaviour = NULL) {
  if (!is.data.table(cancer_long)) {
    cancer_long <- as.data.table(cancer_long)
  }

  result <- cancer_long[grepl(pattern, cancer_long[["cancer_icd10_code"]], perl = TRUE)]

  if (!is.null(histology)) {
    histology <- .extract_integer_code(histology)
    result <- result[result[["cancer_histology"]] %in% histology]
  }

  if (!is.null(behaviour)) {
    behaviour <- .extract_integer_code(behaviour)
    result <- result[result[["cancer_behaviour"]] %in% behaviour]
  }

  set(result, j = "disease", value = disease_label)
  result
}

#' @title Aggregate Earliest Cancer Registry Diagnosis Date
#'
#' @description
#' Computes the earliest cancer registry diagnosis date for each
#' participant-disease combination.
#'
#' @param cancer_filtered A data.table from
#'   \code{\link{filter_cancer_registry}}.
#'
#' @return A data.table with columns: \code{eid}, \code{disease},
#'   \code{earliest_date}, \code{source}.
#'
#' @keywords internal
aggregate_cancer_registry_earliest <- function(cancer_filtered) {
  if (!is.data.table(cancer_filtered)) {
    cancer_filtered <- as.data.table(cancer_filtered)
  }

  cancer_filtered <- cancer_filtered[!is.na(cancer_filtered[["diag_date"]])]
  cancer_filtered[
    ,
    {
      dates <- .SD[["diag_date"]]
      list(earliest_date = min(dates), source = "CancerRegistry")
    },
    by = c("eid", "disease")
  ]
}
