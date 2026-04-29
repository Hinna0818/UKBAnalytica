.ukb_snapshot_store <- new.env(parent = emptyenv())

.ukb_missing_col <- function(x) {
  if (is.character(x) || is.factor(x)) {
    x_chr <- trimws(as.character(x))
    is.na(x) | x_chr == ""
  } else {
    is.na(x)
  }
}

#' @title Clean UK Biobank Missing and Non-response Values
#'
#' @description
#' Converts common UK Biobank non-response labels and numeric missing codes into
#' analysis-ready missing values. Empty strings are always converted to
#' \code{NA}. Informative non-response labels can either be converted to
#' \code{NA} or retained as \code{"Unknown"} for modelling.
#'
#' @param data A data.frame or data.table.
#' @param cols Optional character vector of columns to clean. If NULL, all
#'   columns are considered.
#' @param action How to handle informative character labels:
#'   \code{"na"} converts them to \code{NA}; \code{"unknown"} converts them to
#'   \code{"Unknown"}. Numeric missing codes are always converted to \code{NA}.
#' @param extra_labels Additional character labels to treat as informative
#'   missing.
#' @param numeric_codes Numeric values to treat as missing. Defaults to common
#'   UKB values \code{-1} and \code{-3}.
#' @param trim Logical. Trim leading/trailing whitespace in character columns.
#' @param in_place Logical. If TRUE and \code{data} is a data.table, modify by
#'   reference. Default FALSE returns a cleaned copy.
#' @param verbose Logical. Print a concise cleaning summary.
#'
#' @return A data.table.
#'
#' @export
ukb_clean_missing <- function(data,
                              cols = NULL,
                              action = c("na", "unknown"),
                              extra_labels = NULL,
                              numeric_codes = c(-1, -3),
                              trim = TRUE,
                              in_place = FALSE,
                              verbose = TRUE) {
  action <- match.arg(action)

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table", call. = FALSE)
  }

  if (is.null(cols)) {
    cols <- names(data)
  }
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  informative_labels <- c(
    "Do not know",
    "Prefer not to answer",
    "Prefer not to say",
    extra_labels
  )
  informative_labels <- unique(informative_labels[!is.na(informative_labels) & nzchar(informative_labels)])

  if (data.table::is.data.table(data) && isTRUE(in_place)) {
    dt <- data
  } else {
    dt <- data.table::as.data.table(data.table::copy(data))
  }

  n_cols_affected <- 0L
  n_values_replaced <- 0L
  per_col <- data.table::data.table(
    column = character(0),
    empty_to_na = integer(0),
    label_replaced = integer(0),
    numeric_code_to_na = integer(0)
  )

  for (col in cols) {
    x <- dt[[col]]
    empty_n <- 0L
    label_n <- 0L
    numeric_n <- 0L

    if (is.factor(x)) {
      x <- as.character(x)
    }

    if (is.character(x)) {
      x_chr <- if (isTRUE(trim)) trimws(x) else x
      empty_idx <- !is.na(x_chr) & x_chr == ""
      label_idx <- !is.na(x_chr) & x_chr %in% informative_labels

      empty_n <- sum(empty_idx)
      label_n <- sum(label_idx)

      if (empty_n > 0L) {
        x_chr[empty_idx] <- NA_character_
      }
      if (label_n > 0L) {
        x_chr[label_idx] <- if (action == "unknown") "Unknown" else NA_character_
      }
      if (empty_n + label_n > 0L || is.factor(dt[[col]]) || isTRUE(trim)) {
        data.table::set(dt, j = col, value = x_chr)
      }
    } else if (is.numeric(x) || is.integer(x)) {
      numeric_idx <- !is.na(x) & x %in% numeric_codes
      numeric_n <- sum(numeric_idx)
      if (numeric_n > 0L) {
        x[numeric_idx] <- NA
        data.table::set(dt, j = col, value = x)
      }
    }

    total_n <- empty_n + label_n + numeric_n
    if (total_n > 0L) {
      n_cols_affected <- n_cols_affected + 1L
      n_values_replaced <- n_values_replaced + total_n
      per_col <- data.table::rbindlist(
        list(
          per_col,
          data.table::data.table(
            column = col,
            empty_to_na = empty_n,
            label_replaced = label_n,
            numeric_code_to_na = numeric_n
          )
        ),
        use.names = TRUE
      )
    }
  }

  attr(dt, "ukb_clean_missing_summary") <- per_col

  if (isTRUE(verbose)) {
    message(sprintf(
      "[ukb_clean_missing] Replaced %d value(s) across %d column(s) (action = %s)",
      n_values_replaced,
      n_cols_affected,
      action
    ))
  }

  dt
}

#' @title Record or Retrieve UKB Cohort Snapshots
#'
#' @description
#' Records lightweight cohort checkpoints during an analysis pipeline. Each
#' snapshot stores row count, column count, number of columns containing missing
#' values, complete row count, object size, and deltas from the previous
#' snapshot. Calling \code{ukb_snapshot()} without \code{data} returns the
#' current snapshot history.
#'
#' @param data Optional data.frame or data.table. If supplied, records a new
#'   snapshot.
#' @param label Snapshot label. Required when recording a new snapshot.
#' @param id Snapshot stream identifier. Use separate IDs for independent
#'   pipelines in the same R session.
#' @param reset Logical. If TRUE, clears the snapshot history for \code{id}.
#' @param verbose Logical. Print a concise snapshot summary.
#'
#' @return A data.table snapshot history.
#'
#' @export
ukb_snapshot <- function(data = NULL,
                         label = NULL,
                         id = "default",
                         reset = FALSE,
                         verbose = TRUE) {
  if (!is.character(id) || length(id) != 1L || !nzchar(id)) {
    stop("'id' must be a non-empty character scalar", call. = FALSE)
  }

  if (isTRUE(reset)) {
    if (exists(id, envir = .ukb_snapshot_store, inherits = FALSE)) {
      rm(list = id, envir = .ukb_snapshot_store)
    }
    empty <- data.table::data.table()
    if (isTRUE(verbose)) {
      message(sprintf("[ukb_snapshot] Cleared snapshot history: %s", id))
    }
    return(empty)
  }

  if (is.null(data)) {
    if (exists(id, envir = .ukb_snapshot_store, inherits = FALSE)) {
      return(data.table::copy(get(id, envir = .ukb_snapshot_store, inherits = FALSE)))
    }
    return(data.table::data.table())
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table", call. = FALSE)
  }
  if (is.null(label) || !is.character(label) || length(label) != 1L || !nzchar(label)) {
    stop("'label' must be a non-empty character scalar when recording a snapshot", call. = FALSE)
  }

  dt <- data.table::as.data.table(data)
  missing_cols <- vapply(dt, function(x) any(.ukb_missing_col(x)), logical(1))
  complete_rows <- if (ncol(dt) == 0L) {
    nrow(dt)
  } else {
    row_missing <- Reduce(`|`, lapply(dt, .ukb_missing_col))
    sum(!row_missing)
  }

  history <- if (exists(id, envir = .ukb_snapshot_store, inherits = FALSE)) {
    get(id, envir = .ukb_snapshot_store, inherits = FALSE)
  } else {
    data.table::data.table()
  }

  size_mb <- as.numeric(utils::object.size(dt)) / 1024^2
  idx <- nrow(history) + 1L

  previous <- if (nrow(history) > 0L) history[nrow(history)] else NULL
  row_delta <- if (is.null(previous)) NA_integer_ else nrow(dt) - previous[["nrow"]]
  col_delta <- if (is.null(previous)) NA_integer_ else ncol(dt) - previous[["ncol"]]
  size_delta_mb <- if (is.null(previous)) NA_real_ else size_mb - previous[["size_mb"]]

  snapshot <- data.table::data.table(
    idx = idx,
    label = label,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    nrow = nrow(dt),
    ncol = ncol(dt),
    n_missing_cols = sum(missing_cols),
    n_complete_rows = complete_rows,
    size_mb = round(size_mb, 3),
    row_delta = row_delta,
    col_delta = col_delta,
    size_delta_mb = round(size_delta_mb, 3)
  )

  history <- data.table::rbindlist(list(history, snapshot), use.names = TRUE, fill = TRUE)
  assign(id, history, envir = .ukb_snapshot_store)

  if (isTRUE(verbose)) {
    message(sprintf(
      "[ukb_snapshot] %s: rows=%d%s, cols=%d%s, missing_cols=%d, complete_rows=%d",
      label,
      nrow(dt),
      ifelse(is.na(row_delta), "", sprintf(" (%+d)", row_delta)),
      ncol(dt),
      ifelse(is.na(col_delta), "", sprintf(" (%+d)", col_delta)),
      sum(missing_cols),
      complete_rows
    ))
  }

  data.table::copy(history)
}
