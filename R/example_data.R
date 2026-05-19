#' Load the simulated UK Biobank-style demo dataset
#'
#' @description
#' Returns a fully simulated demo dataset generated for documentation and
#' smoke-test workflows. The packaged object preserves the original column names
#' of a RAP-style cohort table, including `eid`, raw field columns, date
#' columns, and diagnosis-code columns, and includes additional raw-style air
#' pollution fields used by package examples. Values were created by column-wise
#' resampling or independent simulation, so the object does not contain UK
#' Biobank participant-level source records.
#'
#' @param n Optional number of rows to return. If `NULL`, all 5,000 simulated
#'   rows are returned.
#'
#' @return A data.frame of simulated cohort variables with missing values
#'   retained.
#' @export
#'
#' @examples
#' demo <- ukb_demo(100)
#' dim(demo)
#' names(demo)
ukb_demo <- function(n = NULL) {
  out <- ukb_demo_data
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
      stop("`n` must be a positive integer or NULL.", call. = FALSE)
    }
    n <- min(as.integer(n), nrow(out))
    out <- out[seq_len(n), , drop = FALSE]
  }
  message("This dataset is simulated for documentation and testing purposes. It does not contain real UK Biobank participant data.")
  return(out)
}

#' Get column names of the simulated UK Biobank-style demo dataset
#'
#' @description
#' Returns the original column names used by `ukb_demo()`. This is useful for
#' documentation examples that need realistic RAP-style column names without
#' loading the full simulated dataset.
#'
#' @return A character vector of original demo-data column names.
#' @export
#'
#' @examples
#' get_ukb_demo_colnames()
get_ukb_demo_colnames <- function() {
  ukb_demo_colnames
}
