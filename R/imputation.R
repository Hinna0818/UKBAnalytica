## Imputation helpers
#' Multiple imputation and merge back to full data
#'
#' Run multiple imputation with the CRAN package \pkg{mice} on a subset of variables,
#' then merge the imputed columns back to the original dataset by an ID column.
#'
#' This function is designed for workflows where you want to keep a set of
#' "static" columns (exposures, outcomes, follow-up time, etc.) untouched while
#' imputing a selected set of covariates.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Subsets the input data to the requested variables.
#'   \item Runs \code{mice::mice()}.
#'   \item Creates \code{m} completed datasets and merges imputed columns back.
#'   \item Optionally merges additional datasets (e.g., omics) by ID.
#' }
#'
#' Factor handling: for variables listed in \code{factor_vars}, the function will
#' coerce them to factors before imputation. All other variables in \code{vars}
#' are coerced to numeric.
#'
#' @references \url{https://github.com/amices/mice}
#'
#' @param data A data.frame/data.table containing the cohort.
#' @param id_col Name of the ID column. Default is \code{"eid"}.
#' @param vars Character vector of column names to impute.
#' @param factor_vars Optional character vector of variables (subset of
#'   \code{vars}) to treat as categorical (factors).
#' @param method Imputation method passed to \code{mice::mice()}. Default is
#'   \code{"pmm"}.
#' @param m Number of multiple imputations. Default is 5.
#' @param maxit Maximum number of iterations. Default is 10.
#' @param seed Random seed for reproducibility.
#' @param print Logical. If TRUE, show mice iteration logs.
#' @param additional_data Optional named list of extra datasets to merge after
#'   imputation. Each element must contain \code{id_col}. Example:
#'   \code{list(protein = protein_df, metabolomics = meta_df)}.
#' @param additional_join Join type for additional datasets. One of
#'   \code{"inner"} or \code{"left"}. Default is \code{"inner"}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{imp}: the \code{mice} \code{mids} object
#'   \item \code{data_list}: a list of length \code{m} containing completed and
#'     merged datasets
#' }
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' dt <- data.table(eid = 1:100, x = rnorm(100), y = rnorm(100))
#' dt[sample.int(100, 10), x := NA]
#'
#' res <- impute_mice_merge(dt, vars = c("x", "y"), m = 2, maxit = 2)
#' imputed_dt1 <- res$data_list[[1]]
#' }
#'
#' @export
run_imputation <- function(data,
                              id_col = "eid",
                              vars,
                              factor_vars = NULL,
                              method = "pmm",
                              m = 5,
                              maxit = 10,
                              seed = 1234,
                              print = TRUE,
                              additional_data = NULL,
                              additional_join = c("inner", "left")) {
  additional_join <- match.arg(additional_join)

  if (!requireNamespace("mice", quietly = TRUE)) {
    stop(
      "Package 'mice' is required but not installed. ",
      "Install it from CRAN with install.packages('mice')."
    )
  }

  if (missing(vars) || length(vars) == 0) {
    stop("`vars` must be a non-empty character vector.")
  }
  if (!is.character(vars)) {
    stop("`vars` must be a character vector.")
  }

  if (!id_col %in% names(data)) {
    stop(sprintf("ID column '%s' not found in `data`.", id_col))
  }

  vars <- unique(vars)
  vars_present <- vars[vars %in% names(data)]
  missing_vars <- setdiff(vars, vars_present)
  if (length(vars_present) == 0) {
    stop("None of the requested `vars` exist in `data`.")
  }
  if (length(missing_vars) > 0) {
    warning(
      "The following variables were not found in `data` and will be ignored: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # Prepare skeleton data with all non-imputed columns.
  static_data <- data
  static_data[vars_present] <- NULL

  # Build imputation input.
  imp_df <- data.frame(data[[id_col]], stringsAsFactors = FALSE)
  names(imp_df) <- id_col
  for (v in vars_present) {
    imp_df[[v]] <- data[[v]]
    if (is.character(imp_df[[v]])) {
      imp_df[[v]][imp_df[[v]] %in% c("")] <- NA
    }
  }

  factor_vars <- intersect(factor_vars %||% character(0), vars_present)
  num_vars <- setdiff(vars_present, factor_vars)
  for (v in factor_vars) imp_df[[v]] <- as.factor(imp_df[[v]])
  for (v in num_vars) imp_df[[v]] <- suppressWarnings(as.numeric(imp_df[[v]]))

  imp_input <- imp_df
  imp_input[[id_col]] <- NULL

  imp <- mice::mice(
    data = imp_input,
    m = m,
    method = method,
    maxit = maxit,
    seed = seed,
    printFlag = isTRUE(print)
  )

  eid_vec <- imp_df[[id_col]]

  out_list <- vector("list", length = imp$m)
  for (i in seq_len(imp$m)) {
    completed <- mice::complete(imp, i)
    completed[[id_col]] <- eid_vec

    # Merge imputed covariates back to static skeleton.
    merged <- merge(
      x = static_data,
      y = completed,
      by = id_col,
      all.x = TRUE,
      sort = FALSE
    )

    # Optionally merge extra datasets.
    if (!is.null(additional_data)) {
      if (!is.list(additional_data) || is.null(names(additional_data))) {
        stop("`additional_data` must be a named list of data.frames/data.tables.")
      }
      for (nm in names(additional_data)) {
        extra <- additional_data[[nm]]
        if (!id_col %in% names(extra)) {
          stop(sprintf("additional_data$%s does not contain '%s'.", nm, id_col))
        }
        merged <- merge(
          x = merged,
          y = extra,
          by = id_col,
          all = (additional_join == "left"),
          all.x = (additional_join == "left"),
          sort = FALSE
        )
      }
    }

    out_list[[i]] <- merged
  }

  return(list(imp = imp, data_list = out_list))
}


#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x