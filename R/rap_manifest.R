#' Check the UK Biobank RAP execution environment
#'
#' @description
#' Inspect whether the current R session is running inside a UK Biobank
#' Research Analysis Platform (RAP)-like environment and return reproducible
#' diagnostics for RAP-aware workflows. The function only checks environment
#' variables, local paths, and the availability of the `dx` command-line tool;
#' it does not read or export participant-level data.
#'
#' @param output_dir Optional output directory to assess.
#' @param require_rap Logical. If `TRUE`, mark the check as failed when the
#'   session does not appear to be running on RAP.
#' @param require_dx Logical. If `TRUE`, mark the check as failed when the
#'   `dx` command-line tool is unavailable.
#' @param check_write Logical. If `TRUE` and `output_dir` is provided, test
#'   whether a small temporary file can be written and removed.
#' @param verbose Logical. If `TRUE`, print a compact summary.
#'
#' @return A list with class `ukb_rap_env` containing RAP environment metadata
#'   and a check table.
#' @export
#'
#' @examples
#' env <- ukb_check_rap_env(verbose = FALSE)
ukb_check_rap_env <- function(output_dir = NULL,
                              require_rap = FALSE,
                              require_dx = FALSE,
                              check_write = FALSE,
                              verbose = TRUE) {
  if (!is.logical(require_rap) || length(require_rap) != 1 || is.na(require_rap)) {
    stop("`require_rap` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(require_dx) || length(require_dx) != 1 || is.na(require_dx)) {
    stop("`require_dx` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(check_write) || length(check_write) != 1 || is.na(check_write)) {
    stop("`check_write` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  env_vars <- Sys.getenv(
    c("DX_PROJECT_CONTEXT_ID", "DX_WORKSPACE_ID", "DX_JOB_ID"),
    unset = NA_character_
  )
  dx_available <- nzchar(Sys.which("dx"))
  rap_path_exists <- dir.exists("/mnt/project")
  is_rap <- any(!is.na(env_vars) & nzchar(env_vars)) || rap_path_exists

  checks <- data.frame(
    check = character(),
    status = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
  add_check <- function(check, ok, message) {
    checks <<- rbind(
      checks,
      data.frame(
        check = check,
        status = if (isTRUE(ok)) "pass" else "fail",
        message = message,
        stringsAsFactors = FALSE
      )
    )
  }

  add_check(
    "RAP environment",
    is_rap || !isTRUE(require_rap),
    if (is_rap) "RAP-like environment detected." else "RAP environment variables or /mnt/project were not detected."
  )
  add_check(
    "dx command",
    dx_available || !isTRUE(require_dx),
    if (dx_available) "dx command-line tool is available." else "dx command-line tool is not available on PATH."
  )

  output_info <- NULL
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1 || is.na(output_dir) || !nzchar(output_dir)) {
      stop("`output_dir` must be a single non-empty path.", call. = FALSE)
    }
    output_dir <- path.expand(output_dir)
    output_exists <- dir.exists(output_dir)
    output_safe <- .ukb_output_path_is_rap_scoped(output_dir)
    output_writable <- NA
    if (output_exists && isTRUE(check_write)) {
      output_writable <- .ukb_test_write_access(output_dir)
    }
    add_check(
      "output directory",
      output_exists,
      if (output_exists) "Output directory exists." else "Output directory does not exist."
    )
    add_check(
      "RAP-scoped output path",
      output_safe,
      if (output_safe) "Output path is under /mnt/project or a DNAnexus project URI." else "Output path is not RAP-scoped."
    )
    if (isTRUE(check_write)) {
      add_check(
        "output write access",
        isTRUE(output_writable),
        if (isTRUE(output_writable)) "Output directory is writable." else "Output directory was not writable."
      )
    }
    output_info <- list(
      path = output_dir,
      exists = output_exists,
      rap_scoped = output_safe,
      writable = output_writable
    )
  }

  out <- list(
    is_rap = is_rap,
    dx_available = dx_available,
    project_context_id = unname(env_vars[["DX_PROJECT_CONTEXT_ID"]]),
    workspace_id = unname(env_vars[["DX_WORKSPACE_ID"]]),
    job_id = unname(env_vars[["DX_JOB_ID"]]),
    output = output_info,
    checks = checks
  )
  class(out) <- c("ukb_rap_env", class(out))

  if (isTRUE(verbose)) {
    print(out)
  }
  out
}

#' @export
print.ukb_rap_env <- function(x, ...) {
  cat("UKB RAP environment check\n")
  cat("  RAP detected: ", if (isTRUE(x$is_rap)) "yes" else "no", "\n", sep = "")
  cat("  dx available: ", if (isTRUE(x$dx_available)) "yes" else "no", "\n", sep = "")
  if (!is.null(x$output)) {
    cat("  output path: ", x$output$path, "\n", sep = "")
  }
  print(x$checks, row.names = FALSE)
  invisible(x)
}

#' Create a RAP extraction manifest
#'
#' @description
#' Build a compact manifest describing the UKB fields intended for RAP
#' extraction. This is designed as an auditable planning object that can be
#' stored with analysis scripts before running `rap_plan_extract()` or
#' `rap_extract_pheno()`.
#'
#' @param field_id Optional numeric or character vector of UKB field IDs.
#' @param variable_set Optional curated variable-set name from
#'   [get_variable_sets()].
#' @param variables Optional predefined variable names from
#'   [get_variable_info()].
#' @param dataset Optional RAP dataset name.
#' @param entity RAP entity name, usually `"participant"`.
#' @param output Optional planned extraction output path.
#' @param include_eid Logical. Whether participant ID is expected in the
#'   extraction.
#' @param purpose Optional short description of the analysis purpose.
#' @param notes Optional free-text notes.
#'
#' @return A list with class `ukb_extraction_manifest`.
#' @export
#'
#' @examples
#' manifest <- ukb_create_extraction_manifest(
#'   field_id = c(31, 21022),
#'   variable_set = "clinical_core",
#'   purpose = "demo"
#' )
ukb_create_extraction_manifest <- function(field_id = NULL,
                                           variable_set = NULL,
                                           variables = NULL,
                                           dataset = NULL,
                                           entity = "participant",
                                           output = NULL,
                                           include_eid = TRUE,
                                           purpose = NULL,
                                           notes = NULL) {
  if (is.null(field_id) && is.null(variable_set) && is.null(variables)) {
    stop("Provide at least one of `field_id`, `variable_set`, or `variables`.", call. = FALSE)
  }
  if (!is.null(entity) && (!is.character(entity) || length(entity) != 1 || is.na(entity))) {
    stop("`entity` must be a single character string.", call. = FALSE)
  }
  if (!is.logical(include_eid) || length(include_eid) != 1 || is.na(include_eid)) {
    stop("`include_eid` must be TRUE or FALSE.", call. = FALSE)
  }

  rows <- list()

  if (!is.null(field_id)) {
    ids <- unique(as.integer(field_id))
    if (anyNA(ids)) {
      stop("`field_id` must contain numeric field IDs.", call. = FALSE)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      source = "field_id",
      set = NA_character_,
      variable = NA_character_,
      field_id = ids,
      ukb_col = paste0("p", ids),
      label = NA_character_,
      role = "requested_field",
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(variable_set)) {
    set_rows <- get_variable_sets(set = variable_set)
    rows[[length(rows) + 1L]] <- data.frame(
      source = "variable_set",
      set = set_rows$set,
      variable = set_rows$variable,
      field_id = set_rows$field_id,
      ukb_col = set_rows$ukb_col,
      label = set_rows$label,
      role = set_rows$role,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(variables)) {
    var_rows <- get_variable_info()
    missing_vars <- setdiff(variables, var_rows$variable)
    if (length(missing_vars) > 0) {
      stop(
        "Unknown predefined variable(s): ",
        paste(missing_vars, collapse = ", "),
        ". Use get_variable_info() to inspect available variables.",
        call. = FALSE
      )
    }
    var_rows <- var_rows[var_rows$variable %in% variables, , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      source = "predefined_variable",
      set = NA_character_,
      variable = var_rows$variable,
      field_id = as.integer(var_rows$field_id),
      ukb_col = var_rows$ukb_column,
      label = var_rows$description,
      role = "predefined_variable",
      stringsAsFactors = FALSE
    )
  }

  fields <- do.call(rbind, rows)
  fields <- fields[!duplicated(paste(fields$field_id, fields$ukb_col, fields$variable, sep = "::")), , drop = FALSE]
  rownames(fields) <- NULL

  env <- ukb_check_rap_env(verbose = FALSE)
  manifest <- list(
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    package_version = as.character(utils::packageVersion("UKBAnalytica")),
    dataset = dataset,
    entity = entity,
    output = output,
    include_eid = include_eid,
    purpose = purpose,
    notes = notes,
    rap = list(
      is_rap = env$is_rap,
      dx_available = env$dx_available,
      project_context_id = env$project_context_id,
      workspace_id = env$workspace_id,
      job_id = env$job_id
    ),
    fields = fields
  )
  class(manifest) <- c("ukb_extraction_manifest", class(manifest))
  manifest
}

#' @export
print.ukb_extraction_manifest <- function(x, ...) {
  cat("UKB RAP extraction manifest\n")
  cat("  created_at: ", x$created_at, "\n", sep = "")
  cat("  dataset: ", if (is.null(x$dataset)) "<not set>" else x$dataset, "\n", sep = "")
  cat("  entity: ", x$entity, "\n", sep = "")
  cat("  include_eid: ", if (isTRUE(x$include_eid)) "yes" else "no", "\n", sep = "")
  cat("  fields: ", nrow(x$fields), " rows, ", length(unique(x$fields$field_id)), " unique field IDs\n", sep = "")
  invisible(x)
}

#' Write a RAP extraction manifest
#'
#' @param manifest A `ukb_extraction_manifest` object.
#' @param path Output path.
#' @param format Output format: `"csv"` writes the field table and a sidecar
#'   summary CSV, while `"rds"` writes the full manifest object.
#'
#' @return The output path, invisibly.
#' @export
#'
#' @examples
#' manifest <- ukb_create_extraction_manifest(field_id = c(31, 21022))
#' tmp <- tempfile(fileext = ".csv")
#' ukb_write_extraction_manifest(manifest, tmp)
ukb_write_extraction_manifest <- function(manifest,
                                          path,
                                          format = c("csv", "rds")) {
  format <- match.arg(format)
  if (!inherits(manifest, "ukb_extraction_manifest")) {
    stop("`manifest` must be created by ukb_create_extraction_manifest().", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty file path.", call. = FALSE)
  }

  if (identical(format, "rds")) {
    saveRDS(manifest, path)
    return(invisible(path))
  }

  utils::write.csv(manifest$fields, path, row.names = FALSE)
  summary_path <- paste0(tools::file_path_sans_ext(path), "_summary.csv")
  summary_df <- data.frame(
    item = c("created_at", "package_version", "dataset", "entity", "output", "include_eid", "purpose", "notes"),
    value = c(
      manifest$created_at,
      manifest$package_version,
      .ukb_null_to_na(manifest$dataset),
      manifest$entity,
      .ukb_null_to_na(manifest$output),
      as.character(manifest$include_eid),
      .ukb_null_to_na(manifest$purpose),
      .ukb_null_to_na(manifest$notes)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary_df, summary_path, row.names = FALSE)
  invisible(path)
}

.ukb_null_to_na <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  as.character(x)
}

.ukb_output_path_is_rap_scoped <- function(path) {
  path <- path.expand(path)
  startsWith(path, "/mnt/project") || startsWith(path, "project-") || startsWith(path, "dx://")
}

.ukb_test_write_access <- function(path) {
  tf <- tempfile(pattern = ".ukbanalytica-write-test-", tmpdir = path)
  ok <- tryCatch({
    writeLines("ok", tf)
    file.exists(tf)
  }, error = function(e) FALSE)
  if (file.exists(tf)) {
    unlink(tf)
  }
  isTRUE(ok)
}
