#' RAP Phenotype Extraction Helpers
#'
#' @description
#' R-native wrappers around DNAnexus `dx extract_dataset` and the RAP
#' `table-exporter` app. These functions are intended for use inside approved
#' UK Biobank RAP sessions or RAP-controlled execution environments.
#'
#' @name rap_extract
#' @keywords internal
NULL

.rap_extract_cache <- new.env(parent = emptyenv())

.rap_check_logical <- function(x, name) {
  if (!isTRUE(x) && !identical(x, FALSE)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
}

.rap_require_dx <- function() {
  if (!nzchar(Sys.which("dx"))) {
    stop(
      "DNAnexus CLI 'dx' was not found on PATH. Run this inside RAP or install dx-toolkit.",
      call. = FALSE
    )
  }
}

.rap_dx_env <- function(env = character()) {
  env <- c(PYTHONIOENCODING = "utf-8", env)
  paste0(names(env), "=", unname(env))
}

.rap_dx_run <- function(args,
                        timeout = 300,
                        env = character()) {
  .rap_require_dx()

  stdout_file <- tempfile("ukba_dx_stdout_")
  stderr_file <- tempfile("ukba_dx_stderr_")
  on.exit(unlink(c(stdout_file, stderr_file)), add = TRUE)

  status <- system2(
    command = "dx",
    args = args,
    stdout = stdout_file,
    stderr = stderr_file,
    env = .rap_dx_env(env),
    timeout = timeout
  )

  stdout <- if (file.exists(stdout_file)) paste(readLines(stdout_file, warn = FALSE), collapse = "\n") else ""
  stderr <- if (file.exists(stderr_file)) paste(readLines(stderr_file, warn = FALSE), collapse = "\n") else ""
  if (is.null(status)) {
    status <- 0L
  }

  list(
    stdout = stdout,
    stderr = stderr,
    status = as.integer(status),
    success = identical(as.integer(status), 0L)
  )
}

.rap_write_lf_file <- function(lines, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n")
  invisible(path)
}

.rap_parse_fields <- function(stdout) {
  lines <- strsplit(stdout, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines[nzchar(lines)])

  if (length(lines) == 0) {
    return(data.frame(
      field_name = character(0),
      title = character(0),
      stringsAsFactors = FALSE
    ))
  }

  parts <- strsplit(lines, "\t", fixed = TRUE)
  data.frame(
    field_name = vapply(parts, function(x) x[[1]], character(1)),
    title = vapply(
      parts,
      function(x) if (length(x) >= 2) x[[2]] else NA_character_,
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

.rap_assert_fields_df <- function(fields_df) {
  if (!is.data.frame(fields_df) || !all(c("field_name", "title") %in% names(fields_df))) {
    stop("'fields_df' must be a data.frame with columns 'field_name' and 'title'.", call. = FALSE)
  }
  invisible(fields_df)
}

.rap_normalize_field_ids <- function(field_id) {
  if (is.null(field_id)) {
    return(character(0))
  }

  field_id <- as.character(field_id)
  field_id <- trimws(field_id[nzchar(field_id)])
  if (length(field_id) == 0 || anyNA(field_id) || any(!grepl("^[0-9]+$", field_id))) {
    stop("'field_id' must contain UKB numeric field IDs, e.g. c(31, 53, 21022).", call. = FALSE)
  }
  unique(field_id)
}

.rap_variable_field_names <- function(variables) {
  if (is.null(variables)) {
    return(character(0))
  }
  if (!is.character(variables) || length(variables) == 0 || anyNA(variables)) {
    stop("'variables' must be a non-empty character vector.", call. = FALSE)
  }

  mapping <- .get_variable_mapping()
  missing_vars <- setdiff(variables, names(mapping))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "Unknown predefined variable(s): %s. Use get_variable_info() to inspect available names.",
      paste(missing_vars, collapse = ", ")
    ), call. = FALSE)
  }

  unique(vapply(mapping[variables], function(x) as.character(x$ukb_col), character(1)))
}

.rap_normalize_field_names <- function(field_names,
                                       entity = "participant") {
  if (is.null(field_names) || length(field_names) == 0) {
    return(character(0))
  }
  if (!is.character(field_names) || length(field_names) == 0 || anyNA(field_names)) {
    stop("'field_names' must be a non-empty character vector.", call. = FALSE)
  }

  field_names <- trimws(field_names[nzchar(field_names)])
  field_names <- ifelse(
    grepl("\\.", field_names, fixed = FALSE),
    field_names,
    paste0(entity, ".", field_names)
  )
  unique(field_names)
}

.rap_match_field_ids <- function(field_id,
                                 fields_df,
                                 entity = "participant") {
  matched <- vector("list", length(field_id))
  unmatched <- character(0)
  out_idx <- 0L

  for (fid in field_id) {
    pattern <- paste0("^", entity, "\\.p", fid, "(_|$)")
    hits <- fields_df[grepl(pattern, fields_df$field_name, perl = TRUE), , drop = FALSE]

    if (nrow(hits) == 0) {
      unmatched <- c(unmatched, fid)
      next
    }

    out_idx <- out_idx + 1L
    base_title <- sub("\\s*\\|.*$", "", hits$title[[1]])
    matched[[out_idx]] <- data.frame(
      request_type = "field_id",
      request = fid,
      field_id = fid,
      title = trimws(base_title),
      n_cols = nrow(hits),
      field_names = paste(hits$field_name, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  matched <- matched[seq_len(out_idx)]
  list(
    matched = if (length(matched) > 0) do.call(rbind, matched) else data.frame(
      request_type = character(0),
      request = character(0),
      field_id = character(0),
      title = character(0),
      n_cols = integer(0),
      field_names = character(0),
      stringsAsFactors = FALSE
    ),
    unmatched = unmatched
  )
}

.rap_match_direct_fields <- function(field_names,
                                     fields_df) {
  if (length(field_names) == 0) {
    return(list(
      matched = data.frame(
        request_type = character(0),
        request = character(0),
        field_id = character(0),
        title = character(0),
        n_cols = integer(0),
        field_names = character(0),
        stringsAsFactors = FALSE
      ),
      unmatched = character(0)
    ))
  }

  hits <- fields_df[match(field_names, fields_df$field_name), , drop = FALSE]
  found <- !is.na(hits$field_name)
  matched <- data.frame(
    request_type = "field_name",
    request = field_names[found],
    field_id = sub("^.*\\.p([0-9]+).*$", "\\1", field_names[found]),
    title = hits$title[found],
    n_cols = 1L,
    field_names = field_names[found],
    stringsAsFactors = FALSE
  )

  list(
    matched = matched,
    unmatched = field_names[!found]
  )
}

.rap_plan_fields <- function(field_id = NULL,
                             field_names = NULL,
                             variables = NULL,
                             fields_df,
                             entity = "participant",
                             include_eid = TRUE,
                             table_exporter = FALSE) {
  .rap_assert_fields_df(fields_df)
  .rap_check_logical(include_eid, "include_eid")
  .rap_check_logical(table_exporter, "table_exporter")

  field_id <- unique(c(
    .rap_normalize_field_ids(field_id)
  ))
  field_names <- unique(c(
    field_names,
    .rap_variable_field_names(variables)
  ))
  field_names <- .rap_normalize_field_names(field_names, entity = entity)

  id_match <- .rap_match_field_ids(field_id, fields_df, entity = entity)
  direct_match <- .rap_match_direct_fields(field_names, fields_df)

  matched <- rbind(id_match$matched, direct_match$matched)
  unmatched <- c(id_match$unmatched, direct_match$unmatched)

  extracted_fields <- character(0)
  if (nrow(matched) > 0) {
    extracted_fields <- unique(unlist(strsplit(matched$field_names, ";", fixed = TRUE), use.names = FALSE))
  }

  eid_field <- paste0(entity, ".eid")
  if (isTRUE(include_eid)) {
    extracted_fields <- unique(c(eid_field, extracted_fields))
  }

  if (isTRUE(table_exporter)) {
    extracted_fields <- sub(paste0("^", entity, "\\."), "", extracted_fields)
  }

  list(
    fields = extracted_fields,
    matched = matched,
    unmatched = unmatched,
    n_columns = length(extracted_fields)
  )
}

.rap_auto_instance_type <- function(n_columns) {
  if (n_columns > 500) {
    "mem1_ssd1_v2_x36"
  } else if (n_columns > 100) {
    "mem1_ssd1_v2_x16"
  } else if (n_columns > 20) {
    "mem1_ssd1_v2_x8"
  } else {
    "mem1_ssd1_v2_x4"
  }
}

.rap_upload_file <- function(local_path) {
  result <- .rap_dx_run(c("upload", local_path, "--brief"), timeout = 60)
  if (!result$success) {
    stop("Failed to upload field list: ", result$stderr, call. = FALSE)
  }
  trimws(result$stdout)
}

.rap_run_table_exporter <- function(dataset,
                                    file_id,
                                    output,
                                    instance_type,
                                    priority = "low") {
  result <- .rap_dx_run(
    c(
      "run", "table-exporter",
      paste0("-idataset_or_cohort_or_dashboard=", dataset),
      "-ientity=participant",
      paste0("-ifield_names_file_txt=", file_id),
      paste0("-ioutput=", output),
      "-ioutput_format=CSV",
      paste0("--instance-type=", instance_type),
      "--priority", priority,
      "--brief", "--yes"
    ),
    timeout = 60
  )

  if (!result$success) {
    stop("Failed to submit table-exporter job: ", result$stderr, call. = FALSE)
  }

  job_id <- trimws(result$stdout)
  if (!grepl("^job-", job_id)) {
    stop("Unexpected response from table-exporter: ", job_id, call. = FALSE)
  }
  job_id
}

.rap_inform_plan <- function(plan, max_preview = 8) {
  message(sprintf(
    "Matched %d request(s), %d extraction column(s).",
    nrow(plan$matched), plan$n_columns
  ))

  if (nrow(plan$matched) > 0) {
    preview_n <- min(nrow(plan$matched), max_preview)
    for (i in seq_len(preview_n)) {
      row <- plan$matched[i, , drop = FALSE]
      message(sprintf("  %s: %s [%d col(s)]", row$request, row$title, row$n_cols))
    }
    if (nrow(plan$matched) > preview_n) {
      message(sprintf("  ... plus %d more request(s).", nrow(plan$matched) - preview_n))
    }
  }

  if (length(plan$unmatched) > 0) {
    warning(
      sprintf("Skipped unmatched request(s): %s", paste(plan$unmatched, collapse = ", ")),
      call. = FALSE
    )
  }
}

#' Find the RAP Dataset File in the Current Project
#'
#' @param refresh Logical. If TRUE, ignore the cached dataset name and call
#'   \code{dx ls} again.
#' @param timeout Timeout in seconds for the \code{dx ls} call.
#'
#' @return A character scalar naming the detected \code{.dataset} file.
#' @export
#'
rap_find_dataset <- function(refresh = FALSE, timeout = 30) {
  .rap_check_logical(refresh, "refresh")

  if (!isTRUE(refresh) && !is.null(.rap_extract_cache$dataset)) {
    return(.rap_extract_cache$dataset)
  }

  result <- .rap_dx_run("ls", timeout = timeout)
  if (!result$success) {
    stop("Failed to list RAP project files: ", result$stderr, call. = FALSE)
  }

  lines <- strsplit(result$stdout, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines[nzchar(lines)])
  datasets <- lines[grepl("\\.dataset$", lines)]

  if (length(datasets) == 0) {
    stop("No .dataset file found in the current RAP project root.", call. = FALSE)
  }

  dataset <- datasets[[length(datasets)]]
  .rap_extract_cache$dataset <- dataset
  dataset
}

#' List Approved RAP Dataset Fields
#'
#' @param dataset Dataset file name. If NULL, \code{rap_find_dataset()} is used.
#' @param pattern Optional regular expression applied to field names and titles.
#' @param entity Dataset entity. Defaults to \code{"participant"}.
#' @param refresh Logical. If TRUE, bypass the session cache.
#' @param timeout Timeout in seconds for \code{dx extract_dataset --list-fields}.
#'
#' @return A data.frame with columns \code{field_name} and \code{title}.
#' @export
#'
rap_list_fields <- function(dataset = NULL,
                            pattern = NULL,
                            entity = "participant",
                            refresh = FALSE,
                            timeout = 120) {
  .rap_check_logical(refresh, "refresh")

  if (is.null(dataset)) {
    dataset <- rap_find_dataset()
  }
  if (!is.character(dataset) || length(dataset) != 1 || is.na(dataset)) {
    stop("'dataset' must be a single dataset file name.", call. = FALSE)
  }
  if (!is.character(entity) || length(entity) != 1 || is.na(entity)) {
    stop("'entity' must be a single character string.", call. = FALSE)
  }

  cache_key <- paste(dataset, entity, sep = "::")
  if (!isTRUE(refresh) && !is.null(.rap_extract_cache$fields[[cache_key]])) {
    fields <- .rap_extract_cache$fields[[cache_key]]
  } else {
    result <- .rap_dx_run(
      c("extract_dataset", dataset, "--list-fields", "--entities", entity),
      timeout = timeout
    )
    if (!result$success) {
      stop("Failed to list RAP dataset fields: ", result$stderr, call. = FALSE)
    }
    fields <- .rap_parse_fields(result$stdout)
    if (is.null(.rap_extract_cache$fields)) {
      .rap_extract_cache$fields <- new.env(parent = emptyenv())
    }
    .rap_extract_cache$fields[[cache_key]] <- fields
  }

  if (!is.null(pattern)) {
    keep <- grepl(pattern, fields$field_name, perl = TRUE) |
      grepl(pattern, fields$title, perl = TRUE, ignore.case = TRUE)
    fields <- fields[keep, , drop = FALSE]
  }

  rownames(fields) <- NULL
  fields
}

#' Plan a RAP Phenotype Extraction
#'
#' @param field_id UKB numeric field IDs to extract. All instances and arrays are
#'   included.
#' @param field_names Exact RAP dataset column names, such as
#'   \code{"participant.p31"} or \code{"p31"}.
#' @param variables Optional predefined variable names from
#'   \code{get_variable_info()}.
#' @param dataset Dataset file name. If NULL, \code{rap_find_dataset()} is used.
#' @param fields_df Optional cached field listing from \code{rap_list_fields()}.
#' @param entity Dataset entity. Defaults to \code{"participant"}.
#' @param include_eid Logical. Include participant ID automatically.
#' @param table_exporter Logical. If TRUE, return field names in the format
#'   expected by the RAP table-exporter app.
#' @param manifest Optional manifest CSV path in the current RAP session.
#'
#' @return A list containing extraction field names, matched requests, unmatched
#'   requests, dataset, entity, and column counts.
#' @export
#'
rap_plan_extract <- function(field_id = NULL,
                             field_names = NULL,
                             variables = NULL,
                             dataset = NULL,
                             fields_df = NULL,
                             entity = "participant",
                             include_eid = TRUE,
                             table_exporter = FALSE,
                             manifest = NULL) {
  if (is.null(dataset)) {
    dataset <- rap_find_dataset()
  }
  if (is.null(fields_df)) {
    fields_df <- rap_list_fields(dataset = dataset, entity = entity)
  }

  plan <- .rap_plan_fields(
    field_id = field_id,
    field_names = field_names,
    variables = variables,
    fields_df = fields_df,
    entity = entity,
    include_eid = include_eid,
    table_exporter = table_exporter
  )

  if (length(plan$fields) == 0) {
    stop("No extractable fields were matched.", call. = FALSE)
  }

  plan$dataset <- dataset
  plan$entity <- entity
  plan$field_id <- .rap_normalize_field_ids(field_id)
  plan$variables <- if (is.null(variables)) character(0) else variables
  plan$table_exporter <- table_exporter
  class(plan) <- c("rap_extract_plan", class(plan))

  if (!is.null(manifest)) {
    utils::write.csv(plan$matched, manifest, row.names = FALSE)
  }

  plan
}

#' Extract RAP Phenotype Data Synchronously
#'
#' @description
#' Uses \code{dx extract_dataset --fields-file} and reads the RAP-generated
#' result back into R within the active RAP session. This is intended for small
#' to medium extractions. For large phenotype pulls, use
#' \code{rap_submit_extract()}.
#'
#' @param field_id UKB numeric field IDs to extract.
#' @param field_names Exact RAP dataset column names to extract.
#' @param variables Optional predefined variable names from
#'   \code{get_variable_info()}.
#' @param dataset Dataset file name. If NULL, \code{rap_find_dataset()} is used.
#' @param output Optional CSV output path in the current RAP session. If NULL,
#'   a temporary file is used.
#' @param read Logical. If TRUE, read the CSV into R and return a data.table.
#'   If FALSE, return the output path.
#' @param strip_entity_prefix Logical. If TRUE, remove \code{"participant."}
#'   from returned column names.
#' @param dry_run Logical. If TRUE, return the extraction plan without running
#'   \code{dx extract_dataset}.
#' @param timeout Timeout in seconds for the extraction.
#' @param ... Additional arguments passed to \code{rap_plan_extract()}.
#'
#' @return A data.table when \code{read = TRUE}; otherwise the output CSV path.
#'   In dry-run mode, returns a \code{rap_extract_plan}.
#' @export
#'
rap_extract_pheno <- function(field_id = NULL,
                              field_names = NULL,
                              variables = NULL,
                              dataset = NULL,
                              output = NULL,
                              read = TRUE,
                              strip_entity_prefix = FALSE,
                              dry_run = FALSE,
                              timeout = 300,
                              ...) {
  .rap_check_logical(read, "read")
  .rap_check_logical(strip_entity_prefix, "strip_entity_prefix")
  .rap_check_logical(dry_run, "dry_run")

  plan <- rap_plan_extract(
    field_id = field_id,
    field_names = field_names,
    variables = variables,
    dataset = dataset,
    table_exporter = FALSE,
    ...
  )
  .rap_inform_plan(plan)

  if (isTRUE(dry_run)) {
    return(plan)
  }

  .ukb_assert_rap_env("rap_extract_pheno()")

  output_is_temp <- is.null(output)
  if (output_is_temp) {
    output <- tempfile(fileext = ".csv")
    if (isTRUE(read)) {
      on.exit(unlink(output), add = TRUE)
    }
  }

  fields_file <- tempfile(fileext = ".txt")
  on.exit(unlink(fields_file), add = TRUE)
  .rap_write_lf_file(plan$fields, fields_file)

  result <- .rap_dx_run(
    c("extract_dataset", plan$dataset, "--fields-file", fields_file, "-o", output),
    timeout = timeout
  )
  if (!result$success) {
    stop("RAP phenotype extraction failed: ", result$stderr, call. = FALSE)
  }

  if (!isTRUE(read)) {
    return(normalizePath(output, mustWork = FALSE))
  }

  dt <- data.table::fread(output, data.table = TRUE, integer64 = "double")
  if (isTRUE(strip_entity_prefix)) {
    data.table::setnames(dt, sub("^participant\\.", "", names(dt)))
  }
  attr(dt, "rap_extract_plan") <- plan
  dt
}

#' Submit a RAP Table-Exporter Phenotype Extraction Job
#'
#' @description
#' Submits an asynchronous RAP \code{table-exporter} job. This is the preferred
#' interface for large extraction jobs because the work runs on RAP rather than
#' inside the current R session.
#'
#' @param field_id UKB numeric field IDs to extract.
#' @param field_names Exact RAP dataset column names to extract.
#' @param variables Optional predefined variable names from
#'   \code{get_variable_info()}.
#' @param dataset Dataset file name. If NULL, \code{rap_find_dataset()} is used.
#' @param file Output file stem on RAP. Defaults to
#'   \code{"ukba_pheno_YYYYMMDD_HHMMSS"}.
#' @param instance_type DNAnexus instance type. If NULL, selected from the number
#'   of columns.
#' @param priority Job priority: \code{"low"} or \code{"high"}.
#' @param dry_run Logical. If TRUE, return the planned fields and command
#'   metadata without uploading or submitting.
#' @param manifest Optional manifest CSV path in the current RAP session.
#' @param ... Additional arguments passed to \code{rap_plan_extract()}.
#'
#' @return A list with class \code{rap_extract_job} containing job metadata.
#'   In dry-run mode, returns a \code{rap_extract_plan}.
#' @export
#'
rap_submit_extract <- function(field_id = NULL,
                               field_names = NULL,
                               variables = NULL,
                               dataset = NULL,
                               file = NULL,
                               instance_type = NULL,
                               priority = c("low", "high"),
                               dry_run = FALSE,
                               manifest = NULL,
                               ...) {
  priority <- match.arg(priority)
  .rap_check_logical(dry_run, "dry_run")

  output <- if (!is.null(file)) {
    sub("\\.csv$", "", file)
  } else {
    paste0("ukba_pheno_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  plan <- rap_plan_extract(
    field_id = field_id,
    field_names = field_names,
    variables = variables,
    dataset = dataset,
    table_exporter = TRUE,
    manifest = manifest,
    ...
  )
  .rap_inform_plan(plan)

  instance_type <- if (!is.null(instance_type)) instance_type else .rap_auto_instance_type(plan$n_columns)
  command <- c(
    "dx", "run", "table-exporter",
    paste0("-idataset_or_cohort_or_dashboard=", plan$dataset),
    "-ientity=participant",
    "-ifield_names_file_txt=<uploaded-field-list>",
    paste0("-ioutput=", output),
    "-ioutput_format=CSV",
    paste0("--instance-type=", instance_type),
    "--priority", priority,
    "--brief", "--yes"
  )

  plan$output <- output
  plan$instance_type <- instance_type
  plan$priority <- priority
  plan$command <- command

  if (isTRUE(dry_run)) {
    return(plan)
  }

  .ukb_assert_rap_env("rap_submit_extract()")

  fields_file <- tempfile(fileext = ".txt")
  on.exit(unlink(fields_file), add = TRUE)
  .rap_write_lf_file(plan$fields, fields_file)

  message("Uploading field list to RAP...")
  file_id <- .rap_upload_file(fields_file)
  message("Submitting table-exporter job...")
  job_id <- .rap_run_table_exporter(
    dataset = plan$dataset,
    file_id = file_id,
    output = output,
    instance_type = instance_type,
    priority = priority
  )

  job <- list(
    job_id = job_id,
    dataset = plan$dataset,
    output = paste0(output, ".csv"),
    fields_file_id = file_id,
    instance_type = instance_type,
    priority = priority,
    n_columns = plan$n_columns,
    matched = plan$matched,
    unmatched = plan$unmatched,
    fields = plan$fields
  )
  class(job) <- c("rap_extract_job", class(job))
  message(sprintf("Job submitted: %s", job_id))
  job
}
