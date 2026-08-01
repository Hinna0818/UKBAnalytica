# UK Biobank primary-care (GP) records are queried separately from the
# participant-wide phenotype table. The functions in this file intentionally
# do not modify or call build_survival_dataset().

.gp_provider_labels <- c(
  `1` = "England Vision",
  `2` = "Scotland",
  `3` = "England TPP",
  `4` = "Wales"
)

.gp_special_date_labels <- c(
  `1900-01-01` = "no_event_date",
  `1901-01-01` = "before_birth",
  `1902-02-02` = "birth_date_proxy",
  `1903-03-03` = "birth_year_proxy",
  `1909-09-09` = "future_placeholder",
  `2037-07-07` = "future_placeholder"
)

.gp_normalize_code_system <- function(x) {
  raw <- trimws(as.character(x))
  key <- gsub("[^a-z0-9]", "", tolower(raw))
  out <- rep(NA_character_, length(key))
  out[key %in% c("read2", "readv2")] <- "Read2"
  out[key %in% c("read3", "readv3", "ctv3")] <- "CTV3"

  if (anyNA(out)) {
    bad <- unique(raw[is.na(out)])
    stop(
      "Unsupported GP code system(s): ",
      paste(bad, collapse = ", "),
      ". Use 'Read2' or 'CTV3'.",
      call. = FALSE
    )
  }
  out
}

.gp_scalar_date <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be NULL or one valid date.", name), call. = FALSE)
  }
  out <- suppressWarnings(.safe_as_date(x, col_name = name, warn = FALSE))
  if (is.na(out)) {
    stop(sprintf("'%s' must be coercible to Date.", name), call. = FALSE)
  }
  out
}

.gp_sql_string <- function(x) {
  paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
}

.gp_sql_identifier <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[A-Za-z_][A-Za-z0-9_]*$", x)) {
    stop(
      sprintf("'%s' must be a simple Spark SQL identifier.", name),
      call. = FALSE
    )
  }
  paste0("`", x, "`")
}

#' Parse a GP Disease Query
#'
#' @description
#' Validates a user-supplied Read v2 / CTV3 disease code table before a
#' `gp_clinical` query is sent to UK Biobank RAP Spark SQL. Matching is exact:
#' Read code punctuation and leading zeroes are preserved.
#'
#' @param codes A data.frame with columns `disease`, `code_system`, and `code`.
#'   `code_system` accepts common forms of Read v2 (`"Read2"`, `"read_2"`) and
#'   CTV3 (`"CTV3"`, `"read_3"`).
#' @param date_from Optional earliest GP event date.
#' @param date_to Optional latest GP event date.
#' @param providers Optional integer provider codes (1 to 4).
#' @param eids Optional participant identifiers. This is intended for small
#'   development cohorts; large cohorts should be represented by a Spark-side
#'   cohort join in a future extension.
#'
#' @return An object of class `ukb_gp_query`.
#' @export
#'
#' @examples
#' query <- parse_gp_query(data.frame(
#'   disease = c("Hypertension", "Type 2 diabetes"),
#'   code_system = c("Read2", "CTV3"),
#'   code = c("G20..", "C10E.")
#' ))
parse_gp_query <- function(codes,
                           date_from = NULL,
                           date_to = NULL,
                           providers = NULL,
                           eids = NULL) {
  required <- c("disease", "code_system", "code")
  if (!is.data.frame(codes) || !all(required %in% names(codes))) {
    stop(
      "'codes' must be a data.frame with columns: disease, code_system, code.",
      call. = FALSE
    )
  }

  code_dt <- data.table::as.data.table(codes)[, required, with = FALSE]
  code_dt <- data.table::copy(code_dt)
  code_dt[, disease := trimws(as.character(disease))]
  code_dt[, code_system := .gp_normalize_code_system(code_system)]
  code_dt[, code := trimws(as.character(code))]

  invalid <- is.na(code_dt$disease) | !nzchar(code_dt$disease) |
    is.na(code_dt$code) | !nzchar(code_dt$code)
  if (any(invalid)) {
    stop("'disease' and 'code' cannot contain missing or empty values.", call. = FALSE)
  }
  code_dt <- unique(code_dt)

  date_from <- .gp_scalar_date(date_from, "date_from")
  date_to <- .gp_scalar_date(date_to, "date_to")
  if (!is.null(date_from) && !is.null(date_to) && date_from > date_to) {
    stop("'date_from' cannot be later than 'date_to'.", call. = FALSE)
  }

  if (!is.null(providers)) {
    providers <- suppressWarnings(as.integer(providers))
    if (length(providers) == 0L || anyNA(providers) ||
        any(!providers %in% 1:4)) {
      stop("'providers' must contain UKB GP provider codes 1, 2, 3, or 4.", call. = FALSE)
    }
    providers <- unique(providers)
  }

  if (!is.null(eids)) {
    eids <- trimws(as.character(eids))
    if (length(eids) == 0L || anyNA(eids) || any(!grepl("^[0-9]+$", eids))) {
      stop("'eids' must contain non-missing numeric participant identifiers.", call. = FALSE)
    }
    eids <- unique(eids)
  }

  out <- list(
    codes = code_dt,
    date_from = date_from,
    date_to = date_to,
    providers = providers,
    eids = eids
  )
  class(out) <- "ukb_gp_query"
  out
}

#' @export
print.ukb_gp_query <- function(x, ...) {
  cat("<ukb_gp_query>\n")
  cat("  Diseases:", data.table::uniqueN(x$codes$disease), "\n")
  cat("  Exact codes:", nrow(x$codes), "\n")
  cat(
    "  Date range:",
    if (is.null(x$date_from)) "<open>" else format(x$date_from),
    "to",
    if (is.null(x$date_to)) "<open>" else format(x$date_to),
    "\n"
  )
  invisible(x)
}

#' Plan a RAP Spark SQL GP Query
#'
#' @description
#' Creates auditable Spark SQL that selects only requested Read v2 / CTV3
#' records from the stable RAP `gp_clinical` view. It does not execute the
#' query or load GP records into R.
#'
#' @param query A `ukb_gp_query` from [parse_gp_query()].
#' @param database Optional RAP Spark database name. If `NULL`, the active
#'   database of the supplied Spark connection is used.
#' @param table GP clinical table name. Defaults to the stable
#'   `"gp_clinical"` view.
#' @param cohort_table Optional Spark table/view containing an `eid` column.
#'   When supplied, a left-semi join restricts GP records before collection.
#' @param cohort_database Optional database containing `cohort_table`. Defaults
#'   to `database`.
#' @param collect Return record-level SQL (`"records"`) or push participant by
#'   disease aggregation into Spark (`"summary"`).
#'
#' @return A `ukb_gp_query_plan` containing the SQL and normalized query.
#' @export
rap_plan_gp_query <- function(query,
                              database = NULL,
                              table = "gp_clinical",
                              cohort_table = NULL,
                              cohort_database = database,
                              collect = c("records", "summary")) {
  if (!inherits(query, "ukb_gp_query")) {
    stop("'query' must be created by parse_gp_query().", call. = FALSE)
  }
  collect <- match.arg(collect)

  table_ref <- .gp_sql_identifier(table, "table")
  if (!is.null(database)) {
    table_ref <- paste0(
      .gp_sql_identifier(database, "database"),
      ".",
      table_ref
    )
  }
  cohort_ref <- NULL
  if (!is.null(cohort_table)) {
    cohort_ref <- .gp_sql_identifier(cohort_table, "cohort_table")
    if (!is.null(cohort_database)) {
      cohort_ref <- paste0(
        .gp_sql_identifier(cohort_database, "cohort_database"),
        ".",
        cohort_ref
      )
    }
  }

  read2 <- unique(query$codes[code_system == "Read2", code])
  ctv3 <- unique(query$codes[code_system == "CTV3", code])
  code_filters <- character(0)
  if (length(read2) > 0L) {
    code_filters <- c(
      code_filters,
      paste0("`g`.`read_2` IN (", paste(.gp_sql_string(read2), collapse = ", "), ")")
    )
  }
  if (length(ctv3) > 0L) {
    code_filters <- c(
      code_filters,
      paste0("`g`.`read_3` IN (", paste(.gp_sql_string(ctv3), collapse = ", "), ")")
    )
  }

  filters <- paste0("(", paste(code_filters, collapse = " OR "), ")")
  if (!is.null(query$date_from)) {
    filters <- c(
      filters,
      paste0("CAST(`g`.`event_dt` AS DATE) >= DATE '", format(query$date_from), "'")
    )
  }
  if (!is.null(query$date_to)) {
    filters <- c(
      filters,
      paste0("CAST(`g`.`event_dt` AS DATE) <= DATE '", format(query$date_to), "'")
    )
  }
  if (!is.null(query$providers)) {
    filters <- c(
      filters,
      paste0("CAST(`g`.`data_provider` AS INT) IN (", paste(query$providers, collapse = ", "), ")")
    )
  }
  if (!is.null(query$eids)) {
    filters <- c(
      filters,
      paste0("CAST(`g`.`eid` AS STRING) IN (", paste(.gp_sql_string(query$eids), collapse = ", "), ")")
    )
  }

  from_sql <- paste0("FROM ", table_ref, " AS `g`")
  if (!is.null(cohort_ref)) {
    from_sql <- paste0(
      from_sql,
      "\nLEFT SEMI JOIN ", cohort_ref, " AS `c`",
      "\n  ON CAST(`g`.`eid` AS STRING) = CAST(`c`.`eid` AS STRING)"
    )
  }
  record_sql <- paste(
    "SELECT",
    paste0(
      "  `g`.`eid` AS `eid`, `g`.`data_provider` AS `data_provider`,\n",
      "  `g`.`event_dt` AS `event_dt`, `g`.`read_2` AS `read_2`,\n",
      "  `g`.`read_3` AS `read_3`, `g`.`value1` AS `value1`,\n",
      "  `g`.`value2` AS `value2`, `g`.`value3` AS `value3`"
    ),
    from_sql,
    paste0("WHERE ", paste(filters, collapse = "\n  AND ")),
    sep = "\n"
  )

  sql <- record_sql
  if (identical(collect, "summary")) {
    code_rows <- apply(
      as.data.frame(query$codes),
      1L,
      function(row) {
        paste0(
          "SELECT ", .gp_sql_string(row[["disease"]]), " AS `disease`, ",
          .gp_sql_string(row[["code_system"]]), " AS `code_system`, ",
          .gp_sql_string(row[["code"]]), " AS `code`"
        )
      }
    )
    special_dates <- paste(
      .gp_sql_string(names(.gp_special_date_labels)),
      collapse = ", "
    )
    valid_date_sql <- paste0(
      "CAST(`r`.`event_dt` AS DATE) IS NOT NULL AND ",
      "CAST(`r`.`event_dt` AS DATE) NOT IN (", special_dates, ")"
    )
    sql <- paste0(
      "WITH `ukba_gp_records` AS (\n",
      record_sql,
      "\n),\n`ukba_gp_code_rows` AS (\n",
      "  SELECT `eid`, `data_provider`, `event_dt`, 'Read2' AS `code_system`, `read_2` AS `code`\n",
      "  FROM `ukba_gp_records`\n",
      "  WHERE `read_2` IS NOT NULL AND TRIM(`read_2`) <> ''\n",
      "  UNION ALL\n",
      "  SELECT `eid`, `data_provider`, `event_dt`, 'CTV3' AS `code_system`, `read_3` AS `code`\n",
      "  FROM `ukba_gp_records`\n",
      "  WHERE `read_3` IS NOT NULL AND TRIM(`read_3`) <> ''\n",
      "),\n`ukba_gp_codes` AS (\n  ",
      paste(code_rows, collapse = "\n  UNION ALL\n  "),
      "\n)\n",
      "SELECT\n",
      "  `r`.`eid` AS `eid`, `m`.`disease` AS `disease`,\n",
      "  CAST(1 AS INT) AS `gp_case`,\n",
      "  MIN(CASE WHEN ", valid_date_sql,
      " THEN CAST(`r`.`event_dt` AS DATE) ELSE NULL END) AS `first_gp_date`,\n",
      "  CAST(COUNT(*) AS INT) AS `n_gp_records`,\n",
      "  CAST(COUNT(DISTINCT CONCAT_WS(':', `r`.`code_system`, `r`.`code`)) AS INT) AS `n_gp_codes`,\n",
      "  CAST(SUM(CASE WHEN ", valid_date_sql,
      " THEN 1 ELSE 0 END) AS INT) AS `n_usable_dates`,\n",
      "  'GP' AS `source`\n",
      "FROM `ukba_gp_code_rows` AS `r`\n",
      "INNER JOIN `ukba_gp_codes` AS `m`\n",
      "  ON `r`.`code_system` = `m`.`code_system` AND `r`.`code` = `m`.`code`\n",
      "GROUP BY `r`.`eid`, `m`.`disease`"
    )
  }

  out <- list(
    query = query,
    database = database,
    table = table,
    cohort_table = cohort_table,
    cohort_database = cohort_database,
    collect = collect,
    sql = sql
  )
  class(out) <- "ukb_gp_query_plan"
  out
}

#' @export
print.ukb_gp_query_plan <- function(x, ...) {
  cat("<ukb_gp_query_plan>\n")
  cat("  Collect:", x$collect, "\n")
  cat(x$sql, "\n")
  invisible(x)
}

#' Run a RAP Spark SQL GP Query
#'
#' @description
#' Executes a planned GP query through an existing Spark DBI connection. The
#' function adds a row guard before collecting results into the R worker.
#' It is intended for a Spark-enabled UK Biobank RAP R notebook.
#'
#' @param plan A plan from [rap_plan_gp_query()].
#' @param connection An active Spark DBI connection, commonly created with
#'   `sparklyr::spark_connect()`.
#' @param output Optional `.csv` or `.tsv` path for the filtered records.
#' @param max_rows Maximum rows allowed to be collected into R. The executed
#'   SQL requests at most `max_rows + 1` rows and stops if the limit is exceeded.
#' @param dry_run If `TRUE`, return the SQL without executing it.
#'
#' @return A data.table containing filtered GP records or a participant-disease
#'   summary, according to `plan$collect`. If `output` is supplied, the same
#'   result is also written to disk.
#' @export
rap_run_gp_query <- function(plan,
                             connection = NULL,
                             output = NULL,
                             max_rows = 1000000L,
                             dry_run = FALSE) {
  if (!inherits(plan, "ukb_gp_query_plan")) {
    stop("'plan' must be created by rap_plan_gp_query().", call. = FALSE)
  }
  .rap_check_logical(dry_run, "dry_run")
  if (isTRUE(dry_run)) {
    return(plan)
  }

  .ukb_assert_rap_env("rap_run_gp_query()")
  if (is.null(connection)) {
    stop(
      "'connection' is required. Supply an active Spark DBI connection.",
      call. = FALSE
    )
  }
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required to execute GP Spark SQL.", call. = FALSE)
  }
  if (length(max_rows) != 1L || is.na(max_rows) ||
      max_rows < 1 || max_rows != as.integer(max_rows)) {
    stop("'max_rows' must be one positive integer.", call. = FALSE)
  }
  max_rows <- as.integer(max_rows)

  guarded_sql <- paste0(
    "SELECT * FROM (\n",
    sub(";[[:space:]]*$", "", plan$sql),
    "\n) AS `ukba_gp_filtered`\nLIMIT ",
    max_rows + 1L
  )
  result <- DBI::dbGetQuery(connection, guarded_sql)
  if (nrow(result) > max_rows) {
    stop(
      sprintf(
        paste0(
          "GP query matched more than %s rows. Refine the codes/date range or ",
          "increase 'max_rows' explicitly."
        ),
        format(max_rows, big.mark = ",", scientific = FALSE)
      ),
      call. = FALSE
    )
  }

  result <- data.table::as.data.table(result)
  attr(result, "gp_query") <- plan$query
  attr(result, "gp_query_plan") <- plan

  if (!is.null(output)) {
    if (!is.character(output) || length(output) != 1L || is.na(output)) {
      stop("'output' must be one CSV or TSV path.", call. = FALSE)
    }
    ext <- tolower(tools::file_ext(output))
    if (!ext %in% c("csv", "tsv")) {
      stop("'output' must end in .csv or .tsv.", call. = FALSE)
    }
    data.table::fwrite(result, output, sep = if (ext == "tsv") "\t" else ",")
    attr(result, "gp_output") <- normalizePath(output, mustWork = FALSE)
  }

  result
}

.gp_empty_parsed <- function(with_disease = FALSE) {
  out <- data.table::data.table(
    gp_record_id = integer(0),
    eid = integer(0),
    event_date = as.Date(character(0)),
    code_system = character(0),
    code = character(0),
    data_provider = integer(0),
    provider_label = character(0),
    value1 = character(0),
    value2 = character(0),
    value3 = character(0),
    date_quality = character(0),
    usable_for_timing = logical(0),
    source = character(0)
  )
  if (isTRUE(with_disease)) {
    out[, disease := character(0)]
    data.table::setcolorder(out, c("gp_record_id", "eid", "disease", setdiff(names(out), c("gp_record_id", "eid", "disease"))))
  }
  out
}

#' Parse Filtered UKB GP Clinical Records
#'
#' @description
#' Standardizes filtered `gp_clinical` records into one row per Read v2 / CTV3
#' code. This function does not read the full RAP table and does not call Spark.
#' If a GP query is supplied, exact code matching adds the disease label.
#'
#' @param data Filtered `gp_clinical` data.frame/data.table.
#' @param query Optional `ukb_gp_query`. If omitted, the query attached by
#'   [rap_run_gp_query()] is used when available.
#' @param special_date_action Keep and flag UKB placeholder dates (`"flag"`) or
#'   remove those records (`"exclude"`). Missing dates remain because they can
#'   still represent coded events, but are never usable for timing.
#' @param keep_empty_code Retain records where Read v2 / CTV3 code is empty.
#'
#' @return A standardized long-format data.table.
#' @export
parse_gp_clinical <- function(data,
                              query = NULL,
                              special_date_action = c("flag", "exclude"),
                              keep_empty_code = FALSE) {
  special_date_action <- match.arg(special_date_action)
  .rap_check_logical(keep_empty_code, "keep_empty_code")
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }

  if (is.null(query)) {
    query <- attr(data, "gp_query", exact = TRUE)
  }
  if (!is.null(query) && !inherits(query, "ukb_gp_query")) {
    stop("'query' must be NULL or created by parse_gp_query().", call. = FALSE)
  }

  required <- c("eid", "data_provider", "event_dt", "read_2", "read_3")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required gp_clinical column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dt <- data.table::copy(data.table::as.data.table(data))
  for (col in c("value1", "value2", "value3")) {
    if (!col %in% names(dt)) {
      dt[, (col) := NA_character_]
    }
  }
  if (nrow(dt) == 0L) {
    return(.gp_empty_parsed(with_disease = !is.null(query)))
  }

  dt[, gp_record_id := .I]
  dt[, event_date := .safe_as_date(event_dt, col_name = "event_dt")]
  dt[, read_2 := as.character(read_2)]
  dt[, read_3 := as.character(read_3)]

  long <- data.table::melt(
    dt,
    id.vars = c(
      "gp_record_id", "eid", "data_provider", "event_date",
      "value1", "value2", "value3"
    ),
    measure.vars = c("read_2", "read_3"),
    variable.name = "code_field",
    value.name = "code",
    variable.factor = FALSE
  )
  long[, code_system := ifelse(code_field == "read_2", "Read2", "CTV3")]
  long[, code_field := NULL]
  long[, code := trimws(as.character(code))]

  empty_code <- is.na(long$code) | !nzchar(long$code)
  if (!isTRUE(keep_empty_code)) {
    long <- long[!empty_code]
  }

  date_chr <- as.character(long$event_date)
  long[, date_quality := "valid"]
  long[is.na(event_date), date_quality := "missing"]
  for (special_date in names(.gp_special_date_labels)) {
    long[date_chr == special_date, date_quality := .gp_special_date_labels[[special_date]]]
  }
  long[, usable_for_timing := date_quality == "valid"]

  provider_key <- as.character(long$data_provider)
  long[, provider_label := unname(.gp_provider_labels[provider_key])]
  long[is.na(provider_label), provider_label := "Unknown"]
  long[, source := "GP"]

  if (identical(special_date_action, "exclude")) {
    special_labels <- unname(.gp_special_date_labels)
    long <- long[!date_quality %in% special_labels]
  }

  if (!is.null(query)) {
    long <- merge(
      long,
      query$codes,
      by = c("code_system", "code"),
      all = FALSE,
      allow.cartesian = TRUE,
      sort = FALSE
    )
  }

  preferred <- c(
    "gp_record_id", "eid",
    if (!is.null(query)) "disease",
    "event_date", "code_system", "code", "data_provider", "provider_label",
    "value1", "value2", "value3", "date_quality", "usable_for_timing",
    "source"
  )
  data.table::setcolorder(long, preferred)
  data.table::setorder(long, eid, event_date, na.last = TRUE)
  long
}

#' Summarise GP Diagnoses to Participant Level
#'
#' @description
#' Collapses matched, standardized GP records to one row per participant and
#' disease. Only participants with at least one matching GP code are returned;
#' absence from this result must not be interpreted as a confirmed control
#' without separately evaluating GP registration and coverage.
#'
#' @param gp_records Output from [parse_gp_clinical()] with a `disease` column.
#'
#' @return A data.table with GP case indicator, earliest usable GP event date,
#'   record counts, and source.
#' @export
summarise_gp_diagnoses <- function(gp_records) {
  if (!is.data.frame(gp_records)) {
    stop("'gp_records' must be a data.frame or data.table.", call. = FALSE)
  }
  required <- c(
    "eid", "disease", "gp_record_id", "event_date",
    "usable_for_timing", "code_system", "code"
  )
  missing_cols <- setdiff(required, names(gp_records))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing standardized GP column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dt <- data.table::copy(data.table::as.data.table(gp_records))
  if (nrow(dt) == 0L) {
    return(data.table::data.table(
      eid = integer(0),
      disease = character(0),
      gp_case = integer(0),
      first_gp_date = as.Date(character(0)),
      n_gp_records = integer(0),
      n_gp_codes = integer(0),
      n_usable_dates = integer(0),
      source = character(0)
    ))
  }

  out <- dt[
    ,
    {
      valid_dates <- event_date[usable_for_timing & !is.na(event_date)]
      list(
        gp_case = 1L,
        first_gp_date = if (length(valid_dates) > 0L) min(valid_dates) else as.Date(NA),
        n_gp_records = data.table::uniqueN(gp_record_id),
        n_gp_codes = data.table::uniqueN(paste(code_system, code, sep = ":")),
        n_usable_dates = length(valid_dates),
        source = "GP"
      )
    },
    by = .(eid, disease)
  ]
  data.table::setorder(out, eid, disease)
  out
}

.gp_date_quality <- function(x) {
  x_chr <- as.character(x)
  out <- rep("valid", length(x))
  out[is.na(x)] <- "missing"
  for (special_date in names(.gp_special_date_labels)) {
    out[!is.na(x_chr) & x_chr == special_date] <-
      .gp_special_date_labels[[special_date]]
  }
  out
}

.gp_empty_registrations <- function() {
  out <- data.table::data.table(
    gp_registration_id = integer(0),
    eid = integer(0),
    data_provider = integer(0),
    provider_label = character(0),
    reg_date = as.Date(character(0)),
    deduct_date = as.Date(character(0)),
    reg_date_quality = character(0),
    deduct_date_quality = character(0),
    open_ended = logical(0),
    reversed_interval = logical(0),
    interval_quality = character(0),
    usable_for_coverage = logical(0),
    source = character(0)
  )
  class(out) <- c("ukb_gp_registrations", class(out))
  out
}

#' Parse UKB GP Registration Records
#'
#' @description
#' Standardizes the four official columns in the UK Biobank
#' `gp_registrations` table and flags the pseudo dates defined by UKB
#' Data-Coding 819. Registration records remain separate from clinical disease
#' records.
#'
#' @param data A data.frame/data.table with `eid`, `data_provider`, `reg_date`,
#'   and `deduct_date`.
#' @param special_date_action Keep and flag UKB pseudo dates (`"flag"`) or
#'   remove registration rows containing a pseudo date (`"exclude"`).
#'
#' @return A long-format data.table with one row per GP registration record.
#' @export
parse_gp_registrations <- function(data,
                                   special_date_action = c("flag", "exclude")) {
  special_date_action <- match.arg(special_date_action)
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }

  required <- c("eid", "data_provider", "reg_date", "deduct_date")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required gp_registrations column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dt <- data.table::copy(data.table::as.data.table(data))
  if (nrow(dt) == 0L) {
    return(.gp_empty_registrations())
  }

  dt[, gp_registration_id := .I]
  dt[, reg_date := .safe_as_date(reg_date, col_name = "reg_date")]
  dt[, deduct_date := .safe_as_date(deduct_date, col_name = "deduct_date")]
  dt[, reg_date_quality := .gp_date_quality(reg_date)]
  dt[, deduct_date_quality := .gp_date_quality(deduct_date)]

  provider_key <- as.character(dt$data_provider)
  dt[, provider_label := unname(.gp_provider_labels[provider_key])]
  dt[is.na(provider_label), provider_label := "Unknown"]

  dt[, open_ended := deduct_date_quality == "missing"]
  dt[, reversed_interval :=
    reg_date_quality == "valid" &
      deduct_date_quality == "valid" &
      deduct_date < reg_date
  ]
  dt[, interval_quality := "valid"]
  dt[open_ended == TRUE, interval_quality := "open_ended"]
  dt[reg_date_quality != "valid", interval_quality := "invalid_start"]
  dt[
    reg_date_quality == "valid" &
      !deduct_date_quality %in% c("valid", "missing"),
    interval_quality := "invalid_end"
  ]
  dt[reversed_interval == TRUE, interval_quality := "reversed"]
  dt[, usable_for_coverage :=
    reg_date_quality == "valid" &
      deduct_date_quality %in% c("valid", "missing") &
      !reversed_interval
  ]
  dt[, source := "GPRegistration"]

  if (identical(special_date_action, "exclude")) {
    special_labels <- unique(unname(.gp_special_date_labels))
    dt <- dt[
      !reg_date_quality %in% special_labels &
        !deduct_date_quality %in% special_labels
    ]
  }

  preferred <- c(
    "gp_registration_id", "eid", "data_provider", "provider_label",
    "reg_date", "deduct_date", "reg_date_quality", "deduct_date_quality",
    "open_ended", "reversed_interval", "interval_quality",
    "usable_for_coverage", "source"
  )
  data.table::setcolorder(dt, preferred)
  data.table::setorder(dt, eid, reg_date, na.last = TRUE)
  class(dt) <- unique(c("ukb_gp_registrations", class(dt)))
  dt
}

.gp_coverage_days <- function(start, end) {
  keep <- !is.na(start) & !is.na(end) & end >= start
  start <- start[keep]
  end <- end[keep]
  if (length(start) == 0L) {
    return(NA_integer_)
  }

  ord <- order(start, end)
  start <- start[ord]
  end <- end[ord]
  current_start <- start[[1]]
  current_end <- end[[1]]
  total <- 0L

  if (length(start) > 1L) {
    for (i in 2:length(start)) {
      if (start[[i]] <= current_end + 1) {
        if (end[[i]] > current_end) {
          current_end <- end[[i]]
        }
      } else {
        total <- total + as.integer(current_end - current_start) + 1L
        current_start <- start[[i]]
        current_end <- end[[i]]
      }
    }
  }
  as.integer(total + as.integer(current_end - current_start) + 1L)
}

.gp_extract_eids <- function(x, name, allow_empty = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    if (!"eid" %in% names(x)) {
      stop(sprintf("'%s' data.frame must contain an 'eid' column.", name), call. = FALSE)
    }
    x <- x[["eid"]]
  }
  if (length(x) == 0L) {
    if (isTRUE(allow_empty)) {
      return(x)
    }
    stop(sprintf("'%s' must contain participant identifiers.", name), call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("'%s' must contain non-missing participant identifiers.", name), call. = FALSE)
  }
  unique(x)
}

.gp_empty_coverage <- function(eid_template = integer(0)) {
  data.table::data.table(
    eid = eid_template,
    coverage_start = as.Date(rep(NA_character_, length(eid_template))),
    coverage_end = as.Date(rep(NA_character_, length(eid_template))),
    coverage_days = rep(NA_integer_, length(eid_template)),
    n_registration_records = integer(length(eid_template)),
    n_usable_registration_records = integer(length(eid_template)),
    n_bounded_intervals = integer(length(eid_template)),
    has_open_registration = rep(FALSE, length(eid_template)),
    gp_providers = rep(NA_character_, length(eid_template)),
    has_gp_clinical_record = rep(FALSE, length(eid_template)),
    coverage_status = rep("unknown", length(eid_template))
  )
}

#' Summarise UKB GP Registration Coverage
#'
#' @description
#' Summarises parsed GP registration intervals to one row per participant.
#' Open-ended intervals are bounded only when `observation_end` is supplied.
#' Filtered clinical records may be supplied as secondary evidence, but they
#' produce `"clinical_only"` rather than confirmed registration coverage.
#'
#' @param registrations Output from [parse_gp_registrations()].
#' @param clinical_records Optional filtered `gp_clinical` records (raw or
#'   parsed). Only the `eid` column is used.
#' @param participants Optional participant vector or data.frame with `eid`.
#'   Requested participants without registration/clinical evidence are retained
#'   with `coverage_status = "unknown"`.
#' @param observation_end Optional study/data-release cut-off used to close
#'   missing deduction dates and cap later deduction dates.
#'
#' @return A participant-level data.table containing coverage bounds, days,
#'   registration counts, providers, and coverage status.
#' @export
summarise_gp_coverage <- function(registrations,
                                  clinical_records = NULL,
                                  participants = NULL,
                                  observation_end = NULL) {
  if (!is.data.frame(registrations)) {
    stop("'registrations' must be a data.frame or data.table.", call. = FALSE)
  }
  required <- c(
    "eid", "data_provider", "provider_label", "reg_date", "deduct_date",
    "reg_date_quality", "deduct_date_quality", "open_ended",
    "usable_for_coverage"
  )
  missing_cols <- setdiff(required, names(registrations))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing parsed GP registration column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  observation_end <- .gp_scalar_date(observation_end, "observation_end")
  participant_eids <- .gp_extract_eids(participants, "participants")
  clinical_eids <- NULL
  if (!is.null(clinical_records)) {
    clinical_eids <- .gp_extract_eids(
      clinical_records,
      "clinical_records",
      allow_empty = TRUE
    )
  }

  dt <- data.table::copy(data.table::as.data.table(registrations))
  if (nrow(dt) == 0L) {
    all_eids <- unique(c(participant_eids, clinical_eids))
    out <- .gp_empty_coverage(all_eids)
    if (length(clinical_eids) > 0L) {
      out[eid %in% clinical_eids, `:=`(
        has_gp_clinical_record = TRUE,
        coverage_status = "clinical_only"
      )]
    }
    return(out)
  }

  dt[, coverage_start_i := as.Date(NA)]
  dt[usable_for_coverage == TRUE, coverage_start_i := reg_date]
  dt[, coverage_end_i := as.Date(NA)]
  dt[
    usable_for_coverage & deduct_date_quality == "valid",
    coverage_end_i := deduct_date
  ]
  if (!is.null(observation_end)) {
    dt[
      usable_for_coverage & open_ended,
      coverage_end_i := observation_end
    ]
    dt[
      !is.na(coverage_end_i) & coverage_end_i > observation_end,
      coverage_end_i := observation_end
    ]
    dt[
      !is.na(coverage_start_i) & coverage_start_i > observation_end,
      `:=`(coverage_start_i = as.Date(NA), coverage_end_i = as.Date(NA))
    ]
  }
  dt[
    !is.na(coverage_start_i) &
      !is.na(coverage_end_i) &
      coverage_end_i < coverage_start_i,
    coverage_end_i := as.Date(NA)
  ]
  dt[, usable_in_window := usable_for_coverage & !is.na(coverage_start_i)]

  out <- dt[
    ,
    {
      usable_start <- coverage_start_i[!is.na(coverage_start_i)]
      bounded_end <- coverage_end_i[
        !is.na(coverage_start_i) & !is.na(coverage_end_i)
      ]
      bounded_start <- coverage_start_i[
        !is.na(coverage_start_i) & !is.na(coverage_end_i)
      ]
      unresolved_open <- any(
        usable_in_window & open_ended & is.na(coverage_end_i)
      )
      provider_values <- sort(unique(provider_label[usable_in_window]))
      provider_values <- provider_values[!is.na(provider_values)]

      list(
        coverage_start = if (length(usable_start) > 0L) min(usable_start) else as.Date(NA),
        coverage_end = if (unresolved_open) {
          as.Date(NA)
        } else if (length(bounded_end) > 0L) {
          max(bounded_end)
        } else {
          as.Date(NA)
        },
        coverage_days = if (unresolved_open) {
          NA_integer_
        } else {
          .gp_coverage_days(bounded_start, bounded_end)
        },
        n_registration_records = .N,
        n_usable_registration_records = sum(usable_in_window),
        n_bounded_intervals = length(bounded_end),
        has_open_registration = any(usable_in_window & open_ended),
        gp_providers = if (length(provider_values) > 0L) {
          paste(provider_values, collapse = "; ")
        } else {
          NA_character_
        }
      )
    },
    by = eid
  ]

  out[, has_gp_clinical_record := eid %in% clinical_eids]
  out[, coverage_status := ifelse(
    n_usable_registration_records > 0L,
    "confirmed",
    ifelse(has_gp_clinical_record, "clinical_only", "unknown")
  )]

  all_eids <- unique(c(participant_eids, out$eid, clinical_eids))
  missing_eids <- setdiff(all_eids, out$eid)
  if (length(missing_eids) > 0L) {
    extra <- .gp_empty_coverage(missing_eids)
    extra[eid %in% clinical_eids, `:=`(
      has_gp_clinical_record = TRUE,
      coverage_status = "clinical_only"
    )]
    out <- data.table::rbindlist(list(out, extra), use.names = TRUE, fill = TRUE)
  }

  data.table::setorder(out, eid)
  out
}

.gp_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      x < 0 || x != as.integer(x)) {
    stop(sprintf("'%s' must be one non-negative integer.", name), call. = FALSE)
  }
  as.integer(x)
}

.gp_fraction <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("'%s' must be one number between 0 and 1.", name), call. = FALSE)
  }
  as.numeric(x)
}

.gp_resolve_participant_date <- function(participants, x, name, required = TRUE) {
  if (is.null(x)) {
    if (isTRUE(required)) {
      stop(sprintf("'%s' is required.", name), call. = FALSE)
    }
    return(as.Date(rep(NA_character_, nrow(participants))))
  }

  if (is.character(x) && length(x) == 1L && x %in% names(participants)) {
    values <- participants[[x]]
  } else if (length(x) == 1L) {
    values <- rep(x, nrow(participants))
  } else if (length(x) == nrow(participants)) {
    values <- x
  } else {
    stop(
      sprintf(
        "'%s' must be one date, a participant date-column name, or one value per participant.",
        name
      ),
      call. = FALSE
    )
  }

  parsed <- .safe_as_date(values, col_name = name)
  if (isTRUE(required) && anyNA(parsed)) {
    stop(sprintf("'%s' cannot contain missing or invalid dates.", name), call. = FALSE)
  }
  parsed
}

.gp_merge_intervals <- function(start, end) {
  keep <- !is.na(start) & !is.na(end) & end >= start
  start_i <- as.integer(start[keep])
  end_i <- as.integer(end[keep])
  if (length(start_i) == 0L) {
    return(data.table::data.table(
      interval_start = as.Date(character(0)),
      interval_end = as.Date(character(0))
    ))
  }

  ord <- order(start_i, end_i)
  start_i <- start_i[ord]
  end_i <- end_i[ord]
  merged_start <- integer(0)
  merged_end <- integer(0)
  current_start <- start_i[[1L]]
  current_end <- end_i[[1L]]

  if (length(start_i) > 1L) {
    for (i in 2:length(start_i)) {
      if (start_i[[i]] <= current_end + 1L) {
        current_end <- max(current_end, end_i[[i]])
      } else {
        merged_start <- c(merged_start, current_start)
        merged_end <- c(merged_end, current_end)
        current_start <- start_i[[i]]
        current_end <- end_i[[i]]
      }
    }
  }
  merged_start <- c(merged_start, current_start)
  merged_end <- c(merged_end, current_end)

  data.table::data.table(
    interval_start = as.Date(merged_start, origin = "1970-01-01"),
    interval_end = as.Date(merged_end, origin = "1970-01-01")
  )
}

#' Assess GP Observability in a Participant-Specific Study Window
#'
#' @description
#' Uses parsed `gp_registrations` intervals to determine whether each
#' participant has sufficient primary-care observation for classification as a
#' non-case. Overlapping and adjacent registration intervals are merged before
#' calculating coverage, so duplicated registrations do not inflate observed
#' days. Open registrations are usable only when `observation_end` is supplied.
#'
#' @param registrations Output from [parse_gp_registrations()].
#' @param participants Optional participant vector or data.frame containing
#'   `eid`. A data.frame may also contain participant-specific study dates.
#' @param window_start Study-window start. Supply one date, a column name in
#'   `participants`, or one value per participant.
#' @param window_end Study-window end, using the same forms as `window_start`.
#' @param index_date Optional index/baseline date. When supplied, the
#'   participant must be covered on this date and continuous lookback/follow-up
#'   around the date can be required.
#' @param observation_end Optional GP release/provider cut-off used to close
#'   open registrations and cap later deduction dates.
#' @param min_lookback_days Minimum continuous covered days immediately before
#'   `index_date`.
#' @param min_followup_days Minimum continuous covered days immediately after
#'   `index_date`.
#' @param min_coverage_fraction Minimum fraction of the complete study window
#'   covered by valid registration intervals.
#'
#' @return A participant-level data.table with observed days, coverage fraction,
#'   continuous lookback/follow-up, `control_eligible`, and exclusion reason.
#' @export
assess_gp_observability <- function(registrations,
                                    participants = NULL,
                                    window_start,
                                    window_end,
                                    index_date = NULL,
                                    observation_end = NULL,
                                    min_lookback_days = 0L,
                                    min_followup_days = 0L,
                                    min_coverage_fraction = 1) {
  if (!is.data.frame(registrations)) {
    stop("'registrations' must be a data.frame or data.table.", call. = FALSE)
  }
  required <- c(
    "eid", "reg_date", "deduct_date", "deduct_date_quality",
    "open_ended", "usable_for_coverage"
  )
  missing_cols <- setdiff(required, names(registrations))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing parsed GP registration column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  min_lookback_days <- .gp_nonnegative_integer(
    min_lookback_days,
    "min_lookback_days"
  )
  min_followup_days <- .gp_nonnegative_integer(
    min_followup_days,
    "min_followup_days"
  )
  min_coverage_fraction <- .gp_fraction(
    min_coverage_fraction,
    "min_coverage_fraction"
  )
  if (is.null(index_date) &&
      (min_lookback_days > 0L || min_followup_days > 0L)) {
    stop(
      "'index_date' is required when lookback or follow-up is requested.",
      call. = FALSE
    )
  }
  observation_end <- .gp_scalar_date(observation_end, "observation_end")

  if (is.null(participants)) {
    participant_dt <- data.table::data.table(
      eid = unique(registrations[["eid"]])
    )
  } else if (is.data.frame(participants)) {
    if (!"eid" %in% names(participants)) {
      stop("'participants' data.frame must contain an 'eid' column.", call. = FALSE)
    }
    participant_dt <- data.table::copy(data.table::as.data.table(participants))
  } else {
    participant_dt <- data.table::data.table(eid = participants)
  }
  if (nrow(participant_dt) == 0L || anyNA(participant_dt$eid)) {
    stop("'participants' must contain non-missing participant identifiers.", call. = FALSE)
  }
  if (anyDuplicated(participant_dt$eid)) {
    stop("'participants' must contain one row per participant.", call. = FALSE)
  }

  windows <- data.table::data.table(eid = participant_dt$eid)
  windows[, window_start := .gp_resolve_participant_date(
    participant_dt,
    window_start,
    "window_start"
  )]
  windows[, window_end := .gp_resolve_participant_date(
    participant_dt,
    window_end,
    "window_end"
  )]
  windows[, index_date := .gp_resolve_participant_date(
    participant_dt,
    index_date,
    "index_date",
    required = !is.null(index_date)
  )]
  windows[, valid_window := window_end >= window_start]
  windows[, index_in_window := is.na(index_date) |
    (index_date >= window_start & index_date <= window_end)]

  registration_dt <- data.table::copy(data.table::as.data.table(registrations))
  registration_dt <- registration_dt[eid %in% windows$eid]
  registration_dt[, interval_start := as.Date(NA)]
  registration_dt[usable_for_coverage == TRUE, interval_start := reg_date]
  registration_dt[, interval_end := as.Date(NA)]
  registration_dt[
    usable_for_coverage == TRUE & deduct_date_quality == "valid",
    interval_end := deduct_date
  ]
  if (!is.null(observation_end)) {
    registration_dt[
      usable_for_coverage == TRUE & open_ended == TRUE,
      interval_end := observation_end
    ]
    registration_dt[
      !is.na(interval_end) & interval_end > observation_end,
      interval_end := observation_end
    ]
  }
  registration_dt <- registration_dt[
    !is.na(interval_start) & !is.na(interval_end) & interval_end >= interval_start
  ]

  if (nrow(registration_dt) > 0L) {
    intervals <- registration_dt[
      ,
      .gp_merge_intervals(interval_start, interval_end),
      by = eid
    ]
  } else {
    intervals <- data.table::data.table(
      eid = windows$eid[0],
      interval_start = as.Date(character(0)),
      interval_end = as.Date(character(0))
    )
  }

  joined <- merge(windows, intervals, by = "eid", all.x = TRUE, sort = FALSE)
  joined[, overlap_start_i := pmax(
    as.integer(window_start),
    as.integer(interval_start),
    na.rm = FALSE
  )]
  joined[, overlap_end_i := pmin(
    as.integer(window_end),
    as.integer(interval_end),
    na.rm = FALSE
  )]
  joined[, overlap_days := data.table::fifelse(
    !is.na(overlap_start_i) & !is.na(overlap_end_i) &
      overlap_end_i >= overlap_start_i,
    overlap_end_i - overlap_start_i + 1L,
    0L
  )]
  joined[, contains_index := !is.na(index_date) &
    !is.na(interval_start) &
    index_date >= interval_start & index_date <= interval_end]
  joined[, continuous_lookback_i := data.table::fifelse(
    contains_index,
    as.integer(index_date - pmax(interval_start, window_start)),
    NA_integer_
  )]
  joined[, continuous_followup_i := data.table::fifelse(
    contains_index,
    as.integer(pmin(interval_end, window_end) - index_date),
    NA_integer_
  )]

  out <- joined[
    ,
    .(
      window_start = window_start[[1L]],
      window_end = window_end[[1L]],
      index_date = index_date[[1L]],
      valid_window = valid_window[[1L]],
      index_in_window = index_in_window[[1L]],
      target_window_days = as.integer(window_end[[1L]] - window_start[[1L]]) + 1L,
      covered_window_days = as.integer(sum(overlap_days, na.rm = TRUE)),
      coverage_at_index = if (all(is.na(index_date))) NA else any(contains_index),
      continuous_lookback_days = if (all(is.na(index_date))) {
        NA_integer_
      } else if (any(contains_index)) {
        max(continuous_lookback_i, na.rm = TRUE)
      } else {
        0L
      },
      continuous_followup_days = if (all(is.na(index_date))) {
        NA_integer_
      } else if (any(contains_index)) {
        max(continuous_followup_i, na.rm = TRUE)
      } else {
        0L
      }
    ),
    by = eid
  ]
  out[, coverage_fraction := data.table::fifelse(
    valid_window & target_window_days > 0L,
    covered_window_days / target_window_days,
    NA_real_
  )]
  index_required <- !all(is.na(out$index_date))
  out[, control_eligible :=
    valid_window & index_in_window &
      covered_window_days > 0L &
      coverage_fraction >= min_coverage_fraction
  ]
  if (index_required) {
    out[, control_eligible := control_eligible &
      coverage_at_index == TRUE &
      continuous_lookback_days >= min_lookback_days &
      continuous_followup_days >= min_followup_days
    ]
  }
  out[, gp_observability_reason := data.table::fcase(
    !valid_window, "invalid_window",
    !index_in_window, "index_outside_window",
    covered_window_days == 0L, "no_coverage_in_window",
    !is.na(coverage_at_index) & coverage_at_index == FALSE,
      "not_covered_at_index",
    !is.na(continuous_lookback_days) &
      continuous_lookback_days < min_lookback_days,
      "insufficient_lookback",
    !is.na(continuous_followup_days) &
      continuous_followup_days < min_followup_days,
      "insufficient_followup",
    coverage_fraction < min_coverage_fraction,
      "insufficient_window_coverage",
    control_eligible == TRUE, "eligible",
    default = "ineligible"
  )]
  data.table::setorder(out, eid)
  out
}

#' Integrate GP Diagnoses and Registration Coverage
#'
#' @description
#' Creates one row per participant and disease. A matched GP code is assigned
#' `gp_case = 1`; a non-match is assigned `0` only when registration coverage
#' is confirmed; otherwise it remains `NA`.
#'
#' @param gp_diagnoses Output from [summarise_gp_diagnoses()].
#' @param gp_coverage Output from [summarise_gp_coverage()].
#' @param participants Optional participant vector or data.frame with `eid`.
#' @param diseases Optional disease names. Required when no diagnosis rows are
#'   present.
#' @param gp_observability Optional output from [assess_gp_observability()].
#'   When supplied, unmatched participants receive `gp_case = 0` only when
#'   `control_eligible` is `TRUE`.
#'
#' @return A participant-by-disease GP phenotype data.table.
#' @export
integrate_gp_results <- function(gp_diagnoses,
                                 gp_coverage,
                                 participants = NULL,
                                 diseases = NULL,
                                 gp_observability = NULL) {
  if (!is.data.frame(gp_diagnoses) || !is.data.frame(gp_coverage)) {
    stop("'gp_diagnoses' and 'gp_coverage' must be data.frames.", call. = FALSE)
  }
  if (!all(c("eid", "disease", "gp_case", "first_gp_date") %in% names(gp_diagnoses))) {
    stop("'gp_diagnoses' is missing required diagnosis summary columns.", call. = FALSE)
  }
  if (!all(c("eid", "coverage_status") %in% names(gp_coverage))) {
    stop("'gp_coverage' is missing required coverage summary columns.", call. = FALSE)
  }
  if (!is.null(gp_observability) &&
      (!is.data.frame(gp_observability) ||
       !all(c("eid", "control_eligible") %in% names(gp_observability)))) {
    stop(
      "'gp_observability' must contain 'eid' and 'control_eligible'.",
      call. = FALSE
    )
  }

  diagnoses <- data.table::copy(data.table::as.data.table(gp_diagnoses))
  coverage <- data.table::copy(data.table::as.data.table(gp_coverage))
  participant_eids <- .gp_extract_eids(participants, "participants")
  if (is.null(participant_eids)) {
    participant_eids <- unique(c(diagnoses$eid, coverage$eid))
  }
  if (is.null(diseases)) {
    diseases <- unique(as.character(diagnoses$disease))
  }
  diseases <- unique(trimws(as.character(diseases)))
  if (length(participant_eids) == 0L || length(diseases) == 0L ||
      anyNA(diseases) || any(!nzchar(diseases))) {
    stop("At least one participant and disease are required.", call. = FALSE)
  }

  grid <- data.table::CJ(
    eid = participant_eids,
    disease = diseases,
    unique = TRUE,
    sorted = FALSE
  )
  if ("source" %in% names(diagnoses)) {
    data.table::setnames(diagnoses, "source", "diagnosis_source")
  }
  out <- merge(
    grid,
    diagnoses,
    by = c("eid", "disease"),
    all.x = TRUE,
    sort = FALSE
  )
  out <- merge(out, coverage, by = "eid", all.x = TRUE, sort = FALSE)
  out[is.na(coverage_status), coverage_status := "unknown"]
  if (!is.null(gp_observability)) {
    observability <- data.table::copy(data.table::as.data.table(gp_observability))
    if (anyDuplicated(observability$eid)) {
      stop("'gp_observability' must contain one row per participant.", call. = FALSE)
    }
    out <- merge(out, observability, by = "eid", all.x = TRUE, sort = FALSE)
  }

  matched <- !is.na(out$gp_case) & out$gp_case == 1L
  eligible_control <- if (is.null(gp_observability)) {
    out$coverage_status == "confirmed"
  } else {
    !is.na(out$control_eligible) & out$control_eligible
  }
  out[, gp_case := as.integer(NA)]
  out[matched, gp_case := 1L]
  out[!matched & eligible_control, gp_case := 0L]
  out[, gp_case_reason := ifelse(
    !is.na(gp_case) & gp_case == 1L,
    "matched_code",
    ifelse(
      !is.na(gp_case) & gp_case == 0L,
      "registered_no_match",
      ifelse(
        !is.null(gp_observability),
        "ineligible_observation_window",
        "insufficient_coverage"
      )
    )
  )]

  count_cols <- intersect(
    c("n_gp_records", "n_gp_codes", "n_usable_dates"),
    names(out)
  )
  for (col in count_cols) {
    data.table::set(out, which(is.na(out[[col]])), col, 0L)
  }
  data.table::setorder(out, eid, disease)
  out
}

.gp_registration_sql <- function(query,
                                 database = NULL,
                                 table = "gp_registrations",
                                 cohort_table = NULL,
                                 cohort_database = database) {
  table_ref <- .gp_sql_identifier(table, "registration_table")
  if (!is.null(database)) {
    table_ref <- paste0(
      .gp_sql_identifier(database, "database"),
      ".",
      table_ref
    )
  }
  cohort_ref <- NULL
  if (!is.null(cohort_table)) {
    cohort_ref <- .gp_sql_identifier(cohort_table, "cohort_table")
    if (!is.null(cohort_database)) {
      cohort_ref <- paste0(
        .gp_sql_identifier(cohort_database, "cohort_database"),
        ".",
        cohort_ref
      )
    }
  }

  filters <- character(0)
  if (!is.null(query$providers)) {
    filters <- c(
      filters,
      paste0(
        "CAST(`g`.`data_provider` AS INT) IN (",
        paste(query$providers, collapse = ", "),
        ")"
      )
    )
  }
  if (!is.null(query$eids)) {
    filters <- c(
      filters,
      paste0(
        "CAST(`g`.`eid` AS STRING) IN (",
        paste(.gp_sql_string(query$eids), collapse = ", "),
        ")"
      )
    )
  }

  sql <- paste0(
    "SELECT\n",
    "  `g`.`eid` AS `eid`, `g`.`data_provider` AS `data_provider`,\n",
    "  `g`.`reg_date` AS `reg_date`, `g`.`deduct_date` AS `deduct_date`\n",
    "FROM ", table_ref, " AS `g`"
  )
  if (!is.null(cohort_ref)) {
    sql <- paste0(
      sql,
      "\nLEFT SEMI JOIN ", cohort_ref, " AS `c`",
      "\n  ON CAST(`g`.`eid` AS STRING) = CAST(`c`.`eid` AS STRING)"
    )
  }
  if (length(filters) > 0L) {
    sql <- paste0(sql, "\nWHERE ", paste(filters, collapse = "\n  AND "))
  }
  sql
}

.gp_collect_sql <- function(connection,
                            sql,
                            max_rows,
                            label,
                            output = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required to execute GP Spark SQL.", call. = FALSE)
  }
  if (length(max_rows) != 1L || is.na(max_rows) ||
      max_rows < 1 || max_rows != as.integer(max_rows)) {
    stop("'max_rows' must be one positive integer.", call. = FALSE)
  }
  max_rows <- as.integer(max_rows)
  guarded_sql <- paste0(
    "SELECT * FROM (\n",
    sub(";[[:space:]]*$", "", sql),
    "\n) AS `ukba_gp_filtered`\nLIMIT ",
    max_rows + 1L
  )
  result <- DBI::dbGetQuery(connection, guarded_sql)
  if (nrow(result) > max_rows) {
    stop(
      sprintf(
        "%s query matched more than %s rows. Refine the query or increase the row limit explicitly.",
        label,
        format(max_rows, big.mark = ",", scientific = FALSE)
      ),
      call. = FALSE
    )
  }

  result <- data.table::as.data.table(result)
  if (!is.null(output)) {
    if (!is.character(output) || length(output) != 1L || is.na(output)) {
      stop("'output' must be one CSV or TSV path.", call. = FALSE)
    }
    ext <- tolower(tools::file_ext(output))
    if (!ext %in% c("csv", "tsv")) {
      stop("'output' must end in .csv or .tsv.", call. = FALSE)
    }
    data.table::fwrite(result, output, sep = if (ext == "tsv") "\t" else ",")
  }
  result
}

.gp_is_diagnosis_summary <- function(x) {
  is.data.frame(x) && all(c(
    "eid", "disease", "gp_case", "first_gp_date",
    "n_gp_records", "n_gp_codes", "n_usable_dates"
  ) %in% names(x))
}

.gp_standardize_diagnosis_summary <- function(x) {
  if (!.gp_is_diagnosis_summary(x)) {
    stop("GP diagnosis summary columns are incomplete.", call. = FALSE)
  }
  out <- data.table::copy(data.table::as.data.table(x))
  out[, disease := as.character(disease)]
  out[, gp_case := as.integer(gp_case)]
  out[, first_gp_date := .safe_as_date(
    first_gp_date,
    col_name = "first_gp_date"
  )]
  for (col in c("n_gp_records", "n_gp_codes", "n_usable_dates")) {
    data.table::set(out, j = col, value = as.integer(out[[col]]))
  }
  if (!"source" %in% names(out)) {
    out[, source := "GP"]
  }
  data.table::setorder(out, eid, disease)
  out
}

#' Run an End-to-End UKB GP Phenotype Workflow
#'
#' @description
#' Runs the independent GP workflow from a parsed disease query through RAP
#' Spark SQL extraction, clinical/registration parsing, coverage summary, and
#' participant-level phenotype integration. The query is the first argument so
#' the function can follow [parse_gp_query()] with the base R pipe (`|>`).
#'
#' @param query A `ukb_gp_query` from [parse_gp_query()].
#' @param connection An active RAP Spark DBI connection. Not required when both
#'   `clinical_data` and `registration_data` are supplied for local testing.
#' @param database Optional RAP Spark database name.
#' @param participants Optional participant vector or data.frame with `eid`.
#' @param observation_end Optional study/data-release cut-off used to close
#'   open registration intervals.
#' @param clinical_data Optional pre-extracted `gp_clinical` records.
#' @param registration_data Optional pre-extracted `gp_registrations` records.
#' @param clinical_output Optional local CSV/TSV path for filtered clinical
#'   records.
#' @param registration_output Optional local CSV/TSV path for registration
#'   records.
#' @param max_clinical_rows Maximum clinical records collected into R.
#' @param max_registration_rows Maximum registration records collected into R.
#' @param dry_run Return both SQL plans without executing them.
#' @param cohort_table Optional Spark cohort table/view containing `eid`. Both
#'   clinical and registration queries use a left-semi join to this table.
#' @param cohort_database Optional database containing `cohort_table`. Defaults
#'   to `database`.
#' @param collect RAP clinical collection mode. `"summary"` performs disease
#'   aggregation in Spark; `"records"` collects matched record-level rows.
#' @param window_start Optional study-window start passed to
#'   [assess_gp_observability()]. Supply together with `window_end` to enable
#'   strict control eligibility.
#' @param window_end Optional study-window end.
#' @param index_date Optional participant index/baseline date.
#' @param min_lookback_days Minimum continuous coverage before `index_date`.
#' @param min_followup_days Minimum continuous coverage after `index_date`.
#' @param min_coverage_fraction Minimum covered fraction of the study window.
#'
#' @return A participant-by-disease GP phenotype data.table, or a dry-run plan.
#' @export
run_gp_workflow <- function(query,
                            connection = NULL,
                            database = NULL,
                            participants = NULL,
                            observation_end = NULL,
                            clinical_data = NULL,
                            registration_data = NULL,
                            clinical_output = NULL,
                            registration_output = NULL,
                            max_clinical_rows = 1000000L,
                            max_registration_rows = 1000000L,
                            dry_run = FALSE,
                            cohort_table = NULL,
                            cohort_database = database,
                            collect = c("summary", "records"),
                            window_start = NULL,
                            window_end = NULL,
                            index_date = NULL,
                            min_lookback_days = 0L,
                            min_followup_days = 0L,
                            min_coverage_fraction = 1) {
  if (!inherits(query, "ukb_gp_query")) {
    stop("'query' must be created by parse_gp_query().", call. = FALSE)
  }
  .rap_check_logical(dry_run, "dry_run")
  collect <- match.arg(collect)
  has_window <- !is.null(window_start) || !is.null(window_end)
  if (xor(is.null(window_start), is.null(window_end))) {
    stop("Supply both 'window_start' and 'window_end', or neither.", call. = FALSE)
  }

  clinical_plan <- rap_plan_gp_query(
    query,
    database = database,
    cohort_table = cohort_table,
    cohort_database = cohort_database,
    collect = collect
  )
  registration_sql <- .gp_registration_sql(
    query,
    database = database,
    cohort_table = cohort_table,
    cohort_database = cohort_database
  )
  if (isTRUE(dry_run)) {
    out <- list(
      query = query,
      collect = collect,
      cohort_table = cohort_table,
      clinical_sql = clinical_plan$sql,
      registration_sql = registration_sql
    )
    class(out) <- "ukb_gp_workflow_plan"
    return(out)
  }

  supplied_data <- !is.null(clinical_data) || !is.null(registration_data)
  if (supplied_data && (is.null(clinical_data) || is.null(registration_data))) {
    stop(
      "Supply both 'clinical_data' and 'registration_data', or neither.",
      call. = FALSE
    )
  }
  if (!supplied_data) {
    .ukb_assert_rap_env("run_gp_workflow()")
    if (is.null(connection)) {
      stop("'connection' is required for RAP Spark execution.", call. = FALSE)
    }
    clinical_data <- rap_run_gp_query(
      clinical_plan,
      connection = connection,
      output = clinical_output,
      max_rows = max_clinical_rows
    )
    registration_data <- .gp_collect_sql(
      connection = connection,
      sql = registration_sql,
      max_rows = max_registration_rows,
      label = "GP registration",
      output = registration_output
    )
  }

  diagnoses <- if (.gp_is_diagnosis_summary(clinical_data)) {
    .gp_standardize_diagnosis_summary(clinical_data)
  } else {
    summarise_gp_diagnoses(
      parse_gp_clinical(clinical_data, query = query)
    )
  }
  parsed_registrations <- parse_gp_registrations(registration_data)
  coverage <- summarise_gp_coverage(
    parsed_registrations,
    clinical_records = clinical_data,
    participants = participants,
    observation_end = observation_end
  )
  observability <- NULL
  if (has_window) {
    observability <- assess_gp_observability(
      parsed_registrations,
      participants = participants,
      window_start = window_start,
      window_end = window_end,
      index_date = index_date,
      observation_end = observation_end,
      min_lookback_days = min_lookback_days,
      min_followup_days = min_followup_days,
      min_coverage_fraction = min_coverage_fraction
    )
  }
  result <- integrate_gp_results(
    diagnoses,
    gp_coverage = coverage,
    participants = participants,
    diseases = unique(query$codes$disease),
    gp_observability = observability
  )
  attr(result, "gp_query") <- query
  attr(result, "clinical_sql") <- clinical_plan$sql
  attr(result, "registration_sql") <- registration_sql
  result
}

#' @export
print.ukb_gp_workflow_plan <- function(x, ...) {
  cat("<ukb_gp_workflow_plan>\n")
  cat("  Collect:", x$collect, "\n")
  if (!is.null(x$cohort_table)) {
    cat("  Cohort table:", x$cohort_table, "\n")
  }
  cat("\nClinical SQL:\n", x$clinical_sql, "\n", sep = "")
  cat("\nRegistration SQL:\n", x$registration_sql, "\n", sep = "")
  invisible(x)
}
