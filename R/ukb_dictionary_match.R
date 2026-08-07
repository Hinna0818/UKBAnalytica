#' Download the official RAP data dictionary
#'
#' @description
#' Runs \code{dx extract_dataset -ddd} inside UK Biobank RAP and returns the
#' generated official data dictionary CSV path. This function checks that it is
#' being executed in a RAP-like environment before calling \code{dx}.
#'
#' @param dataset RAP \code{.dataset} file or record identifier. If \code{NULL},
#'   \code{\link{rap_find_dataset}} is used.
#' @param output_dir Directory where the dictionary files should be written.
#' @param delimiter Output delimiter passed to \code{dx extract_dataset}.
#' @param timeout Timeout in seconds for the \code{dx} command.
#' @param require_rap Logical. If \code{TRUE}, require a RAP-like environment.
#'
#' @return Path to the generated \code{*.data_dictionary.csv} file.
#' @export
#'
ukb_download_rap_dictionary <- function(dataset = NULL,
                                        output_dir = ".",
                                        delimiter = ",",
                                        timeout = 600,
                                        require_rap = TRUE) {
  if (isTRUE(require_rap)) {
    .ukb_assert_rap_env("ukb_download_rap_dictionary()")
  }
  if (!nzchar(Sys.which("dx"))) {
    stop("The dx command-line tool is required to download the RAP data dictionary.", call. = FALSE)
  }
  if (is.null(dataset)) {
    dataset <- rap_find_dataset()
  }
  if (!is.character(dataset) || length(dataset) != 1L || is.na(dataset) || !nzchar(dataset)) {
    stop("'dataset' must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("'output_dir' must be a single non-empty path.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  before <- list.files(output_dir, pattern = "data_dictionary.*\\.csv$", full.names = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)

  result <- .rap_dx_run(
    c("extract_dataset", dataset, "-ddd", "--delimiter", delimiter),
    timeout = timeout
  )
  if (!result$success) {
    stop("Failed to download RAP data dictionary: ", result$stderr, call. = FALSE)
  }

  after <- list.files(output_dir, pattern = "data_dictionary.*\\.csv$|\\.data_dictionary\\.csv$", full.names = TRUE)
  candidates <- setdiff(after, before)
  if (length(candidates) == 0L) {
    candidates <- after
  }
  if (length(candidates) == 0L) {
    stop("No data_dictionary CSV file was created by dx extract_dataset -ddd.", call. = FALSE)
  }

  candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
  normalizePath(candidates[[1]], mustWork = TRUE)
}

#' Query UK Biobank dictionary metadata
#'
#' @description
#' Searches UK Biobank variable metadata using a RAP-generated official data
#' dictionary and the UKBAnalytica Chinese dictionary. Chinese queries are first
#' matched against the built-in Chinese dictionary and translated into English
#' candidate terms before matching the official dictionary. English queries,
#' field IDs, and RAP/UKB column names are searched directly in the official
#' dictionary.
#'
#' This function is intended for RAP use. By default it requires a RAP-like
#' environment; set \code{require_rap = FALSE} only for package development or
#' tests using simulated dictionaries.
#'
#' @param query Character vector of query terms, Chinese variable names, English
#'   names, UKB field IDs, or UKB/RAP column names.
#' @param official_dict Optional official RAP data dictionary CSV. If
#'   \code{NULL}, \code{\link{ukb_download_rap_dictionary}} is called.
#' @param zh_dict Optional Chinese dictionary CSV. Defaults to the
#'   UKBAnalytica built-in Chinese dictionary.
#' @param dataset RAP \code{.dataset} file used when \code{official_dict} is
#'   \code{NULL}.
#' @param output_dir Directory used when downloading the official dictionary.
#' @param language Query language. \code{"auto"} detects Chinese, field IDs, and
#'   column names.
#' @param translation_map Optional data.frame with columns \code{zh} and
#'   \code{en}, or a named character vector mapping Chinese terms to English
#'   query terms.
#' @param max_results Maximum official dictionary matches returned per query.
#' @param min_score Minimum internal matching score for official dictionary
#'   matches.
#' @param require_rap Logical. Require a RAP-like environment before querying.
#' @param timeout Timeout in seconds when downloading the official dictionary.
#'
#' @return A list of class \code{ukb_dictionary_query} with official matches,
#'   Chinese matches, query metadata, and source paths.
#' @export
#'
ukb_query_dictionary <- function(query,
                                 official_dict = NULL,
                                 zh_dict = NULL,
                                 dataset = NULL,
                                 output_dir = tempdir(),
                                 language = c("auto", "zh", "en", "field_id", "column"),
                                 translation_map = NULL,
                                 max_results = 20,
                                 min_score = 0.35,
                                 require_rap = TRUE,
                                 timeout = 600) {
  language <- match.arg(language)
  if (isTRUE(require_rap)) {
    .ukb_assert_rap_env("ukb_query_dictionary()")
  }
  if (!is.character(query) || length(query) == 0L || anyNA(query) || any(!nzchar(trimws(query)))) {
    stop("'query' must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.numeric(max_results) || length(max_results) != 1L || is.na(max_results) || max_results < 1) {
    stop("'max_results' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(min_score) || length(min_score) != 1L || is.na(min_score)) {
    stop("'min_score' must be a numeric scalar.", call. = FALSE)
  }

  if (is.null(official_dict)) {
    official_dict <- ukb_download_rap_dictionary(
      dataset = dataset,
      output_dir = output_dir,
      timeout = timeout,
      require_rap = require_rap
    )
  }
  official <- .ukb_read_official_dictionary(official_dict)
  zh <- .ukb_read_zh_dictionary(zh_dict)
  trans <- .ukb_dictionary_translation_map(translation_map)

  official_results <- list()
  zh_results <- list()
  query_info <- data.frame(
    query = trimws(query),
    detected_language = character(length(query)),
    english_candidate = character(length(query)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(query)) {
    q <- trimws(query[[i]])
    q_lang <- .ukb_dictionary_detect_query(q, language)
    query_info$detected_language[[i]] <- q_lang

    zh_match <- if (q_lang == "zh") .ukb_search_zh_dictionary(q, zh, max_results = max_results) else .ukb_empty_zh_result()
    zh_results[[i]] <- zh_match

    candidate <- if (q_lang == "zh") {
      .ukb_make_english_candidate(q, zh_match, trans)
    } else {
      q
    }
    query_info$english_candidate[[i]] <- candidate

    official_match <- .ukb_search_official_dictionary(
      query = q,
      english_candidate = candidate,
      official = official,
      query_language = q_lang,
      max_results = max_results,
      min_score = min_score
    )
    official_results[[i]] <- official_match
  }

  official_matches <- do.call(rbind, official_results)
  zh_matches <- do.call(rbind, zh_results)
  rownames(official_matches) <- NULL
  rownames(zh_matches) <- NULL

  out <- list(
    query = trimws(query),
    query_info = query_info,
    official_matches = official_matches,
    chinese_matches = zh_matches,
    official_dict = normalizePath(official_dict, mustWork = FALSE),
    zh_dict = if (is.null(zh_dict)) "ukb_dictionary_zh" else normalizePath(zh_dict, mustWork = FALSE),
    require_rap = require_rap
  )
  class(out) <- c("ukb_dictionary_query", class(out))
  out
}

#' @export
print.ukb_dictionary_query <- function(x, ...) {
  cat("UKB dictionary query\n")
  cat("  query: ", paste(x$query, collapse = "; "), "\n", sep = "")
  cat("  official matches: ", nrow(x$official_matches), "\n", sep = "")
  cat("  Chinese matches: ", nrow(x$chinese_matches), "\n", sep = "")
  if (nrow(x$official_matches) > 0L) {
    print(head(x$official_matches, 10), row.names = FALSE)
  }
  invisible(x)
}

#' Validate requested columns against a data object
#'
#' @description
#' Checks whether requested UKB/RAP columns are present in a data.frame or a
#' character vector of available column names. The function can optionally treat
#' \code{participant.p31} and \code{p31} as equivalent.
#'
#' @param data A data.frame/data.table or a character vector of available column
#'   names.
#' @param columns Character vector of requested column names.
#' @param ignore_entity_prefix Logical. If \code{TRUE}, compare both original
#'   names and names with a leading \code{"participant."} prefix removed.
#' @param error Logical. If \code{TRUE}, stop when any requested column is
#'   missing.
#'
#' @return A data.frame of class \code{ukb_column_validation}.
#' @export
#'
#' @examples
#' dat <- data.frame(eid = 1:3, p31 = c(0, 1, 0))
#' ukb_validate_columns(dat, c("eid", "p31", "p21022"))
ukb_validate_columns <- function(data,
                                 columns,
                                 ignore_entity_prefix = TRUE,
                                 error = FALSE) {
  if (is.data.frame(data)) {
    available <- names(data)
  } else if (is.character(data)) {
    available <- data
  } else {
    stop("'data' must be a data.frame/data.table or a character vector of column names.", call. = FALSE)
  }
  if (!is.character(columns) || length(columns) == 0L || anyNA(columns)) {
    stop("'columns' must be a non-empty character vector.", call. = FALSE)
  }

  available_cmp <- .ukb_normalize_column_names(available, ignore_entity_prefix)
  requested_cmp <- .ukb_normalize_column_names(columns, ignore_entity_prefix)
  present <- requested_cmp %in% available_cmp
  matched <- rep(NA_character_, length(columns))
  matched[present] <- available[match(requested_cmp[present], available_cmp)]

  suggestions <- vapply(seq_along(columns), function(i) {
    if (present[[i]]) {
      return(NA_character_)
    }
    hits <- adist(requested_cmp[[i]], available_cmp)
    if (length(hits) == 0L) {
      return(NA_character_)
    }
    ord <- order(as.numeric(hits))[seq_len(min(3L, length(hits)))]
    paste(available[ord], collapse = "; ")
  }, character(1))

  out <- data.frame(
    requested = columns,
    present = present,
    matched_column = matched,
    suggestion = suggestions,
    stringsAsFactors = FALSE
  )
  class(out) <- c("ukb_column_validation", class(out))
  attr(out, "all_present") <- all(present)

  if (isTRUE(error) && any(!present)) {
    stop("Missing column(s): ", paste(columns[!present], collapse = ", "), call. = FALSE)
  }
  out
}

.ukb_read_zh_dictionary <- function(path = NULL) {
  if (is.null(path)) {
    data_env <- new.env(parent = emptyenv())
    data("ukb_dictionary_zh", package = "UKBAnalytica", envir = data_env)
    if (!exists("ukb_dictionary_zh", envir = data_env, inherits = FALSE)) {
      if (exists("ukb_dictionary_zh", envir = parent.frame(), inherits = TRUE)) {
        x <- get("ukb_dictionary_zh", envir = parent.frame(), inherits = TRUE)
      } else {
        stop("Built-in data object 'ukb_dictionary_zh' was not found.", call. = FALSE)
      }
    } else {
      x <- get("ukb_dictionary_zh", envir = data_env, inherits = FALSE)
    }
  } else {
    if (!file.exists(path)) {
      stop("Chinese dictionary file does not exist: ", path, call. = FALSE)
    }
    x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  }
  expected <- paste0("\u7c7b\u522b", 1:6)
  missing <- setdiff(expected, names(x))
  if (length(missing) > 0L) {
    stop("Chinese dictionary is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x[] <- lapply(x, function(z) {
    z <- as.character(z)
    z[is.na(z)] <- ""
    trimws(z)
  })
  x$zh_path <- apply(x[expected], 1, function(row) paste(row[nzchar(row)], collapse = " > "))
  x$zh_variable <- x[["\u7c7b\u522b6"]]
  x
}

.ukb_read_official_dictionary <- function(path) {
  if (!file.exists(path)) {
    stop("Official dictionary file does not exist: ", path, call. = FALSE)
  }
  raw <- fread(path, data.table = FALSE, encoding = "UTF-8")
  id_col <- .field_metadata_pick_col(raw, c("ID", "FieldID", "Field ID", "field_id", "name"))
  title_col <- .field_metadata_pick_col(raw, c("title", "Title", "Field", "field"))
  name_col <- .field_metadata_pick_col(raw, c("name", "Name", "field_name", "Field"))
  category_col <- .field_metadata_pick_col(raw, c("folder_path", "Category", "Path", "Category Path", "folder"))
  type_col <- .field_metadata_pick_col(raw, c("type", "ValueType", "Value Type", "value_type"))
  coding_col <- .field_metadata_pick_col(raw, c("coding_name", "Coding", "Coding ID", "coding"))
  units_col <- .field_metadata_pick_col(raw, c("units", "Units", "Unit"))
  instance_col <- .field_metadata_pick_col(raw, c("Instance", "instance"))
  array_col <- .field_metadata_pick_col(raw, c("Array", "array"))
  link_col <- .field_metadata_pick_col(raw, c("linkout", "Link", "URL", "link"))

  if (is.null(id_col) || is.null(title_col)) {
    stop("Official dictionary must contain field/column IDs and titles.", call. = FALSE)
  }

  id <- as.character(raw[[id_col]])
  official <- data.frame(
    FieldID = .ukb_parse_field_id_from_column(id),
    UKBColumn = id,
    Field = .ukb_pick_text(raw, name_col, title_col),
    Title = .ukb_pick_text(raw, title_col, name_col),
    Category = .ukb_pick_text(raw, category_col, NULL),
    ValueType = .ukb_pick_text(raw, type_col, NULL),
    Units = .ukb_pick_text(raw, units_col, NULL),
    Coding = .ukb_pick_text(raw, coding_col, NULL),
    Instance = .ukb_pick_text(raw, instance_col, NULL),
    Array = .ukb_pick_text(raw, array_col, NULL),
    Link = .ukb_pick_text(raw, link_col, NULL),
    stringsAsFactors = FALSE
  )
  official$search_text <- paste(
    official$FieldID,
    official$UKBColumn,
    official$Field,
    official$Title,
    official$Category,
    official$ValueType,
    official$Coding,
    sep = " "
  )
  official
}

.ukb_pick_text <- function(raw, primary, secondary = NULL) {
  out <- if (!is.null(primary) && primary %in% names(raw)) as.character(raw[[primary]]) else rep("", nrow(raw))
  if (!is.null(secondary) && secondary %in% names(raw)) {
    idx <- is.na(out) | !nzchar(trimws(out))
    out[idx] <- as.character(raw[[secondary]])[idx]
  }
  out[is.na(out)] <- ""
  trimws(out)
}

.ukb_parse_field_id_from_column <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  out[x == "eid" | x == "participant.eid"] <- "eid"
  is_field <- grepl("(^|\\.)p[0-9]+", x)
  out[is_field] <- sub("^.*\\.?p([0-9]+).*$", "\\1", x[is_field], perl = TRUE)
  numeric_id <- grepl("^[0-9]+$", x)
  out[numeric_id] <- x[numeric_id]
  out
}

.ukb_dictionary_detect_query <- function(query, language) {
  if (language != "auto") {
    return(language)
  }
  if (grepl("[\u4e00-\u9fff]", query)) {
    return("zh")
  }
  if (grepl("^[0-9]+$", query)) {
    return("field_id")
  }
  if (grepl("^(participant\\.)?p[0-9]+|^eid$", query, ignore.case = TRUE)) {
    return("column")
  }
  "en"
}

.ukb_empty_zh_result <- function() {
  data.frame(
    query = character(0),
    zh_variable = character(0),
    zh_path = character(0),
    zh_score = numeric(0),
    stringsAsFactors = FALSE
  )
}

.ukb_search_zh_dictionary <- function(query, zh, max_results) {
  text <- paste(zh$zh_variable, zh$zh_path)
  score <- .ukb_text_match_score(query, text)
  keep <- score > 0
  if (!any(keep)) {
    return(.ukb_empty_zh_result())
  }
  out <- data.frame(
    query = query,
    zh_variable = zh$zh_variable[keep],
    zh_path = zh$zh_path[keep],
    zh_score = score[keep],
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$zh_score, out$zh_path), , drop = FALSE]
  head(out, max_results)
}

.ukb_dictionary_translation_map <- function(translation_map = NULL) {
  base <- setNames(
    c(
      "sex",
      "age",
      "age at recruitment",
      "body mass index",
      "body mass index",
      "standing height",
      "weight",
      "waist circumference",
      "hip circumference",
      "systolic blood pressure",
      "diastolic blood pressure",
      "smoking",
      "smoking",
      "alcohol",
      "qualifications education",
      "ethnic background",
      "Townsend deprivation index",
      "asthma",
      "diabetes",
      "chronic obstructive pulmonary disease COPD",
      "hypertension",
      "date first reported",
      "source of report",
      "diagnosed diagnosis",
      "date",
      "source",
      "Instance 0",
      "Instance 1",
      "Instance 2",
      "Instance 3"
    ),
    c(
      "\u6027\u522b",
      "\u5e74\u9f84",
      "\u62db\u52df\u5e74\u9f84",
      "\u8eab\u4f53\u8d28\u91cf\u6307\u6570",
      "\u4f53\u91cd\u6307\u6570",
      "\u8eab\u9ad8",
      "\u4f53\u91cd",
      "\u8170\u56f4",
      "\u81c0\u56f4",
      "\u6536\u7f29\u538b",
      "\u8212\u5f20\u538b",
      "\u5438\u70df",
      "\u62bd\u70df",
      "\u996e\u9152",
      "\u6559\u80b2",
      "\u79cd\u65cf",
      "\u591a\u91cd\u5265\u593a\u6307\u6570",
      "\u54ee\u5598",
      "\u7cd6\u5c3f\u75c5",
      "\u6162\u6027\u963b\u585e\u6027\u80ba\u75be\u75c5",
      "\u9ad8\u8840\u538b",
      "\u9996\u6b21\u62a5\u544a\u7684\u65e5\u671f",
      "\u62a5\u544a\u6765\u6e90",
      "\u8bca\u65ad",
      "\u65e5\u671f",
      "\u6765\u6e90",
      "\u5b9e\u4f8b0",
      "\u5b9e\u4f8b1",
      "\u5b9e\u4f8b2",
      "\u5b9e\u4f8b3"
    )
  )
  if (is.null(translation_map)) {
    return(base)
  }
  extra <- if (is.data.frame(translation_map)) {
    if (!all(c("zh", "en") %in% names(translation_map))) {
      stop("'translation_map' data.frame must contain columns 'zh' and 'en'.", call. = FALSE)
    }
    setNames(as.character(translation_map$en), as.character(translation_map$zh))
  } else if (is.character(translation_map) && !is.null(names(translation_map))) {
    translation_map
  } else {
    stop("'translation_map' must be NULL, a named character vector, or a data.frame with zh/en columns.", call. = FALSE)
  }
  c(base, extra)
}

.ukb_make_english_candidate <- function(query, zh_match, trans) {
  terms <- c(query)
  if (nrow(zh_match) > 0L) {
    terms <- c(terms, zh_match$zh_variable[[1]], zh_match$zh_path[[1]])
  }
  terms <- unique(terms[nzchar(terms)])

  candidates <- character()
  for (zh_term in names(trans)) {
    if (any(grepl(zh_term, terms, fixed = TRUE))) {
      candidates <- c(candidates, unname(trans[[zh_term]]))
    }
  }

  codes <- unique(unlist(regmatches(terms, gregexpr("\\b[A-Z][0-9]{2}(?:\\.[0-9A-Z]+)?\\b", terms, perl = TRUE))))
  if (length(codes) > 0L) {
    candidates <- c(candidates, codes, paste(codes, "first reported"), paste("source of report", codes))
  }

  if (length(candidates) == 0L) {
    candidates <- terms
  }
  paste(unique(candidates), collapse = " ")
}

.ukb_search_official_dictionary <- function(query,
                                            english_candidate,
                                            official,
                                            query_language,
                                            max_results,
                                            min_score) {
  if (query_language == "field_id") {
    score <- ifelse(official$FieldID == query, 1, 0)
  } else if (query_language == "column") {
    q <- sub("^participant\\.", "", query)
    cols <- sub("^participant\\.", "", official$UKBColumn)
    score <- ifelse(cols == q, 1, .ukb_text_match_score(q, cols))
  } else {
    score <- .ukb_text_match_score(english_candidate, official$search_text)
  }
  keep <- score >= min_score
  if (!any(keep)) {
    return(.ukb_empty_official_result())
  }
  out <- official[keep, c("FieldID", "UKBColumn", "Field", "Title", "Category",
                          "ValueType", "Units", "Coding", "Instance", "Array", "Link"), drop = FALSE]
  out$query <- query
  out$english_candidate <- english_candidate
  out$match_score <- score[keep]
  out$match_method <- if (query_language == "zh") "zh_dictionary_rule_translation" else query_language
  out$review_status <- ifelse(out$match_score >= 0.8, "high_confidence", "needs_review")
  out <- out[order(-out$match_score, out$FieldID, out$UKBColumn), , drop = FALSE]
  out <- head(out, max_results)
  out[, c("query", "english_candidate", "FieldID", "UKBColumn", "Field", "Title",
          "Category", "ValueType", "Units", "Coding", "Instance", "Array", "Link",
          "match_method", "match_score", "review_status")]
}

.ukb_empty_official_result <- function() {
  data.frame(
    query = character(0),
    english_candidate = character(0),
    FieldID = character(0),
    UKBColumn = character(0),
    Field = character(0),
    Title = character(0),
    Category = character(0),
    ValueType = character(0),
    Units = character(0),
    Coding = character(0),
    Instance = character(0),
    Array = character(0),
    Link = character(0),
    match_method = character(0),
    match_score = numeric(0),
    review_status = character(0),
    stringsAsFactors = FALSE
  )
}

.ukb_text_match_score <- function(query, text) {
  q <- tolower(trimws(query))
  text_lower <- tolower(trimws(text))
  text_lower[is.na(text_lower)] <- ""
  exact <- text_lower == q
  contains <- grepl(q, text_lower, fixed = TRUE)

  q_terms <- unique(strsplit(gsub("[^a-z0-9]+", " ", q), "\\s+")[[1]])
  q_terms <- q_terms[nzchar(q_terms)]
  term_score <- rep(0, length(text_lower))
  if (length(q_terms) > 0L) {
    term_hits <- vapply(q_terms, function(term) grepl(term, text_lower, fixed = TRUE), logical(length(text_lower)))
    if (is.null(dim(term_hits))) {
      term_hits <- matrix(term_hits, ncol = 1L)
    }
    term_score <- rowMeans(term_hits)
  }

  approx_score <- rep(0, length(text_lower))
  candidate <- substr(text_lower, 1, pmin(nchar(text_lower), max(nchar(q), 1L) + 20L))
  ok <- nzchar(q) & nzchar(candidate)
  if (isTRUE(ok)) {
    dist <- adist(q, candidate)
    approx_score <- pmax(0, 1 - as.numeric(dist) / pmax(nchar(q), nchar(candidate), 1L))
  }

  pmax(
    ifelse(exact, 1, 0),
    ifelse(contains, 0.9, 0),
    0.75 * term_score,
    0.55 * approx_score
  )
}

.ukb_normalize_column_names <- function(x, ignore_entity_prefix) {
  x <- trimws(as.character(x))
  if (isTRUE(ignore_entity_prefix)) {
    x <- sub("^participant\\.", "", x)
  }
  x
}
