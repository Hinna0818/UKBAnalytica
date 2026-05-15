#' Get structured UK Biobank field metadata
#'
#' @description
#' Returns a structured data.frame of UK Biobank field metadata. When
#' `ukb_data_dict` is supplied, the function reads a UK Biobank data dictionary
#' metadata file available in the current session and standardizes common
#' metadata columns. When `fields_df` or
#' a RAP dataset is supplied, the function also records the approved RAP field
#' names available in the current project.
#'
#' This is intended to be a simple entry point for users who want to inspect UKB
#' field metadata in R before planning an extraction.
#'
#' @param field_id Optional UKB numeric field IDs to keep.
#' @param query Optional keyword used to filter the metadata table. The keyword
#'   is matched against the title, description, category, and RAP field names.
#' @param ukb_data_dict Optional path to a `Data_Dictionary_Showcase.tsv`
#'   file or equivalent UKB metadata export available in the current session.
#' @param dataset Optional RAP `.dataset` file name. Used only when `fields_df`
#'   is `NULL` and RAP field metadata should be retrieved live.
#' @param fields_df Optional data.frame returned by `rap_list_fields()`. This is
#'   useful for offline testing or when you already cached the RAP field list.
#' @param entity RAP dataset entity. Defaults to `"participant"`.
#'
#' @return A data.frame with one row per UKB field and standardized metadata
#'   columns. When RAP field metadata is available, the result also includes the
#'   matching RAP column names and the number of approved RAP columns per field.
#'
#' @examples
#' \dontrun{
#' # Use a UKB data dictionary metadata file available in the current session
#' meta <- get_field_metadata(
#'   query = "blood pressure",
#'   ukb_data_dict = "Data_Dictionary_Showcase.tsv"
#' )
#'
#' # Combine the data dictionary with the current RAP-approved field list
#' meta <- get_field_metadata(
#'   field_id = c(31, 53, 21022, 4080),
#'   ukb_data_dict = "Data_Dictionary_Showcase.tsv"
#' )
#' }
#'
#' @export
get_field_metadata <- function(field_id = NULL,
                               query = NULL,
                               ukb_data_dict = NULL,
                               dataset = NULL,
                               fields_df = NULL,
                               entity = "participant") {
  if (is.null(ukb_data_dict) && is.null(fields_df) && is.null(dataset)) {
    stop(
      "Provide `ukb_data_dict`, `fields_df`, or a RAP `dataset` to retrieve field metadata.",
      call. = FALSE
    )
  }

  id_filter <- .field_metadata_normalize_ids(field_id)

  showcase_meta <- NULL
  if (!is.null(ukb_data_dict)) {
    showcase_meta <- .field_metadata_read_showcase(ukb_data_dict)
  }

  rap_meta <- NULL
  if (!is.null(fields_df) || is.null(ukb_data_dict) || !is.null(dataset)) {
    if (is.null(fields_df)) {
      fields_df <- rap_list_fields(dataset = dataset, entity = entity)
    }
    rap_meta <- .field_metadata_from_rap(fields_df, entity = entity)
  }

  out <- .field_metadata_merge(showcase_meta, rap_meta)
  out <- .field_metadata_filter(out, field_id = id_filter, query = query)

  rownames(out) <- NULL
  out
}

#' Get one UK Biobank field's metadata
#'
#' @description
#' Convenience wrapper around `get_field_metadata()` for a single UKB `field_id`.
#' This is the simplest way to ask “what does field 4080 correspond to?” and get
#' a one-row metadata table back in R.
#'
#' @param field_id A single UKB numeric field ID.
#' @param ukb_data_dict Optional path to a `Data_Dictionary_Showcase.tsv`
#'   file or equivalent UKB metadata export available in the current session.
#' @param dataset Optional RAP `.dataset` file name. Used only when `fields_df`
#'   is `NULL` and RAP field metadata should be retrieved live.
#' @param fields_df Optional data.frame returned by `rap_list_fields()`. This is
#'   useful for offline testing or when you already cached the RAP field list.
#' @param entity RAP dataset entity. Defaults to `"participant"`.
#' @param live Logical. If `TRUE`, fetch the field page from the public UK
#'   Biobank Showcase website and parse the displayed metadata for this
#'   `field_id`.
#' @param timeout Timeout in seconds used for the live web request.
#'
#' @return A one-row data.frame when the field is found.
#' @export
get_field_info <- function(field_id,
                           ukb_data_dict = NULL,
                           dataset = NULL,
                           fields_df = NULL,
                           entity = "participant",
                           live = FALSE,
                           timeout = 30) {
  normalized_id <- .field_metadata_normalize_ids(field_id)
  if (length(normalized_id) != 1L) {
    stop("`field_id` must be a single UKB numeric field ID.", call. = FALSE)
  }

  out <- NULL
  if (!is.null(ukb_data_dict) || !is.null(dataset) || !is.null(fields_df)) {
    out <- get_field_metadata(
      field_id = normalized_id,
      ukb_data_dict = ukb_data_dict,
      dataset = dataset,
      fields_df = fields_df,
      entity = entity
    )
  }

  if (isTRUE(live)) {
    live_info <- .field_metadata_fetch_live(normalized_id, timeout = timeout)
    out <- .field_metadata_overlay_one(out, live_info)
  }

  if (is.null(out) || nrow(out) == 0) {
    stop(sprintf("No metadata found for field_id %s.", normalized_id), call. = FALSE)
  }

  out[1, , drop = FALSE]
}

.field_metadata_empty <- function() {
  data.frame(
    field_id = integer(0),
    title = character(0),
    description = character(0),
    category = character(0),
    value_type = character(0),
    units = character(0),
    coding = character(0),
    instances = character(0),
    array = character(0),
    participants = character(0),
    item_count = character(0),
    sexed = character(0),
    debut = character(0),
    version = character(0),
    stability = character(0),
    strata = character(0),
    cost_tier = character(0),
    notes = character(0),
    link = character(0),
    field_name = character(0),
    rap_field_names = character(0),
    n_rap_columns = integer(0),
    source_backend = character(0),
    stringsAsFactors = FALSE
  )
}

.field_metadata_normalize_ids <- function(field_id) {
  if (is.null(field_id)) {
    return(NULL)
  }

  field_id <- as.character(field_id)
  field_id <- trimws(field_id[nzchar(field_id)])
  if (length(field_id) == 0 || anyNA(field_id) || any(!grepl("^[0-9]+$", field_id))) {
    stop("`field_id` must contain UKB numeric field IDs, e.g. c(31, 53, 21022).", call. = FALSE)
  }

  as.integer(unique(field_id))
}

.field_metadata_pick_col <- function(df, aliases) {
  nm <- names(df)
  clean <- gsub("[^a-z0-9]+", "_", tolower(nm))
  alias_clean <- gsub("[^a-z0-9]+", "_", tolower(aliases))
  hit <- match(alias_clean, clean)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) {
    return(NULL)
  }
  nm[[hit[[1]]]]
}

.field_metadata_char_col <- function(df, col) {
  if (is.null(col) || !col %in% names(df)) {
    return(rep(NA_character_, nrow(df)))
  }
  as.character(df[[col]])
}

.field_metadata_read_showcase <- function(path) {
  if (!file.exists(path)) {
    stop("`ukb_data_dict` does not exist: ", path, call. = FALSE)
  }

  raw <- data.table::fread(path, data.table = FALSE, encoding = "UTF-8")
  if (!is.data.frame(raw) || nrow(raw) == 0) {
    stop("`ukb_data_dict` could not be parsed into a non-empty table.", call. = FALSE)
  }

  field_id_col <- .field_metadata_pick_col(raw, c("FieldID", "Field ID", "field_id"))
  title_col <- .field_metadata_pick_col(raw, c("Field", "Title", "field"))

  if (is.null(field_id_col) || is.null(title_col)) {
    stop(
      "`ukb_data_dict` must include columns equivalent to `FieldID` and `Field`.",
      call. = FALSE
    )
  }

  out <- data.frame(
    field_id = suppressWarnings(as.integer(raw[[field_id_col]])),
    title = .field_metadata_char_col(raw, title_col),
    description = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Description", "Notes", "description"))),
    category = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Category", "Main Category", "Path", "Category Path"))),
    value_type = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("ValueType", "Value Type", "value_type"))),
    units = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Units", "Unit", "units"))),
    coding = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Coding", "Coding ID", "coding"))),
    instances = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Instances", "Instance", "instances"))),
    array = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Array", "Arrays", "array"))),
    participants = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Participants", "participants"))),
    notes = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Notes", "notes"))),
    link = .field_metadata_char_col(raw, .field_metadata_pick_col(raw, c("Link", "URL", "link"))),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$field_id), , drop = FALSE]
  if (nrow(out) == 0) {
    stop("No valid UKB field IDs were found in `ukb_data_dict`.", call. = FALSE)
  }

  keep <- !duplicated(out$field_id)
  out <- out[keep, , drop = FALSE]
  out$field_name <- NA_character_
  out$rap_field_names <- NA_character_
  out$n_rap_columns <- NA_integer_
  out$item_count <- NA_character_
  out$sexed <- NA_character_
  out$debut <- NA_character_
  out$version <- NA_character_
  out$stability <- NA_character_
  out$strata <- NA_character_
  out$cost_tier <- NA_character_
  out$source_backend <- "showcase"
  out
}

.field_metadata_from_rap <- function(fields_df,
                                     entity = "participant") {
  if (!is.data.frame(fields_df) || !all(c("field_name", "title") %in% names(fields_df))) {
    stop("`fields_df` must be a data.frame with columns `field_name` and `title`.", call. = FALSE)
  }

  dt <- as.data.frame(fields_df, stringsAsFactors = FALSE)
  field_pattern <- paste0("^", entity, "\\.p[0-9]+(?:_|$)")
  id_pattern <- paste0("^", entity, "\\.p([0-9]+)(?:_|$).*")
  dt <- dt[grepl(field_pattern, dt$field_name, perl = TRUE), , drop = FALSE]

  if (nrow(dt) == 0) {
    return(.field_metadata_empty())
  }

  dt$field_id <- suppressWarnings(as.integer(sub(id_pattern, "\\1", dt$field_name, perl = TRUE)))
  dt <- dt[!is.na(dt$field_id), , drop = FALSE]
  dt$title_clean <- trimws(sub("\\s*\\|.*$", "", dt$title))

  split_rows <- split(dt, dt$field_id)
  out <- do.call(rbind, lapply(split_rows, function(x) {
    data.frame(
      field_id = x$field_id[[1]],
      field_name = x$field_name[[1]],
      title = x$title_clean[[1]],
      rap_field_names = paste(unique(x$field_name), collapse = ";"),
      n_rap_columns = length(unique(x$field_name)),
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out$description <- NA_character_
  out$category <- NA_character_
  out$value_type <- NA_character_
  out$units <- NA_character_
  out$coding <- NA_character_
  out$instances <- NA_character_
  out$array <- NA_character_
  out$participants <- NA_character_
  out$item_count <- NA_character_
  out$sexed <- NA_character_
  out$debut <- NA_character_
  out$version <- NA_character_
  out$stability <- NA_character_
  out$strata <- NA_character_
  out$cost_tier <- NA_character_
  out$notes <- NA_character_
  out$link <- NA_character_
  out$source_backend <- "rap"

  out[, c(
    "field_id", "title", "description", "category", "value_type",
    "units", "coding", "instances", "array", "participants",
    "item_count", "sexed", "debut", "version", "stability", "strata",
    "cost_tier", "notes", "link", "field_name", "rap_field_names",
    "n_rap_columns", "source_backend"
  )]
}

.field_metadata_merge <- function(showcase_meta,
                                  rap_meta) {
  if (is.null(showcase_meta) && is.null(rap_meta)) {
    return(.field_metadata_empty())
  }
  if (is.null(showcase_meta)) {
    return(rap_meta)
  }
  if (is.null(rap_meta) || nrow(rap_meta) == 0) {
    return(showcase_meta)
  }

  merged <- merge(showcase_meta, rap_meta, by = "field_id", all = TRUE, suffixes = c(".showcase", ".rap"))

  pick <- function(showcase_col, rap_col = NULL) {
    showcase_vals <- if (showcase_col %in% names(merged)) merged[[showcase_col]] else NULL
    rap_vals <- if (!is.null(rap_col) && rap_col %in% names(merged)) merged[[rap_col]] else NULL

    if (is.null(rap_vals)) {
      return(showcase_vals)
    }

    out <- showcase_vals
    missing_idx <- is.na(out) | (!nzchar(out) & !is.na(rap_vals))
    out[missing_idx] <- rap_vals[missing_idx]
    out
  }

  out <- data.frame(
    field_id = merged$field_id,
    title = pick("title.showcase", "title.rap"),
    description = pick("description.showcase"),
    category = pick("category.showcase"),
    value_type = pick("value_type.showcase"),
    units = pick("units.showcase"),
    coding = pick("coding.showcase"),
    instances = pick("instances.showcase"),
    array = pick("array.showcase"),
    participants = pick("participants.showcase"),
    item_count = pick("item_count.showcase"),
    sexed = pick("sexed.showcase"),
    debut = pick("debut.showcase"),
    version = pick("version.showcase"),
    stability = pick("stability.showcase"),
    strata = pick("strata.showcase"),
    cost_tier = pick("cost_tier.showcase"),
    notes = pick("notes.showcase"),
    link = pick("link.showcase"),
    field_name = pick("field_name.showcase", "field_name.rap"),
    rap_field_names = pick("rap_field_names.showcase", "rap_field_names.rap"),
    n_rap_columns = if ("n_rap_columns.rap" %in% names(merged)) merged[["n_rap_columns.rap"]] else NA_integer_,
    source_backend = ifelse(
      !is.na(merged$source_backend.showcase) & !is.na(merged$source_backend.rap),
      "showcase+rap",
      ifelse(!is.na(merged$source_backend.showcase), "showcase", "rap")
    ),
    stringsAsFactors = FALSE
  )

  out
}

.field_metadata_match_query <- function(x,
                                        query) {
  if (is.na(x) || !nzchar(x)) {
    return(FALSE)
  }
  x <- tolower(as.character(x))
  query <- tolower(as.character(query))
  any(vapply(query, function(q) grepl(q, x, fixed = TRUE), logical(1)))
}

.field_metadata_filter <- function(x,
                                   field_id = NULL,
                                   query = NULL) {
  if (nrow(x) == 0) {
    return(x)
  }

  if (!is.null(field_id)) {
    x <- x[x$field_id %in% field_id, , drop = FALSE]
  }

  if (!is.null(query)) {
    if (!is.character(query) || length(query) == 0 || anyNA(query)) {
      stop("`query` must be a non-empty character vector.", call. = FALSE)
    }
    query <- trimws(query[nzchar(query)])
    keep <- vapply(seq_len(nrow(x)), function(i) {
      row <- x[i, , drop = FALSE]
      any(vapply(
        c("title", "description", "category", "field_name", "rap_field_names"),
        function(col) .field_metadata_match_query(row[[col]], query),
        logical(1)
      ))
    }, logical(1))
    x <- x[keep, , drop = FALSE]
  }

  x <- x[order(x$field_id), , drop = FALSE]
  rownames(x) <- NULL
  x
}

.field_metadata_require_xml2 <- function() {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop(
      "Live UKB field lookup requires the `xml2` package. Install it or use `ukb_data_dict` / `fields_df` instead.",
      call. = FALSE
    )
  }
}

.field_metadata_parse_value_pairs <- function(cells) {
  cells <- trimws(cells)
  cells <- cells[nzchar(cells)]
  if (length(cells) < 2) {
    return(list())
  }

  out <- list()
  idx <- seq.int(1L, length(cells) - 1L, by = 2L)
  for (i in idx) {
    key <- sub(":$", "", cells[[i]])
    val <- cells[[i + 1L]]
    if (nzchar(key) && nzchar(val)) {
      out[[key]] <- val
    }
  }
  out
}

.field_metadata_compact_ws <- function(x) {
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

.field_metadata_parse_live_html <- function(html,
                                            field_id,
                                            url) {
  .field_metadata_require_xml2()

  doc <- xml2::read_html(html)
  rows <- xml2::xml_find_all(doc, ".//tr")
  pairs <- list()

  for (row in rows) {
    cells <- xml2::xml_text(xml2::xml_find_all(row, ".//th|.//td"), trim = TRUE)
    row_pairs <- .field_metadata_parse_value_pairs(cells)
    if (length(row_pairs) > 0) {
      pairs[names(row_pairs)] <- row_pairs
    }
  }

  page_text <- .field_metadata_compact_ws(xml2::xml_text(doc))
  title <- unname(pairs[["Description"]])
  category <- unname(pairs[["Category"]])
  participants <- unname(pairs[["Participants"]])
  item_count <- unname(pairs[["Item count"]])
  value_type_raw <- unname(pairs[["Value Type"]])
  instances <- unname(pairs[["Instances"]])
  array <- unname(pairs[["Array"]])
  sexed <- unname(pairs[["Sexed"]])
  debut <- unname(pairs[["Debut"]])
  version <- unname(pairs[["Version"]])
  stability <- unname(pairs[["Stability"]])
  strata <- unname(pairs[["Strata"]])
  cost_tier <- unname(pairs[["Cost Tier"]])

  units <- NA_character_
  value_type <- value_type_raw
  if (!is.null(value_type_raw) && nzchar(value_type_raw) && grepl(",", value_type_raw, fixed = TRUE)) {
    parts <- strsplit(value_type_raw, ",", fixed = TRUE)[[1]]
    value_type <- trimws(parts[[1]])
    units <- trimws(paste(parts[-1], collapse = ","))
  }

  if ((is.na(units) || !nzchar(units)) && grepl("Units of measurement are ", page_text, fixed = TRUE)) {
    units <- sub(".*Units of measurement are ([^.]+)\\..*", "\\1", page_text)
    units <- trimws(units)
  }

  notes <- NA_character_
  if (grepl("Defined-instances run from ", page_text, fixed = TRUE)) {
    notes <- sub(
      ".*(Defined-instances run from .*?)(Statistics|Enabling scientific discoveries.*|$)",
      "\\1",
      page_text
    )
    notes <- .field_metadata_compact_ws(notes)
  }

  data.frame(
    field_id = as.integer(field_id),
    title = if (!is.null(title) && nzchar(title)) title else NA_character_,
    description = NA_character_,
    category = if (!is.null(category) && nzchar(category)) category else NA_character_,
    value_type = if (!is.null(value_type) && nzchar(value_type)) value_type else NA_character_,
    units = if (!is.null(units) && nzchar(units)) units else NA_character_,
    coding = NA_character_,
    instances = if (!is.null(instances) && nzchar(instances)) instances else NA_character_,
    array = if (!is.null(array) && nzchar(array)) array else NA_character_,
    participants = if (!is.null(participants) && nzchar(participants)) participants else NA_character_,
    item_count = if (!is.null(item_count) && nzchar(item_count)) item_count else NA_character_,
    sexed = if (!is.null(sexed) && nzchar(sexed)) sexed else NA_character_,
    debut = if (!is.null(debut) && nzchar(debut)) debut else NA_character_,
    version = if (!is.null(version) && nzchar(version)) version else NA_character_,
    stability = if (!is.null(stability) && nzchar(stability)) stability else NA_character_,
    strata = if (!is.null(strata) && nzchar(strata)) strata else NA_character_,
    cost_tier = if (!is.null(cost_tier) && nzchar(cost_tier)) cost_tier else NA_character_,
    notes = if (!is.null(notes) && nzchar(notes)) notes else NA_character_,
    link = url,
    field_name = NA_character_,
    rap_field_names = NA_character_,
    n_rap_columns = NA_integer_,
    source_backend = "showcase_live",
    stringsAsFactors = FALSE
  )
}

.field_metadata_fetch_live <- function(field_id,
                                       timeout = 30) {
  .field_metadata_require_xml2()

  url <- sprintf("https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=%s", as.integer(field_id))
  old_timeout <- getOption("timeout")
  options(timeout = max(as.integer(timeout), old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  doc <- tryCatch(
    xml2::read_html(url),
    error = function(e) {
      stop("Failed to retrieve UKB field page: ", conditionMessage(e), call. = FALSE)
    }
  )

  .field_metadata_parse_live_html(as.character(doc), field_id = field_id, url = url)
}

.field_metadata_overlay_one <- function(base,
                                        extra) {
  if (is.null(base) || nrow(base) == 0) {
    return(extra)
  }
  if (is.null(extra) || nrow(extra) == 0) {
    return(base)
  }

  all_cols <- union(names(base), names(extra))
  for (col in setdiff(all_cols, names(base))) {
    base[[col]] <- NA
  }
  for (col in setdiff(all_cols, names(extra))) {
    extra[[col]] <- NA
  }

  base <- base[, all_cols, drop = FALSE]
  extra <- extra[, all_cols, drop = FALSE]

  out <- base[1, , drop = FALSE]
  for (col in setdiff(all_cols, c("field_id", "source_backend"))) {
    cur <- out[[col]][[1]]
    new <- extra[[col]][[1]]
    if (is.na(cur) || (is.character(cur) && !nzchar(cur))) {
      out[[col]] <- new
    }
  }

  backends <- unique(c(as.character(base$source_backend[[1]]), as.character(extra$source_backend[[1]])))
  backends <- backends[!is.na(backends) & nzchar(backends)]
  out$source_backend <- paste(backends, collapse = "+")
  out
}
