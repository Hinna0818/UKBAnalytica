#' Load the bundled UK Biobank non-ratio metabolite panel
#'
#' @description
#' Read the metabolite metadata table bundled in `inst/extdata`. The current
#' file contains UK Biobank Nightingale non-ratio metabolite names, field IDs,
#' and RAP-style column names. It is mainly intended as a helper for examples,
#' tests, and metabolite-name checking.
#'
#' @param file Optional path to a metabolite panel file. If `NULL`, the bundled
#'   `metabolites_non_ratio.txt` file is used.
#' @param file_encoding Character file encoding. Default is `"UTF-16LE"` because
#'   the bundled table is UTF-16 little-endian encoded.
#'
#' @return A data.frame with columns such as `Description`, `UKB_ID`, and
#'   `meta_ID`.
#'
#' @examples
#' panel <- load_ukb_metabolite_panel()
#' head(panel)
#'
#' @export
load_ukb_metabolite_panel <- function(file = NULL,
                                      file_encoding = "UTF-16LE") {
  if (is.null(file)) {
    file <- system.file("extdata", "metabolites_non_ratio.txt", package = "UKBAnalytica")
    if (!nzchar(file)) {
      file <- file.path("inst", "extdata", "metabolites_non_ratio.txt")
    }
  }
  if (!file.exists(file)) {
    stop("Metabolite panel file not found: ", file)
  }

  utils::read.delim(
    file,
    sep = "\t",
    fileEncoding = file_encoding,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


#' Classify UK Biobank metabolite names
#'
#' @description
#' Classify metabolite-like names into broad groups used by the metabolomics
#' ORA workflow. Small molecules can be mapped to MetaboAnalyst-compatible names,
#' whereas lipoprotein subclass measures, lipid aggregate measures, and proteins
#' are retained in the mapping table but are not passed to small-molecule ORA by
#' default.
#'
#' @param metabolites Character vector of metabolite names.
#'
#' @return A data.frame with `metabolite`, `category`, and `metaboanalyst_name`.
#'
#' @examples
#' classify_metabolites(c("Alanine", "LDL Cholesterol", "Apolipoprotein B"))
#'
#' @export
classify_metabolites <- function(metabolites) {
  metabolites <- .metab_clean_input(metabolites)
  mapped <- metabolite_to_metaboanalyst_name(metabolites, drop_unmapped = FALSE)

  category <- vapply(metabolites, .metab_classify_one, character(1))
  mapped$category <- category[match(mapped$metabolite, metabolites)]
  mapped[, c("metabolite", "category", "metaboanalyst_name", "mapping_source")]
}


#' Map metabolite names to MetaboAnalyst-compatible names
#'
#' @description
#' Convert common UK Biobank Nightingale metabolite labels to names that are more
#' likely to be recognized by MetaboAnalystR name cross-referencing. Users can
#' provide a custom mapping table to override or extend the built-in map.
#'
#' @param metabolites Character vector of metabolite names.
#' @param mapping_table Optional data.frame with metabolite-to-name mappings.
#' @param mapping_metabolite_col Column in `mapping_table` containing input
#'   metabolite labels. Default `"metabolite"`.
#' @param mapping_name_col Column in `mapping_table` containing mapped names.
#'   Default `"metaboanalyst_name"`.
#' @param drop_unmapped Logical. If `TRUE`, keep only mapped rows. Default
#'   `FALSE`.
#'
#' @return A data.frame with `metabolite`, `metaboanalyst_name`, and
#'   `mapping_source`.
#'
#' @examples
#' metabolite_to_metaboanalyst_name(c("Acetate", "Alanine", "LDL Cholesterol"))
#'
#' @export
metabolite_to_metaboanalyst_name <- function(metabolites,
                                             mapping_table = NULL,
                                             mapping_metabolite_col = "metabolite",
                                             mapping_name_col = "metaboanalyst_name",
                                             drop_unmapped = FALSE) {
  metabolites <- .metab_clean_input(metabolites)
  if (!length(metabolites)) {
    stop("No valid metabolite names were provided.")
  }

  result <- data.frame(
    metabolite = metabolites,
    metaboanalyst_name = NA_character_,
    mapping_source = NA_character_,
    stringsAsFactors = FALSE
  )

  built_in <- .metab_builtin_name_map()
  key <- .metab_normalize_key(metabolites)
  matched <- match(key, .metab_normalize_key(names(built_in)))
  has_match <- !is.na(matched)
  if (any(has_match)) {
    result$metaboanalyst_name[has_match] <- unname(built_in[matched[has_match]])
    result$mapping_source[has_match] <- "built_in"
  }

  if (!is.null(mapping_table)) {
    if (!is.data.frame(mapping_table)) {
      stop("`mapping_table` must be a data.frame.")
    }
    required <- c(mapping_metabolite_col, mapping_name_col)
    if (!all(required %in% names(mapping_table))) {
      stop("`mapping_table` must contain columns: ", paste(required, collapse = ", "))
    }
    custom <- data.frame(
      metabolite = .metab_clean_input(mapping_table[[mapping_metabolite_col]]),
      metaboanalyst_name = .metab_clean_input(mapping_table[[mapping_name_col]]),
      stringsAsFactors = FALSE
    )
    custom <- custom[!is.na(custom$metabolite) & nzchar(custom$metabolite), , drop = FALSE]
    custom <- custom[!duplicated(.metab_normalize_key(custom$metabolite)), , drop = FALSE]
    custom_idx <- match(key, .metab_normalize_key(custom$metabolite))
    custom_match <- !is.na(custom_idx)
    if (any(custom_match)) {
      result$metaboanalyst_name[custom_match] <- custom$metaboanalyst_name[custom_idx[custom_match]]
      result$mapping_source[custom_match] <- "custom_table"
    }
  }

  if (isTRUE(drop_unmapped)) {
    result <- result[!is.na(result$metaboanalyst_name) & nzchar(result$metaboanalyst_name), , drop = FALSE]
  }
  rownames(result) <- NULL
  result
}


#' Run metabolite over-representation analysis
#'
#' @description
#' Run ORA for a metabolite list. The recommended first backend is
#' `backend = "custom"`, where users provide a two-column metabolite pathway
#' library. A `backend = "metaboanalyst"` interface is also provided for users
#' who have installed MetaboAnalystR and want to use its metabolite-set
#' libraries, such as `"smpdb_pathway"`.
#'
#' @param metabolites Character vector of metabolite names.
#' @param pathway_db Optional data.frame for custom ORA. Must contain pathway
#'   and metabolite columns.
#' @param universe Optional background metabolite vector. If `NULL`, the
#'   pathway library metabolites are used for custom ORA.
#' @param backend One of `"custom"` or `"metaboanalyst"`.
#' @param id_type Metabolite identifier type for MetaboAnalystR
#'   cross-referencing. Default `"name"`.
#' @param library MetaboAnalystR metabolite-set library. Default
#'   `"smpdb_pathway"`.
#' @param mapping_table Optional custom mapping table passed to
#'   `metabolite_to_metaboanalyst_name()`.
#' @param pathway_col Column name in `pathway_db` containing pathway names.
#'   Default `"pathway"`.
#' @param metabolite_col Column name in `pathway_db` containing metabolite
#'   names. Default `"metabolite"`.
#' @param min_metabolites Minimum mapped metabolites required for ORA.
#'   Default `3`.
#' @param p_adjust_method Multiple-testing method used by [stats::p.adjust()].
#'   Default `"BH"`.
#' @param run_subprocess Logical. For `backend = "metaboanalyst"`, run
#'   MetaboAnalystR in a clean subprocess to avoid global-state issues. Default
#'   `TRUE`.
#'
#' @return A list of class `ukb_metabolite_ora` with components `input`,
#'   `mapping`, `matched`, `unmatched`, `ora_result`, `backend`, and `library`.
#'
#' @examples
#' panel <- load_ukb_metabolite_panel()
#' hits <- c("Alanine", "Glutamine", "Glycine", "Lactate", "Pyruvate")
#' pathway_db <- data.frame(
#'   pathway = c(rep("Amino acid metabolism", 3), rep("Energy metabolism", 2)),
#'   metabolite = c("L-Alanine", "L-Glutamine", "Glycine", "Lactic acid", "Pyruvic acid")
#' )
#' run_metabolite_ora(hits, pathway_db = pathway_db, backend = "custom")
#'
#' @export
run_metabolite_ora <- function(metabolites,
                               pathway_db = NULL,
                               universe = NULL,
                               backend = c("custom", "metaboanalyst"),
                               id_type = "name",
                               library = "smpdb_pathway",
                               mapping_table = NULL,
                               pathway_col = "pathway",
                               metabolite_col = "metabolite",
                               min_metabolites = 3,
                               p_adjust_method = "BH",
                               run_subprocess = TRUE) {
  backend <- match.arg(backend)
  mapping <- metabolite_to_metaboanalyst_name(
    metabolites,
    mapping_table = mapping_table,
    drop_unmapped = FALSE
  )
  categories <- classify_metabolites(metabolites)
  mapping$category <- categories$category[match(mapping$metabolite, categories$metabolite)]

  matched <- unique(stats::na.omit(mapping$metaboanalyst_name[mapping$category == "small_molecule"]))
  unmatched <- mapping$metabolite[is.na(mapping$metaboanalyst_name) | mapping$category != "small_molecule"]

  if (length(matched) < min_metabolites) {
    stop(
      "Fewer than `min_metabolites` mapped small-molecule metabolites are available: ",
      length(matched)
    )
  }

  if (backend == "custom") {
    ora_df <- .metab_run_custom_ora(
      matched = matched,
      pathway_db = pathway_db,
      universe = universe,
      pathway_col = pathway_col,
      metabolite_col = metabolite_col,
      p_adjust_method = p_adjust_method
    )
  } else {
    ora_df <- .metab_run_metaboanalyst_ora(
      matched = matched,
      id_type = id_type,
      library = library,
      p_adjust_method = p_adjust_method,
      run_subprocess = run_subprocess
    )
  }

  result <- list(
    input = .metab_clean_input(metabolites),
    mapping = mapping,
    matched = matched,
    unmatched = unique(unmatched),
    ora_result = ora_df,
    backend = backend,
    library = if (backend == "metaboanalyst") library else "custom",
    p_adjust_method = p_adjust_method
  )
  class(result) <- "ukb_metabolite_ora"
  result
}


#' Plot metabolite ORA results as a dot plot
#'
#' @param x A data.frame returned by `run_metabolite_ora()$ora_result` or a
#'   `ukb_metabolite_ora` object.
#' @param top_n Number of pathways to show. Default `15`.
#' @param p_col P-value column used for ordering and color. Default `"pvalue"`.
#' @param size_col Column used for point size. Default `"hits"`.
#' @param pathway_col Column containing pathway names. Default `"pathway"`.
#' @param color_low,color_high Colors for the sequential p-value gradient.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn
#'   scale_size_continuous labs theme_classic theme element_blank element_line
#'   element_text
#' @importFrom rlang .data
#' @export
plot_metabolite_ora_dotplot <- function(x,
                                        top_n = 15,
                                        p_col = "pvalue",
                                        size_col = "hits",
                                        pathway_col = "pathway",
                                        color_low = "#2F6FA3",
                                        color_high = "#C74732") {
  plot_df <- .metab_extract_ora_df(x)
  plot_df <- .metab_prepare_plot_df(plot_df, top_n, p_col, pathway_col)
  plot_df$size_value <- as.numeric(plot_df[[size_col]])

  ggplot(plot_df, aes(x = .data$fold_enrichment, y = .data$pathway_label)) +
    geom_point(aes(size = .data$size_value, color = .data$neg_log10_p), alpha = 0.88) +
    scale_color_gradientn(
      colors = c(color_low, "#67A9CF", "#FDD49E", color_high),
      name = expression(-log[10](italic(P)))
    ) +
    scale_size_continuous(range = c(2.3, 7), name = "Hits") +
    labs(x = "Fold enrichment", y = NULL) +
    .metab_theme()
}


#' Plot metabolite ORA results as a bar plot
#'
#' @inheritParams plot_metabolite_ora_dotplot
#' @param fill_color Bar color.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_hline coord_flip labs
#'   theme_classic theme element_blank element_line element_text
#' @importFrom rlang .data
#' @export
plot_metabolite_ora_barplot <- function(x,
                                        top_n = 15,
                                        p_col = "pvalue",
                                        pathway_col = "pathway",
                                        fill_color = "#2F6FA3") {
  plot_df <- .metab_extract_ora_df(x)
  plot_df <- .metab_prepare_plot_df(plot_df, top_n, p_col, pathway_col)

  ggplot(plot_df, aes(x = .data$pathway_label, y = .data$neg_log10_p)) +
    geom_col(fill = fill_color, width = 0.68) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.32, colour = "#767676") +
    coord_flip() +
    labs(x = NULL, y = expression(-log[10](italic(P)))) +
    .metab_theme()
}


.metab_builtin_name_map <- function() {
  c(
    "Acetone" = "Acetone",
    "Acetate" = "Acetic acid",
    "Acetoacetate" = "Acetoacetic acid",
    "3-Hydroxybutyrate" = "3-Hydroxybutyric acid",
    "Lactate" = "Lactic acid",
    "Glucose-lactate" = "Lactic acid",
    "Citrate" = "Citric acid",
    "Glucose" = "D-Glucose",
    "Alanine" = "L-Alanine",
    "Spectrometer-corrected alanine" = "L-Alanine",
    "Valine" = "L-Valine",
    "Leucine" = "L-Leucine",
    "Isoleucine" = "L-Isoleucine",
    "Tyrosine" = "L-Tyrosine",
    "Phenylalanine" = "L-Phenylalanine",
    "Histidine" = "L-Histidine",
    "Glycine" = "Glycine",
    "Creatinine" = "Creatinine",
    "Linoleic Acid" = "Linoleic acid",
    "Linoleic acid" = "Linoleic acid",
    "Docosahexaenoic Acid" = "Docosahexaenoic acid",
    "Docosahexaenoic acid" = "Docosahexaenoic acid",
    "Glutamine" = "L-Glutamine",
    "Pyruvate" = "Pyruvic acid"
  )
}

.metab_classify_one <- function(name) {
  if (.metab_normalize_key(name) %in% .metab_normalize_key(names(.metab_builtin_name_map()))) {
    return("small_molecule")
  }

  if (grepl("Albumin|Apolipoprotein", name, ignore.case = TRUE)) {
    return("protein")
  }

  lipid_terms <- c(
    "VLDL", "LDL", "HDL", "IDL", "Lipoprotein", "Chylomicrons",
    "Cholesterol", "Triglyceride", "Triglycerides", "Phospholipid",
    "Lipid", "Sphingomyelin", "Cholines", "Fatty Acid", "Fatty Acids",
    "Omega", "Saturated", "Monounsaturated", "Polyunsaturated",
    "Unsaturation", "Remnant", "Glycoprotein", "Branched-Chain",
    "Concentration of", "Phosphoglyceride"
  )
  if (any(grepl(paste(lipid_terms, collapse = "|"), name, ignore.case = TRUE))) {
    return("lipoprotein_lipid")
  }

  "unknown"
}

.metab_run_custom_ora <- function(matched,
                                  pathway_db,
                                  universe,
                                  pathway_col,
                                  metabolite_col,
                                  p_adjust_method) {
  if (is.null(pathway_db)) {
    stop("`pathway_db` is required when `backend = \"custom\"`.")
  }
  if (!is.data.frame(pathway_db)) {
    stop("`pathway_db` must be a data.frame.")
  }
  required <- c(pathway_col, metabolite_col)
  if (!all(required %in% names(pathway_db))) {
    stop("`pathway_db` must contain columns: ", paste(required, collapse = ", "))
  }

  db <- data.frame(
    pathway = as.character(pathway_db[[pathway_col]]),
    metabolite = as.character(pathway_db[[metabolite_col]]),
    stringsAsFactors = FALSE
  )
  db <- db[!is.na(db$pathway) & nzchar(db$pathway) & !is.na(db$metabolite) & nzchar(db$metabolite), , drop = FALSE]
  db$key <- .metab_normalize_key(db$metabolite)
  db <- db[!duplicated(db[, c("pathway", "key")]), , drop = FALSE]

  query_key <- unique(.metab_normalize_key(matched))
  if (is.null(universe)) {
    universe_key <- unique(db$key)
  } else {
    universe_map <- metabolite_to_metaboanalyst_name(universe, drop_unmapped = FALSE)
    universe_key <- unique(.metab_normalize_key(stats::na.omit(universe_map$metaboanalyst_name)))
    universe_key <- unique(c(universe_key, query_key))
  }
  universe_key <- universe_key[nzchar(universe_key)]

  query_key <- intersect(query_key, universe_key)
  if (!length(query_key)) {
    stop("No query metabolites overlap with the ORA universe.")
  }

  M <- length(universe_key)
  N <- length(query_key)
  split_db <- split(db, db$pathway)
  rows <- lapply(names(split_db), function(pathway) {
    pathway_key <- intersect(unique(split_db[[pathway]]$key), universe_key)
    K <- length(pathway_key)
    hit_key <- intersect(query_key, pathway_key)
    x <- length(hit_key)
    if (K == 0 || x == 0) {
      return(NULL)
    }
    expected <- N * K / M
    data.frame(
      pathway = pathway,
      hits = x,
      pathway_size = K,
      query_size = N,
      universe_size = M,
      expected = expected,
      fold_enrichment = ifelse(expected > 0, x / expected, NA_real_),
      pvalue = stats::phyper(x - 1, K, M - K, N, lower.tail = FALSE),
      hit_names = paste(split_db[[pathway]]$metabolite[match(hit_key, split_db[[pathway]]$key)], collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) {
    out <- data.frame(
      pathway = character(),
      hits = integer(),
      pathway_size = integer(),
      query_size = integer(),
      universe_size = integer(),
      expected = numeric(),
      fold_enrichment = numeric(),
      pvalue = numeric(),
      hit_names = character(),
      stringsAsFactors = FALSE
    )
  }
  out$p_adjust <- stats::p.adjust(out$pvalue, method = p_adjust_method)
  out$neg_log10_p <- -log10(pmax(out$pvalue, .Machine$double.xmin))
  out <- out[order(out$pvalue, -out$hits), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.metab_run_metaboanalyst_ora <- function(matched,
                                         id_type,
                                         library,
                                         p_adjust_method,
                                         run_subprocess) {
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    stop(
      "Package 'MetaboAnalystR' is required for `backend = \"metaboanalyst\"`. ",
      "Install it from the xia-lab/MetaboAnalystR repository, then rerun."
    )
  }

  if (!isTRUE(run_subprocess)) {
    return(.metab_run_metaboanalyst_current_session(matched, id_type, library, p_adjust_method))
  }

  tmp_script <- tempfile(fileext = ".R")
  tmp_out <- tempfile(fileext = ".csv")
  tmp_dir <- tempfile("ukb_metaboanalyst_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    unlink(tmp_script)
    unlink(tmp_out)
    unlink(tmp_dir, recursive = TRUE)
  }, add = TRUE)

  script <- sprintf(
    'library(MetaboAnalystR)
setwd(%s)
.ukb_prepare_qs_file <- function(file, url) {
  if (!file.exists(file)) {
    utils::download.file(url, file, method = "libcurl", quiet = TRUE)
  }
  ok <- tryCatch({
    qs::qread(file)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop("MetaboAnalystR downloaded a qs2-format library file. Install qs2 and rerun.")
    }
    obj <- qs2::qs_read(file)
    qs::qsave(obj, file)
  }
  invisible(file)
}
.ukb_prepare_qs_file("compound_db.qs", "https://www.metaboanalyst.ca/resources/libs/compound_db.qs")
.ukb_prepare_qs_file("syn_nms.qs", "https://www.metaboanalyst.ca/resources/libs/syn_nms.qs")
.ukb_prepare_qs_file(%s, paste0("https://www.metaboanalyst.ca/resources/libs/msets/", %s))
metabolites <- c(%s)
mSet <- MetaboAnalystR::InitDataObjects("conc", "msetora", FALSE, default.dpi = 300)
mSet <- MetaboAnalystR::Setup.MapData(mSet, metabolites)
mSet <- MetaboAnalystR::CrossReferencing(mSet, %s)
mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)
mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE)
mSet <- MetaboAnalystR::SetCurrentMsetLib(mSet, %s, 0)
mSet <- MetaboAnalystR::CalculateHyperScore(mSet)
if (!is.null(mSet$analSet$ora.mat)) {
  out <- as.data.frame(mSet$analSet$ora.mat)
  out$pathway <- rownames(out)
  write.csv(out, %s, row.names = FALSE)
} else {
  write.csv(data.frame(), %s, row.names = FALSE)
}
',
    shQuote(normalizePath(tmp_dir, mustWork = FALSE)),
    shQuote(paste0(library, ".qs")),
    shQuote(paste0(library, ".qs")),
    paste(shQuote(matched), collapse = ", "),
    shQuote(id_type),
    shQuote(library),
    shQuote(tmp_out),
    shQuote(tmp_out)
  )
  writeLines(script, tmp_script)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(tmp_script), stdout = TRUE, stderr = TRUE)
  if (!file.exists(tmp_out)) {
    stop("MetaboAnalystR subprocess did not return an ORA table: ", paste(status, collapse = "\n"))
  }
  out <- utils::read.csv(tmp_out, stringsAsFactors = FALSE, check.names = FALSE)
  .metab_standardize_metaboanalyst_result(out, p_adjust_method)
}

.metab_run_metaboanalyst_current_session <- function(matched,
                                                     id_type,
                                                     library,
                                                     p_adjust_method) {
  .metab_prepare_metaboanalyst_qs_files(library)
  mSet <- getExportedValue("MetaboAnalystR", "InitDataObjects")("conc", "msetora", FALSE, default.dpi = 300)
  mSet <- getExportedValue("MetaboAnalystR", "Setup.MapData")(mSet, matched)
  mSet <- getExportedValue("MetaboAnalystR", "CrossReferencing")(mSet, id_type)
  mSet <- getExportedValue("MetaboAnalystR", "CreateMappingResultTable")(mSet)
  mSet <- getExportedValue("MetaboAnalystR", "SetMetabolomeFilter")(mSet, FALSE)
  mSet <- getExportedValue("MetaboAnalystR", "SetCurrentMsetLib")(mSet, library, 0)
  mSet <- getExportedValue("MetaboAnalystR", "CalculateHyperScore")(mSet)
  out <- if (!is.null(mSet$analSet$ora.mat)) as.data.frame(mSet$analSet$ora.mat) else data.frame()
  out$pathway <- rownames(out)
  .metab_standardize_metaboanalyst_result(out, p_adjust_method)
}

.metab_prepare_metaboanalyst_qs_files <- function(library) {
  .metab_prepare_qs_file(
    "compound_db.qs",
    "https://www.metaboanalyst.ca/resources/libs/compound_db.qs"
  )
  .metab_prepare_qs_file(
    "syn_nms.qs",
    "https://www.metaboanalyst.ca/resources/libs/syn_nms.qs"
  )
  .metab_prepare_qs_file(
    paste0(library, ".qs"),
    paste0("https://www.metaboanalyst.ca/resources/libs/msets/", library, ".qs")
  )
}

.metab_prepare_qs_file <- function(file, url) {
  if (!file.exists(file)) {
    utils::download.file(url, file, method = "libcurl", quiet = TRUE)
  }
  ok <- tryCatch({
    qs::qread(file)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop("MetaboAnalystR downloaded a qs2-format library file. Install qs2 and rerun.")
    }
    obj <- qs2::qs_read(file)
    qs::qsave(obj, file)
  }
  invisible(file)
}

.metab_standardize_metaboanalyst_result <- function(out, p_adjust_method) {
  if (!nrow(out)) {
    return(data.frame())
  }
  names(out) <- make.names(names(out))
  if (!"pathway" %in% names(out) && "Pathway" %in% names(out)) {
    out$pathway <- out$Pathway
  }
  p_col <- intersect(c("Raw.p", "Raw.p.", "p.value", "P.value"), names(out))[1]
  hits_col <- intersect(c("hits", "Hits"), names(out))[1]
  expected_col <- intersect(c("expected", "Expected"), names(out))[1]
  if (is.na(p_col)) {
    stop("Could not identify a raw p-value column in MetaboAnalystR ORA output.")
  }
  out$pvalue <- as.numeric(out[[p_col]])
  out$hits <- if (!is.na(hits_col)) as.numeric(out[[hits_col]]) else NA_real_
  out$expected <- if (!is.na(expected_col)) as.numeric(out[[expected_col]]) else NA_real_
  out$fold_enrichment <- ifelse(!is.na(out$expected) & out$expected > 0, out$hits / out$expected, NA_real_)
  out$p_adjust <- stats::p.adjust(out$pvalue, method = p_adjust_method)
  out$neg_log10_p <- -log10(pmax(out$pvalue, .Machine$double.xmin))
  out <- out[order(out$pvalue), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.metab_extract_ora_df <- function(x) {
  if (inherits(x, "ukb_metabolite_ora")) {
    return(x$ora_result)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  stop("`x` must be a ukb_metabolite_ora object or a data.frame.")
}

.metab_prepare_plot_df <- function(plot_df, top_n, p_col, pathway_col) {
  if (!nrow(plot_df)) {
    stop("No ORA rows available for plotting.")
  }
  required <- c(p_col, pathway_col, "fold_enrichment")
  if (!all(required %in% names(plot_df))) {
    stop("ORA result must contain columns: ", paste(required, collapse = ", "))
  }
  plot_df$p_for_plot <- as.numeric(plot_df[[p_col]])
  plot_df <- plot_df[is.finite(plot_df$p_for_plot), , drop = FALSE]
  plot_df <- plot_df[order(plot_df$p_for_plot), , drop = FALSE]
  plot_df <- utils::head(plot_df, top_n)
  plot_df$neg_log10_p <- -log10(pmax(plot_df$p_for_plot, .Machine$double.xmin))
  plot_df$pathway_label <- factor(plot_df[[pathway_col]], levels = rev(plot_df[[pathway_col]]))
  plot_df
}

.metab_theme <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = base_size * 0.95),
      legend.text = element_text(size = base_size * 0.85),
      panel.grid = element_blank()
    )
}

.metab_clean_input <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

.metab_normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}
