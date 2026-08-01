#' Combine Harmonized PRS Weights into One PLINK2 Score Matrix
#'
#' @description
#' Combine several harmonized PRS definitions into a wide score file. PLINK2
#' can then calculate all scores during one pass over each genotype file with
#' `--score-col-nums`, instead of rereading the same genotypes for every PRS.
#' Missing variant/score combinations are assigned weight zero.
#'
#' @param weights A named list of `ukb_prs_harmonized` objects.
#' @param score_names Optional output score names. Defaults to `names(weights)`.
#' @param output Optional tab-delimited output path.
#' @param overwrite Overwrite an existing output.
#'
#' @return A `ukb_prs_weight_matrix` data.frame with `ID`, `A1`, one column per
#'   score, and target variant metadata.
#' @export
combine_prs_weights <- function(weights,
                                score_names = names(weights),
                                output = NULL,
                                overwrite = FALSE) {
  .gpp_assert_flag(overwrite, "overwrite")
  if (!is.list(weights) || length(weights) < 2L ||
      !all(vapply(weights, inherits, logical(1), "ukb_prs_harmonized"))) {
    stop("`weights` must be a list of at least two harmonized PRS objects.", call. = FALSE)
  }
  score_names <- .prs_score_names(score_names, length(weights))
  builds <- vapply(weights, function(x) {
    build <- attr(x, "genome_build")
    if (is.null(build)) NA_character_ else build
  }, character(1))
  if (anyNA(builds) || length(unique(builds)) != 1L) {
    stop("All PRS weights must declare the same target genome build.", call. = FALSE)
  }

  required <- c("ID", "A1", "WEIGHT", "A2", "CHR", "POS")
  if (any(!vapply(weights, function(x) all(required %in% names(x)), logical(1)))) {
    stop("Each harmonized PRS must contain ID, A1, WEIGHT, A2, CHR, and POS.", call. = FALSE)
  }
  parts <- lapply(seq_along(weights), function(i) {
    x <- as.data.frame(weights[[i]], stringsAsFactors = FALSE)
    data.frame(
      KEY = paste(x$ID, x$A1, sep = "\r"),
      ID = x$ID, A1 = x$A1, A2 = x$A2, CHR = x$CHR, POS = x$POS,
      SCORE = score_names[[i]], WEIGHT = x$WEIGHT,
      stringsAsFactors = FALSE
    )
  })
  long <- data.table::rbindlist(parts, use.names = TRUE)
  metadata <- unique(long[, .(KEY, ID, A1, A2, CHR, POS)])
  conflict <- metadata[, .N, by = KEY][N > 1L, KEY]
  if (length(conflict) > 0L) {
    stop("PRS definitions contain conflicting target metadata for the same ID/A1.", call. = FALSE)
  }
  metadata <- metadata[!duplicated(metadata$KEY)]
  matrix <- matrix(0, nrow = nrow(metadata), ncol = length(score_names))
  colnames(matrix) <- score_names
  for (i in seq_along(score_names)) {
    part <- long[SCORE == score_names[[i]]]
    matrix[match(part$KEY, metadata$KEY), i] <- part$WEIGHT
  }
  out <- data.frame(
    ID = metadata$ID,
    A1 = metadata$A1,
    as.data.frame(matrix, check.names = FALSE),
    A2 = metadata$A2,
    CHR = metadata$CHR,
    POS = metadata$POS,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(output)) {
    output <- .gpt_path(output, "output")
    if (file.exists(output) && !isTRUE(overwrite)) {
      stop("Output already exists; use `overwrite = TRUE`.", call. = FALSE)
    }
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(out, output, sep = "\t", quote = FALSE, na = "NA")
    attr(out, "path") <- normalizePath(output, mustWork = FALSE)
  }
  attr(out, "score_names") <- score_names
  attr(out, "score_columns") <- seq.int(3L, 2L + length(score_names))
  attr(out, "variant_counts") <- stats::setNames(
    vapply(weights, nrow, integer(1)), score_names
  )
  attr(out, "genome_build") <- builds[[1L]]
  class(out) <- c("ukb_prs_weight_matrix", "data.frame")
  out
}

#' Split PRS Weights by Chromosome
#'
#' @description
#' Write one score file per chromosome. Empty chromosomes are skipped by
#' default, so [ukb_plan_prs()] does not launch jobs which cannot contribute to
#' the score.
#'
#' @param weights A `ukb_prs_harmonized` or `ukb_prs_weight_matrix` object.
#' @param output_dir Directory for chromosome-specific score files.
#' @param prefix File-name prefix.
#' @param chromosomes Optional chromosome subset. Defaults to chromosomes
#'   present in `weights`.
#' @param drop_empty Skip requested chromosomes without weights.
#' @param overwrite Overwrite existing files.
#'
#' @return A `ukb_prs_weight_files` object containing named file paths and row
#'   counts.
#' @export
split_prs_weights <- function(weights,
                              output_dir,
                              prefix = "prs_weights",
                              chromosomes = NULL,
                              drop_empty = TRUE,
                              overwrite = FALSE) {
  .gpp_assert_flag(drop_empty, "drop_empty")
  .gpp_assert_flag(overwrite, "overwrite")
  supported <- inherits(weights, "ukb_prs_harmonized") ||
    inherits(weights, "ukb_prs_weight_matrix")
  if (!supported || !"CHR" %in% names(weights)) {
    stop("`weights` must be harmonized weights with a CHR column.", call. = FALSE)
  }
  output_dir <- .gpt_path(output_dir, "output_dir")
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !grepl("^[A-Za-z0-9_.-]+$", prefix)) {
    stop("`prefix` must be one simple file-name prefix.", call. = FALSE)
  }
  available <- unique(.prs_normalize_chr(weights$CHR))
  requested <- if (is.null(chromosomes)) available else .prs_normalize_chr(chromosomes)
  requested <- unique(requested)
  counts <- vapply(requested, function(chr) sum(weights$CHR == chr), integer(1))
  if (isTRUE(drop_empty)) {
    requested <- requested[counts > 0L]
    counts <- counts[counts > 0L]
  }
  if (length(requested) == 0L) {
    stop("No chromosome-specific PRS weights remain after filtering.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  labels <- paste0("chr", requested)
  paths <- file.path(output_dir, paste0(prefix, "_", labels, ".tsv"))
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L && !isTRUE(overwrite)) {
    stop("Output already exists; use `overwrite = TRUE`: ",
         paste(existing, collapse = ", "), call. = FALSE)
  }
  for (i in seq_along(requested)) {
    part <- as.data.frame(weights[weights$CHR == requested[[i]], , drop = FALSE])
    data.table::fwrite(part, paths[[i]], sep = "\t", quote = FALSE, na = "NA")
  }
  files <- stats::setNames(normalizePath(paths, mustWork = FALSE), labels)
  out <- list(
    files = files,
    counts = stats::setNames(unname(counts), labels),
    score_names = if (inherits(weights, "ukb_prs_weight_matrix")) {
      attr(weights, "score_names")
    } else {
      NULL
    },
    score_columns = if (inherits(weights, "ukb_prs_weight_matrix")) {
      attr(weights, "score_columns")
    } else {
      3L
    },
    genome_build = attr(weights, "genome_build"),
    wide = inherits(weights, "ukb_prs_weight_matrix")
  )
  class(out) <- c("ukb_prs_weight_files", "list")
  out
}

#' @export
print.ukb_prs_weight_matrix <- function(x, ...) {
  cat("Combined PRS weight matrix\n")
  cat("  Variants: ", nrow(x), "\n", sep = "")
  cat("  Scores: ", paste(attr(x, "score_names"), collapse = ", "), "\n", sep = "")
  cat("  File: ", if (is.null(attr(x, "path"))) "not written" else attr(x, "path"), "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_prs_weight_files <- function(x, ...) {
  cat("Chromosome-specific PRS weight files\n")
  cat("  Files: ", length(x$files), "\n", sep = "")
  cat("  Variants: ", sum(x$counts), "\n", sep = "")
  invisible(x)
}

.prs_score_names <- function(x, n = NULL) {
  if (is.null(x) || !is.character(x) || length(x) == 0L || anyNA(x)) {
    stop("`score_names` must contain one or more names.", call. = FALSE)
  }
  if (!is.null(n) && length(x) != n) {
    stop("`score_names` must have one name per PRS definition.", call. = FALSE)
  }
  out <- vapply(x, .prs_score_name, character(1))
  if (anyDuplicated(out)) {
    stop("`score_names` must be unique.", call. = FALSE)
  }
  unname(out)
}
