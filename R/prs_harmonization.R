#' Harmonize PRS Weights to a Target Variant Map
#'
#' @description
#' Align prepared effect alleles to target genotype alleles before PLINK2
#' scoring. Matching can use variant ID, chromosome/position, or ID followed by
#' a position fallback. Direct, swapped, strand, and strand-swapped SNPs are
#' supported. Indels are retained only when their alleles match directly or in
#' swapped order; this function does not normalize or left-align indels.
#'
#' Input and target genome builds must be declared and must match. This
#' function deliberately does not perform liftOver. Palindromic SNPs are
#' dropped by default; frequency inference requires both input effect-allele
#' frequency and target allele-1 frequency.
#'
#' @param weights A `ukb_prs_weights` object from [prepare_prs_weights()].
#' @param variant_map A data.frame or delimited target variant-map path.
#' @param target_id_col,target_chr_col,target_pos_col Target variant columns.
#' @param target_allele1_col,target_allele2_col Target allele columns. For a
#'   REF/ALT map, allele 1 would commonly be ALT and allele 2 REF.
#' @param target_allele1_freq_col Optional frequency column for target allele 1.
#' @param target_build Target genome build, `"GRCh37"` or `"GRCh38"`.
#' @param match_by `"auto"` uses ID first and position as a fallback;
#'   `"id"` or `"position"` restrict matching to one key.
#' @param palindromic `"drop"`, `"infer_by_af"`, or `"keep"`. The last option
#'   trusts literal allele orientation and is not recommended without external
#'   justification.
#' @param af_tolerance Maximum absolute frequency difference for palindromic
#'   inference. Exactly one orientation must pass this threshold.
#' @param output Optional normalized PLINK2 weight-file path.
#' @param overwrite Overwrite an existing output.
#'
#' @return A `ukb_prs_harmonized` data.frame. Attributes
#'   `harmonization_qc` and `harmonization_exclusions` record decisions.
#' @export
harmonize_prs_weights <- function(
    weights,
    variant_map,
    target_id_col = "ID",
    target_chr_col = "CHR",
    target_pos_col = "POS",
    target_allele1_col = "ALT",
    target_allele2_col = "REF",
    target_allele1_freq_col = NULL,
    target_build,
    match_by = c("auto", "id", "position"),
    palindromic = c("drop", "infer_by_af", "keep"),
    af_tolerance = 0.1,
    output = NULL,
    overwrite = FALSE) {
  if (!inherits(weights, "ukb_prs_weights")) {
    stop("`weights` must come from `prepare_prs_weights()`.", call. = FALSE)
  }
  .gpp_assert_flag(overwrite, "overwrite")
  match_by <- match.arg(match_by)
  palindromic <- match.arg(palindromic)
  .gpt_positive_number(af_tolerance, "af_tolerance", allow_zero = TRUE)
  if (af_tolerance > 0.5) {
    stop("`af_tolerance` must be between 0 and 0.5.", call. = FALSE)
  }
  input_build <- .prs_genome_build(attr(weights, "genome_build"))
  target_build <- .prs_genome_build(target_build)
  if (!identical(input_build, target_build)) {
    stop(
      sprintf("Genome-build mismatch: weights use %s and target uses %s.", input_build, target_build),
      call. = FALSE
    )
  }

  target <- .prsh_read_table(variant_map, "variant_map")
  target_columns <- c(
    ID = .prs_resolve_column(target_id_col, target),
    CHR = .prs_resolve_column(target_chr_col, target),
    POS = .prs_resolve_column(target_pos_col, target),
    ALLELE1 = .prs_resolve_column(target_allele1_col, target),
    ALLELE2 = .prs_resolve_column(target_allele2_col, target)
  )
  target_freq_col <- .prs_optional_column(target_allele1_freq_col, target)
  target <- data.frame(
    TARGET_ID = trimws(as.character(target[[target_columns[["ID"]]]])),
    TARGET_CHR = .prs_normalize_chr(target[[target_columns[["CHR"]]]]),
    TARGET_POS = suppressWarnings(as.numeric(target[[target_columns[["POS"]]]])),
    TARGET_A1 = toupper(trimws(as.character(target[[target_columns[["ALLELE1"]]]]))),
    TARGET_A2 = toupper(trimws(as.character(target[[target_columns[["ALLELE2"]]]]))),
    TARGET_AF = if (is.null(target_freq_col)) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(target[[target_freq_col]]))
    },
    stringsAsFactors = FALSE
  )
  invalid_target <- !nzchar(target$TARGET_ID) | !is.finite(target$TARGET_POS) |
    target$TARGET_POS < 1 | target$TARGET_POS != floor(target$TARGET_POS) |
    !nzchar(target$TARGET_A1) | !nzchar(target$TARGET_A2) |
    target$TARGET_A1 == target$TARGET_A2
  if (any(invalid_target)) {
    stop(sprintf("Target variant map contains %d invalid row(s).", sum(invalid_target)), call. = FALSE)
  }
  if (!is.null(target_freq_col) &&
      any(!is.finite(target$TARGET_AF) | target$TARGET_AF < 0 | target$TARGET_AF > 1)) {
    stop("Target allele-1 frequencies must be in [0, 1].", call. = FALSE)
  }

  input <- as.data.frame(weights, stringsAsFactors = FALSE)
  if (identical(match_by, "position") && !all(c("CHR", "POS") %in% names(input))) {
    stop("Position matching requires `CHR` and `POS` in prepared weights.", call. = FALSE)
  }
  index <- .prsh_match_target(input, target, match_by)
  matched <- !is.na(index)
  target_rows <- target[index, , drop = FALSE]

  decision <- .prsh_align_alleles(
    input = input,
    target = target_rows,
    matched = matched,
    palindromic = palindromic,
    af_tolerance = af_tolerance,
    target_frequency_supplied = !is.null(target_freq_col)
  )
  keep <- decision$keep
  output_data <- data.frame(
    ID = target_rows$TARGET_ID[keep],
    A1 = decision$effect_allele[keep],
    WEIGHT = input$WEIGHT[keep],
    A2 = decision$other_allele[keep],
    CHR = target_rows$TARGET_CHR[keep],
    POS = target_rows$TARGET_POS[keep],
    EAF = decision$effect_allele_freq[keep],
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(paste(output_data$ID, output_data$A1, sep = "\r"))) {
    stop("Harmonization produced duplicated target ID/effect-allele rows.", call. = FALSE)
  }

  reason <- decision$reason
  alignment_counts <- as.data.frame(table(
    factor(decision$alignment[keep], levels = c(
      "direct", "swapped", "strand", "strand_swapped",
      "palindromic_direct", "palindromic_swapped"
    ))
  ), stringsAsFactors = FALSE)
  names(alignment_counts) <- c("alignment", "n")
  alignment_counts <- alignment_counts[alignment_counts$n > 0L, , drop = FALSE]
  exclusion_counts <- as.data.frame(table(reason[!keep]), stringsAsFactors = FALSE)
  names(exclusion_counts) <- c("reason", "n")
  qc <- list(
    summary = data.frame(
      metric = c("input", "target_matched", "retained", "excluded"),
      n = c(nrow(input), sum(matched), sum(keep), sum(!keep)),
      stringsAsFactors = FALSE
    ),
    alignment = alignment_counts,
    exclusions = exclusion_counts,
    settings = list(
      input_build = input_build,
      target_build = target_build,
      match_by = match_by,
      palindromic = palindromic,
      af_tolerance = af_tolerance
    )
  )
  exclusions <- data.frame(
    ID = input$ID[!keep],
    reason = reason[!keep],
    stringsAsFactors = FALSE
  )

  if (!is.null(output)) {
    output <- .gpt_path(output, "output")
    if (file.exists(output) && !isTRUE(overwrite)) {
      stop("Output already exists; use `overwrite = TRUE`.", call. = FALSE)
    }
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(output_data, output, sep = "\t", quote = FALSE, na = "NA")
    attr(output_data, "path") <- normalizePath(output, mustWork = FALSE)
  }
  attr(output_data, "source_columns") <- c(ID = "ID", A1 = "A1", WEIGHT = "WEIGHT")
  attr(output_data, "input_weight_type") <- attr(weights, "input_weight_type")
  attr(output_data, "weight_transformation") <- attr(weights, "weight_transformation")
  attr(output_data, "genome_build") <- target_build
  attr(output_data, "harmonization_qc") <- qc
  attr(output_data, "harmonization_exclusions") <- exclusions
  class(output_data) <- c("ukb_prs_harmonized", "ukb_prs_weights", "data.frame")
  output_data
}

#' @export
print.ukb_prs_harmonized <- function(x, ...) {
  qc <- attr(x, "harmonization_qc")
  retained <- nrow(x)
  input <- if (is.null(qc)) retained else qc$summary$n[qc$summary$metric == "input"]
  genome_build <- attr(x, "genome_build")
  cat("Harmonized PRS weights\n")
  cat("  Retained: ", retained, "/", input, "\n", sep = "")
  cat("  Genome build: ", if (is.null(genome_build)) "unknown" else genome_build, "\n", sep = "")
  invisible(x)
}

.prsh_read_table <- function(x, name) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop(sprintf("`%s` file does not exist: %s", name, x), call. = FALSE)
    }
    return(data.table::fread(x, data.table = FALSE, check.names = FALSE))
  }
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a data.frame or one existing file path.", name), call. = FALSE)
  }
  x
}

.prsh_match_target <- function(input, target, match_by) {
  index <- rep(NA_integer_, nrow(input))
  if (match_by %in% c("auto", "id")) {
    duplicated_target <- duplicated(target$TARGET_ID) | duplicated(target$TARGET_ID, fromLast = TRUE)
    unique_ids <- target$TARGET_ID
    unique_ids[duplicated_target] <- NA_character_
    index <- match(input$ID, unique_ids)
  }
  if (match_by %in% c("auto", "position")) {
    needs_position <- if (identical(match_by, "position")) rep(TRUE, nrow(input)) else is.na(index)
    if (any(needs_position) && all(c("CHR", "POS") %in% names(input))) {
      input_key <- paste(input$CHR, input$POS, sep = ":")
      target_key <- paste(target$TARGET_CHR, target$TARGET_POS, sep = ":")
      duplicated_target <- duplicated(target_key) | duplicated(target_key, fromLast = TRUE)
      target_key[duplicated_target] <- NA_character_
      index[needs_position] <- match(input_key[needs_position], target_key)
    }
  }
  index
}

.prsh_align_alleles <- function(input,
                                target,
                                matched,
                                palindromic,
                                af_tolerance,
                                target_frequency_supplied) {
  n <- nrow(input)
  w1 <- toupper(input$A1)
  has_a2 <- "A2" %in% names(input)
  w2 <- if (has_a2) toupper(input$A2) else rep(NA_character_, n)
  t1 <- target$TARGET_A1
  t2 <- target$TARGET_A2
  cw1 <- .prsh_complement(w1)
  cw2 <- .prsh_complement(w2)
  palindrome <- matched & .prsh_palindromic(t1, t2)

  direct <- matched & w1 == t1 & (!has_a2 | w2 == t2)
  swapped <- matched & w1 == t2 & (!has_a2 | w2 == t1)
  strand <- matched & !is.na(cw1) & cw1 == t1 & (!has_a2 | (!is.na(cw2) & cw2 == t2))
  strand_swapped <- matched & !is.na(cw1) & cw1 == t2 & (!has_a2 | (!is.na(cw2) & cw2 == t1))

  alignment <- rep(NA_character_, n)
  effect <- rep(NA_character_, n)
  other <- rep(NA_character_, n)
  effect_freq <- if ("EAF" %in% names(input)) input$EAF else rep(NA_real_, n)
  reason <- rep("allele_mismatch", n)
  reason[!matched] <- "target_not_matched"

  assign_rows <- function(rows, label, effect_value, other_value) {
    rows <- rows & is.na(alignment) & !palindrome
    alignment[rows] <<- label
    effect[rows] <<- effect_value[rows]
    other[rows] <<- other_value[rows]
  }
  assign_rows(direct, "direct", t1, t2)
  assign_rows(swapped, "swapped", t2, t1)
  assign_rows(strand, "strand", t1, t2)
  assign_rows(strand_swapped, "strand_swapped", t2, t1)

  if (any(palindrome)) {
    if (identical(palindromic, "drop")) {
      reason[palindrome] <- "palindromic_dropped"
    } else if (identical(palindromic, "keep")) {
      rows <- palindrome & direct
      alignment[rows] <- "palindromic_direct"
      effect[rows] <- t1[rows]
      other[rows] <- t2[rows]
      rows <- palindrome & is.na(alignment) & swapped
      alignment[rows] <- "palindromic_swapped"
      effect[rows] <- t2[rows]
      other[rows] <- t1[rows]
      reason[palindrome & is.na(alignment)] <- "palindromic_allele_mismatch"
    } else {
      has_input_af <- "EAF" %in% names(input)
      if (!has_input_af || !isTRUE(target_frequency_supplied)) {
        reason[palindrome] <- "palindromic_missing_af"
      } else {
        for (i in which(palindrome)) {
          candidates <- data.frame(
            alignment = c("palindromic_direct", "palindromic_swapped",
                          "palindromic_direct", "palindromic_swapped"),
            valid = c(direct[[i]], swapped[[i]], strand[[i]], strand_swapped[[i]]),
            effect = c(t1[[i]], t2[[i]], t1[[i]], t2[[i]]),
            other = c(t2[[i]], t1[[i]], t2[[i]], t1[[i]]),
            expected_af = c(target$TARGET_AF[[i]], 1 - target$TARGET_AF[[i]],
                            target$TARGET_AF[[i]], 1 - target$TARGET_AF[[i]]),
            stringsAsFactors = FALSE
          )
          candidates <- candidates[candidates$valid, , drop = FALSE]
          if (nrow(candidates) == 0L) {
            reason[[i]] <- "palindromic_allele_mismatch"
            next
          }
          candidates <- candidates[!duplicated(candidates[, c("effect", "expected_af")]), , drop = FALSE]
          passes <- abs(input$EAF[[i]] - candidates$expected_af) <= af_tolerance
          if (sum(passes) != 1L) {
            reason[[i]] <- if (sum(passes) == 0L) {
              "palindromic_af_mismatch"
            } else {
              "palindromic_af_ambiguous"
            }
            next
          }
          chosen <- candidates[which(passes), , drop = FALSE]
          alignment[[i]] <- chosen$alignment[[1L]]
          effect[[i]] <- chosen$effect[[1L]]
          other[[i]] <- chosen$other[[1L]]
          effect_freq[[i]] <- input$EAF[[i]]
        }
      }
    }
  }

  keep <- !is.na(alignment)
  reason[keep] <- "retained"
  if (isTRUE(target_frequency_supplied)) {
    effect_freq[keep & effect == t1] <- target$TARGET_AF[keep & effect == t1]
    effect_freq[keep & effect == t2] <- 1 - target$TARGET_AF[keep & effect == t2]
  }
  list(
    keep = keep,
    alignment = alignment,
    effect_allele = effect,
    other_allele = other,
    effect_allele_freq = effect_freq,
    reason = reason
  )
}

.prsh_complement <- function(x) {
  map <- c(A = "T", T = "A", C = "G", G = "C")
  out <- unname(map[x])
  out[is.na(x)] <- NA_character_
  out
}

.prsh_palindromic <- function(a1, a2) {
  pair <- paste0(pmin(a1, a2), pmax(a1, a2))
  pair %in% c("AT", "CG")
}
