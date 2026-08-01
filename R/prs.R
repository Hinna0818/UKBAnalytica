#' UK Biobank PRS Release V2 Catalogue
#'
#' @description
#' Return the participant-level polygenic risk score fields in UK Biobank
#' Showcase categories 301 (Standard) and 302 (Enhanced). Standard scores use
#' external GWAS training data. Enhanced scores also use within-UKB training
#' data and are available only in the PRS testing subgroup (field 26200).
#'
#' @param type Score set to return: `"standard"`, `"enhanced"`, or `"all"`.
#' @param trait Optional trait name or abbreviation, for example
#'   `"hypertension"`, `"HT"`, or `"T2D"`.
#'
#' @return A data.frame containing score type, trait, abbreviation, field ID,
#'   and testing-subgroup requirement.
#' @export
ukb_prs_catalog <- function(type = c("all", "standard", "enhanced"),
                            trait = NULL) {
  type <- match.arg(type)
  traits <- c(
    AAM = "age at menopause", AMD = "age-related macular degeneration",
    AD = "alzheimer's disease", APOEA = "apolipoprotein a1",
    APOEB = "apolipoprotein b", AST = "asthma", AF = "atrial fibrillation",
    BD = "bipolar disorder", BMI = "body mass index", CRC = "bowel cancer",
    BC = "breast cancer", CAL = "calcium", CVD = "cardiovascular disease",
    CED = "coeliac disease", CAD = "coronary artery disease",
    CD = "crohn's disease", DOA = "docosahexaenoic acid",
    EOC = "epithelial ovarian cancer",
    EBMDT = "estimated bone mineral density t-score",
    EGCR = "estimated glomerular filtration rate (creatinine based)",
    EGCY = "estimated glomerular filtration rate (cystatin based)",
    HBA1C_DF = "glycated haemoglobin", HEIGHT = "height",
    HDL = "high density lipoprotein cholesterol", HT = "hypertension",
    IOP = "intraocular pressure", ISS = "ischaemic stroke",
    LDL_SF = "low density lipoprotein cholesterol", MEL = "melanoma",
    MS = "multiple sclerosis", OTFA = "omega-3 fatty acids",
    OSFA = "omega-6 fatty acids", OP = "osteoporosis",
    PD = "parkinson's disease", PDCL = "phosphatidylcholines",
    PHG = "phosphoglycerides", PFA = "polyunsaturated fatty acids",
    POAG = "primary open angle glaucoma", PC = "prostate cancer",
    PSO = "psoriasis", RMNC = "remnant cholesterol", RHR = "resting heart rate",
    RA = "rheumatoid arthritis", SCZ = "schizophrenia",
    SGM = "sphingomyelins", SLE = "systemic lupus erythematosus",
    TCH = "total cholesterol", TFA = "total fatty acids",
    TTG = "total triglycerides", T1D = "type 1 diabetes",
    T2D = "type 2 diabetes", UC = "ulcerative colitis",
    VTE = "venous thromboembolic disease"
  )
  standard <- c(
    AAM = 26202L, AMD = 26204L, AD = 26206L, AST = 26210L, AF = 26212L,
    BD = 26214L, BMI = 26216L, CRC = 26218L, BC = 26220L, CVD = 26223L,
    CED = 26225L, CAD = 26227L, CD = 26229L, EOC = 26232L,
    EBMDT = 26234L, HBA1C_DF = 26238L, HEIGHT = 26240L, HDL = 26242L,
    HT = 26244L, IOP = 26246L, ISS = 26248L, LDL_SF = 26250L,
    MEL = 26252L, MS = 26254L, OP = 26258L, PD = 26260L,
    POAG = 26265L, PC = 26267L, PSO = 26269L, RHR = 21150L,
    RA = 26273L, SCZ = 26275L, SLE = 26278L, TCH = 21151L,
    TTG = 21152L, T1D = 26283L, T2D = 26285L, UC = 26287L,
    VTE = 26289L
  )
  enhanced <- c(
    AAM = 26203L, AMD = 26205L, AD = 26207L, APOEA = 26208L,
    APOEB = 26209L, AST = 26211L, AF = 26213L, BD = 26215L,
    BMI = 26217L, CRC = 26219L, BC = 26221L, CAL = 26222L,
    CVD = 26224L, CED = 26226L, CAD = 26228L, DOA = 26231L,
    EOC = 26233L, EBMDT = 26235L, EGCR = 26236L, EGCY = 26237L,
    HBA1C_DF = 26239L, HEIGHT = 26241L, HDL = 26243L, HT = 26245L,
    IOP = 26247L, ISS = 26249L, LDL_SF = 26251L, MEL = 26253L,
    MS = 26255L, OTFA = 26256L, OSFA = 26257L, OP = 26259L,
    PD = 26261L, PDCL = 26262L, PHG = 26263L, PFA = 26264L,
    POAG = 26266L, PC = 26268L, PSO = 26270L, RMNC = 26271L,
    RHR = 26272L, RA = 26274L, SCZ = 26276L, SGM = 26277L,
    SLE = 26279L, TCH = 26280L, TFA = 26281L, TTG = 26282L,
    T1D = 26284L, T2D = 26286L, VTE = 26290L
  )
  make_rows <- function(fields, score_type) {
    code <- names(fields)
    data.frame(
      type = score_type,
      trait = unname(traits[code]),
      code = code,
      field_id = unname(as.integer(fields)),
      testing_subgroup_required = identical(score_type, "enhanced"),
      release = "V2 (May 2024)",
      stringsAsFactors = FALSE
    )
  }
  out <- rbind(make_rows(standard, "standard"), make_rows(enhanced, "enhanced"))
  if (!identical(type, "all")) {
    out <- out[out$type == type, , drop = FALSE]
  }
  if (!is.null(trait)) {
    if (!is.character(trait) || length(trait) == 0L || anyNA(trait)) {
      stop("`trait` must be a non-empty character vector.", call. = FALSE)
    }
    query <- tolower(trimws(trait))
    keep <- tolower(out$code) %in% query | tolower(out$trait) %in% query
    missing <- query[!query %in% c(tolower(out$code), tolower(out$trait))]
    if (length(missing) > 0L) {
      stop(sprintf("Unknown PRS trait(s): %s.", paste(unique(missing), collapse = ", ")), call. = FALSE)
    }
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Load Published UK Biobank PRS Fields
#'
#' @description
#' Select already-dispensated PRS fields from a participant table, or extract
#' them from the current RAP dataset through [rap_extract_pheno()].
#' Enhanced requests also include field 26200 by default. No genotype data are
#' read and no participant-level data are stored in the package.
#'
#' @param data Optional participant data.frame. If `NULL`, fields are
#'   extracted from RAP.
#' @param trait Trait names or abbreviations accepted by
#'   [ukb_prs_catalog()].
#' @param field_id Optional explicit published PRS field IDs.
#' @param type `"standard"`, `"enhanced"`, or `"all"`.
#' @param id_col Participant ID column in `data`.
#' @param dataset Optional RAP dataset passed to [rap_extract_pheno()].
#' @param include_testing Include field 26200 for Enhanced PRS requests.
#' @param eligibility Enhanced PRS eligibility handling. `"set_na"` keeps all
#'   participants but clears Enhanced scores outside field 26200;
#'   `"require"` keeps testing-subgroup participants only; `"annotate"` only
#'   appends the indicator.
#' @param standardize Add within-sample z-scores with suffix `_z`.
#' @param dry_run Return a RAP extraction plan without running it.
#' @param ... Additional arguments passed to [rap_extract_pheno()].
#'
#' @return A data.table of participant IDs and named PRS columns, or a
#'   `rap_extract_plan` in dry-run mode.
#' @export
load_ukb_prs <- function(data = NULL,
                         trait = NULL,
                         field_id = NULL,
                         type = c("standard", "enhanced", "all"),
                         id_col = "eid",
                         dataset = NULL,
                         include_testing = TRUE,
                         eligibility = c("set_na", "require", "annotate"),
                         standardize = FALSE,
                         dry_run = FALSE,
                         ...) {
  type <- match.arg(type)
  eligibility <- match.arg(eligibility)
  .gpp_assert_flag(include_testing, "include_testing")
  .gpp_assert_flag(standardize, "standardize")
  .gpp_assert_flag(dry_run, "dry_run")
  if (is.null(trait) && is.null(field_id)) {
    stop("Provide `trait` or `field_id`.", call. = FALSE)
  }

  catalog <- ukb_prs_catalog("all")
  selected <- catalog[FALSE, , drop = FALSE]
  if (!is.null(trait)) {
    selected <- ukb_prs_catalog(type, trait)
  }
  if (!is.null(field_id)) {
    ids <- suppressWarnings(as.integer(field_id))
    if (length(ids) == 0L || anyNA(ids) || any(ids <= 0L)) {
      stop("`field_id` must contain positive UKB field IDs.", call. = FALSE)
    }
    known <- catalog[catalog$field_id %in% ids, , drop = FALSE]
    unknown <- setdiff(ids, known$field_id)
    if (length(unknown) > 0L) {
      stop(
        sprintf("Field ID(s) are not published PRS fields: %s.", paste(unknown, collapse = ", ")),
        call. = FALSE
      )
    }
    selected <- rbind(selected, known)
  }
  selected <- selected[!duplicated(selected$field_id), , drop = FALSE]
  request_ids <- selected$field_id
  needs_testing <- any(selected$testing_subgroup_required)
  if (needs_testing && !isTRUE(include_testing) && !identical(eligibility, "annotate")) {
    stop("Enhanced PRS eligibility requires `include_testing = TRUE`.", call. = FALSE)
  }
  if (isTRUE(include_testing) && needs_testing) {
    request_ids <- c(request_ids, 26200L)
  }

  if (is.null(data)) {
    extracted <- rap_extract_pheno(
      field_id = request_ids,
      dataset = dataset,
      strip_entity_prefix = TRUE,
      dry_run = dry_run,
      ...
    )
    if (isTRUE(dry_run)) {
      return(extracted)
    }
    data <- extracted
  } else if (isTRUE(dry_run)) {
    stop("`dry_run` is only used when `data = NULL`.", call. = FALSE)
  }

  data <- data.table::as.data.table(data.table::copy(data))
  id_name <- .prs_find_id(names(data), id_col)
  out <- data[, id_name, with = FALSE]
  data.table::setnames(out, id_name, id_col)
  score_columns <- character()
  enhanced_columns <- character()
  for (i in seq_len(nrow(selected))) {
    source <- .prs_find_field(names(data), selected$field_id[[i]])
    target <- paste0("prs_", tolower(selected$code[[i]]), "_", selected$type[[i]])
    out[[target]] <- data[[source]]
    score_columns <- c(score_columns, target)
    if (identical(selected$type[[i]], "enhanced")) {
      enhanced_columns <- c(enhanced_columns, target)
    }
  }
  eligibility_summary <- NULL
  if (isTRUE(include_testing) && needs_testing) {
    testing <- .prs_find_field(names(data), 26200L, required = FALSE)
    if (!is.null(testing)) {
      out$prs_testing_subgroup <- data[[testing]]
      eligible <- .prs_testing_eligible(out$prs_testing_subgroup)
      nonmissing_outside <- sum(
        !eligible & rowSums(!is.na(as.data.frame(out[, enhanced_columns, with = FALSE]))) > 0L
      )
      eligibility_summary <- data.frame(
        policy = eligibility,
        n_input = nrow(out),
        n_testing_eligible = sum(eligible),
        n_not_eligible = sum(!eligible),
        n_nonmissing_outside_testing = nonmissing_outside,
        stringsAsFactors = FALSE
      )
      if (identical(eligibility, "set_na")) {
        for (column in enhanced_columns) {
          data.table::set(out, which(!eligible), column, NA_real_)
        }
      } else if (identical(eligibility, "require")) {
        out <- out[eligible]
      }
    } else {
      if (!identical(eligibility, "annotate")) {
        stop("Enhanced PRS eligibility field 26200 was not found.", call. = FALSE)
      }
      warning("Enhanced PRS loaded without testing-subgroup field 26200.", call. = FALSE)
    }
  }
  if (isTRUE(standardize)) {
    for (column in score_columns) {
      out[[paste0(column, "_z")]] <- .prs_zscore(out[[column]])
    }
  }
  attr(out, "ukb_prs_fields") <- selected
  attr(out, "ukb_prs_eligibility") <- eligibility_summary
  class(out) <- unique(c("ukb_prs_data", class(out)))
  out
}

#' Prepare a PLINK2 PRS Weight File
#'
#' @param weights A data.frame or a delimited weight-file path.
#' @param variant_col,effect_allele_col,weight_col Column names or positions.
#' @param other_allele_col,chr_col,pos_col,effect_allele_freq_col Optional
#'   columns used by [harmonize_prs_weights()].
#' @param weight_type Input scale: `"beta"`, `"log_or"`, and `"direct"` are
#'   used unchanged; `"or"` is transformed with `log()`.
#' @param genome_build Optional `"GRCh37"` or `"GRCh38"` metadata.
#' @param output Optional normalized tab-delimited output path. An explicit
#'   output is required before the object can be passed to
#'   [ukb_plan_prs()].
#' @param overwrite Overwrite an existing output.
#'
#' @return A normalized data.frame of class `ukb_prs_weights`, with columns
#'   `ID`, `A1`, and `WEIGHT`.
#' @export
prepare_prs_weights <- function(weights,
                                variant_col = "ID",
                                effect_allele_col = "A1",
                                weight_col = "WEIGHT",
                                other_allele_col = NULL,
                                chr_col = NULL,
                                pos_col = NULL,
                                effect_allele_freq_col = NULL,
                                weight_type = c("beta", "log_or", "or", "direct"),
                                genome_build = NULL,
                                output = NULL,
                                overwrite = FALSE) {
  .gpp_assert_flag(overwrite, "overwrite")
  weight_type <- match.arg(weight_type)
  genome_build <- .prs_genome_build(genome_build, allow_null = TRUE)
  source <- if (is.character(weights) && length(weights) == 1L) {
    if (!file.exists(weights)) {
      stop("PRS weight file does not exist: ", weights, call. = FALSE)
    }
    data.table::fread(weights, data.table = FALSE, check.names = FALSE)
  } else if (is.data.frame(weights)) {
    weights
  } else {
    stop("`weights` must be a data.frame or one existing file path.", call. = FALSE)
  }
  required_columns <- vapply(
    list(variant_col, effect_allele_col, weight_col),
    .prs_resolve_column,
    character(1),
    data = source
  )
  optional_specs <- list(
    A2 = other_allele_col, CHR = chr_col, POS = pos_col,
    EAF = effect_allele_freq_col
  )
  optional_columns <- lapply(optional_specs, .prs_optional_column, data = source)
  raw_weight <- suppressWarnings(as.numeric(source[[required_columns[[3L]]]]))
  if (identical(weight_type, "or")) {
    if (any(!is.finite(raw_weight) | raw_weight <= 0)) {
      stop("OR weights must be finite and greater than zero.", call. = FALSE)
    }
    raw_weight <- log(raw_weight)
  }
  out <- data.frame(
    ID = trimws(as.character(source[[required_columns[[1L]]]])),
    A1 = toupper(trimws(as.character(source[[required_columns[[2L]]]]))),
    WEIGHT = raw_weight,
    stringsAsFactors = FALSE
  )
  invalid <- !nzchar(out$ID) | !nzchar(out$A1) | !is.finite(out$WEIGHT)
  if (any(invalid)) {
    stop(sprintf("PRS weights contain %d invalid row(s).", sum(invalid)), call. = FALSE)
  }
  if (anyDuplicated(paste(out$ID, out$A1, sep = "\r"))) {
    stop("PRS weights contain duplicated ID/effect-allele rows.", call. = FALSE)
  }
  if (!is.null(optional_columns$A2)) {
    out$A2 <- toupper(trimws(as.character(source[[optional_columns$A2]])))
    if (any(!nzchar(out$A2)) || any(out$A1 == out$A2)) {
      stop("Other alleles must be non-empty and differ from effect alleles.", call. = FALSE)
    }
  }
  if (!is.null(optional_columns$CHR)) {
    out$CHR <- .prs_normalize_chr(source[[optional_columns$CHR]])
  }
  if (!is.null(optional_columns$POS)) {
    out$POS <- suppressWarnings(as.numeric(source[[optional_columns$POS]]))
    if (any(!is.finite(out$POS) | out$POS < 1 | out$POS != floor(out$POS))) {
      stop("Variant positions must be positive integers.", call. = FALSE)
    }
  }
  if (!is.null(optional_columns$EAF)) {
    out$EAF <- suppressWarnings(as.numeric(source[[optional_columns$EAF]]))
    if (any(!is.finite(out$EAF) | out$EAF < 0 | out$EAF > 1)) {
      stop("Effect-allele frequencies must be in [0, 1].", call. = FALSE)
    }
  }
  if (xor("CHR" %in% names(out), "POS" %in% names(out))) {
    stop("Provide both chromosome and position columns, or neither.", call. = FALSE)
  }
  if (!is.null(output)) {
    output <- .gpt_path(output, "output")
    if (file.exists(output) && !isTRUE(overwrite)) {
      stop("Output already exists; use `overwrite = TRUE`.", call. = FALSE)
    }
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(out, output, sep = "\t", quote = FALSE, na = "NA")
    attr(out, "path") <- normalizePath(output, mustWork = FALSE)
  }
  attr(out, "source_columns") <- c(
    ID = required_columns[[1L]], A1 = required_columns[[2L]],
    WEIGHT = required_columns[[3L]], unlist(optional_columns, use.names = TRUE)
  )
  attr(out, "input_weight_type") <- weight_type
  attr(out, "weight_transformation") <- if (identical(weight_type, "or")) "log" else "none"
  attr(out, "genome_build") <- genome_build
  class(out) <- c("ukb_prs_weights", "data.frame")
  out
}

#' Plan PRS Calculation with PLINK2
#'
#' @description
#' Create one auditable PLINK2 `--score` command per genotype input. The plan
#' requests score sums so chromosome-level results can be added without
#' averaging artefacts. It does not execute PLINK2 or submit a DNAnexus job.
#'
#' @param weights Normalized output from [prepare_prs_weights()], a combined
#'   matrix from [combine_prs_weights()], chromosome files from
#'   [split_prs_weights()], or a PLINK2 score-file path.
#' @param bgen One or more BGEN paths. Retained for backward compatibility;
#'   new code can use `genotype` with `genotype_format`.
#' @param sample One sample path or one path per BGEN. Must be `NULL` for
#'   PGEN and BED inputs, whose sample information is part of the file set.
#' @param output_prefix Output prefix.
#' @param keep_file Optional PLINK two-column sample inclusion file.
#' @param variant_col,effect_allele_col,weight_col One-based columns in a raw
#'   score file. Normalized weights always use columns 1, 2, and 3.
#' @param header Whether the score file has a header.
#' @param ref_first Add the PLINK2 BGEN `ref-first` modifier.
#' @param no_mean_imputation Add the PLINK2 modifier of the same name.
#' @param frequency_file Optional PLINK2 allele-frequency file passed through
#'   `--read-freq` to make mean imputation independent of the scoring cohort.
#' @param error_on_freq_calc Add `--error-on-freq-calc`. By default this is
#'   enabled when `frequency_file` is supplied.
#' @param list_variants Write the variants actually used by PLINK2.
#' @param score_name Name used for a single-score R output. Combined weight
#'   matrices retain their own score names.
#' @param plink2_args Additional non-conflicting PLINK2 argument tokens.
#' @param genotype One or more BGEN paths, PGEN prefixes/paths, or BED
#'   prefixes/paths.
#' @param genotype_format Input format: `"auto"`, `"bgen"`, `"pgen"`, or
#'   `"bed"`.
#'
#' @return A `ukb_prs_plan` inheriting from `ukb_plink2_plan`.
#' @export
ukb_plan_prs <- function(weights,
                         bgen = NULL,
                         sample = NULL,
                         output_prefix,
                         keep_file = NULL,
                         variant_col = 1L,
                         effect_allele_col = 2L,
                         weight_col = 3L,
                         header = TRUE,
                         ref_first = TRUE,
                         no_mean_imputation = FALSE,
                         frequency_file = NULL,
                         error_on_freq_calc = !is.null(frequency_file),
                         list_variants = TRUE,
                         score_name = "PRS",
                         plink2_args = character(),
                         genotype = NULL,
                         genotype_format = c("auto", "bgen", "pgen", "bed")) {
  .gpp_assert_flag(header, "header")
  .gpp_assert_flag(ref_first, "ref_first")
  .gpp_assert_flag(no_mean_imputation, "no_mean_imputation")
  .gpp_assert_flag(error_on_freq_calc, "error_on_freq_calc")
  .gpp_assert_flag(list_variants, "list_variants")
  genotype_format <- match.arg(genotype_format)
  wide <- FALSE
  split <- inherits(weights, "ukb_prs_weight_files")
  if (split) {
    score_files <- weights$files
    wide <- isTRUE(weights$wide)
    score_columns <- as.integer(weights$score_columns)
    score_names <- if (wide) weights$score_names else .prs_score_names(score_name, 1L)
    header <- TRUE
  } else if (inherits(weights, "ukb_prs_weight_matrix")) {
    score_files <- attr(weights, "path")
    if (is.null(score_files)) {
      stop("Write the combined weights with `output=` before planning PRS.", call. = FALSE)
    }
    wide <- TRUE
    score_columns <- as.integer(attr(weights, "score_columns"))
    score_names <- attr(weights, "score_names")
    header <- TRUE
  } else if (inherits(weights, "ukb_prs_weights")) {
    score_file <- attr(weights, "path")
    if (is.null(score_file)) {
      stop("Write normalized weights with `output=` before planning PRS.", call. = FALSE)
    }
    variant_col <- 1L
    effect_allele_col <- 2L
    weight_col <- 3L
    score_files <- score_file
    score_columns <- 3L
    score_names <- .prs_score_names(score_name, 1L)
    header <- TRUE
  } else {
    score_files <- .gpt_path(weights, "weights")
    score_columns <- vapply(
      list(variant_col, effect_allele_col, weight_col),
      .prs_positive_integer,
      integer(1)
    )
    if (length(unique(score_columns)) != 3L) {
      stop("PRS score columns must be distinct.", call. = FALSE)
    }
    score_names <- .prs_score_names(score_name, 1L)
  }
  if (!wide && !split && inherits(weights, "ukb_prs_weights")) {
    score_columns <- c(1L, 2L, 3L)
  }
  if (!is.null(genotype) && !is.null(bgen)) {
    stop("Provide only one of `genotype` and `bgen`.", call. = FALSE)
  }
  legacy_bgen <- is.null(genotype)
  genotype <- if (legacy_bgen) bgen else genotype
  if (is.null(genotype)) {
    stop("Provide `genotype` or `bgen`.", call. = FALSE)
  }
  if (legacy_bgen && identical(genotype_format, "auto")) {
    genotype_format <- "bgen"
  }
  input <- .prs_genotype_spec(genotype, genotype_format, sample, ref_first)
  output_prefix <- .gpt_path(output_prefix, "output_prefix")
  if (!is.null(keep_file)) {
    keep_file <- .gpt_path(keep_file, "keep_file")
  }
  if (!is.null(frequency_file)) {
    frequency_file <- .gpt_path(frequency_file, "frequency_file")
  }
  if (isTRUE(error_on_freq_calc) && is.null(frequency_file) && !isTRUE(no_mean_imputation)) {
    stop(
      "`error_on_freq_calc = TRUE` requires `frequency_file` unless mean imputation is disabled.",
      call. = FALSE
    )
  }
  plink2_args <- .gpt_cli_args(plink2_args, "plink2_args", allow_empty = TRUE)
  .gpt_reject_managed_flags(
    plink2_args,
    c("--bgen", "--sample", "--keep", "--score", "--score-list",
      "--score-col-nums", "--read-freq", "--error-on-freq-calc",
      "--out", "--bfile", "--pfile", "--vcf", "--bcf"),
    "ukb_plan_prs()"
  )

  labels <- input$labels
  if (split) {
    keep <- labels %in% names(score_files)
    if (!any(keep)) {
      stop("No genotype input has a chromosome label matching split PRS weights.", call. = FALSE)
    }
    labels <- labels[keep]
    input$args <- input$args[keep]
    input$validation <- input$validation[keep]
    input$paths <- input$paths[keep]
    if (!is.null(input$sample)) input$sample <- input$sample[keep]
    score_files <- unname(score_files[labels])
  } else {
    score_files <- rep(score_files, length(labels))
  }
  outputs <- if (length(labels) == 1L) output_prefix else paste0(output_prefix, "_", labels)
  commands <- lapply(seq_along(labels), function(i) {
    args <- input$args[[i]]
    if (!is.null(keep_file)) {
      args <- c(args, "--keep", keep_file)
    }
    if (!is.null(frequency_file)) {
      args <- c(args, "--read-freq", frequency_file)
    }
    if (isTRUE(error_on_freq_calc)) {
      args <- c(args, "--error-on-freq-calc")
    }
    score_args <- if (wide) {
      c("--score", score_files[[i]], "1", "2", "header-read",
        "--score-col-nums", .prs_column_range(score_columns))
    } else {
      columns <- if (split) c(1L, 2L, score_columns[[1L]]) else score_columns
      c("--score", score_files[[i]], as.character(columns),
        if (isTRUE(header)) "header-read" else character())
    }
    score_args <- c(
      score_args,
      if (isTRUE(no_mean_imputation)) "no-mean-imputation" else character(),
      if (isTRUE(list_variants)) "list-variants" else character(),
      "cols=maybefid,nallele,denom,scoresums"
    )
    .gpt_command("plink2", c(args, score_args, plink2_args, "--out", outputs[[i]]))
  })
  names(commands) <- labels
  out <- list(
    tool = "PLINK2",
    commands = commands,
    expected_sscore = stats::setNames(paste0(outputs, ".sscore"), labels),
    expected_variants = if (isTRUE(list_variants)) {
      stats::setNames(paste0(outputs, ".sscore.vars"), labels)
    } else {
      character()
    },
    settings = list(
      genotype = input$paths, genotype_format = input$format,
      genotype_validation = input$validation,
      weights = score_files, bgen = if (identical(input$format, "bgen")) input$paths else NULL,
      sample = input$sample,
      keep_file = keep_file, score_columns = score_columns, header = header,
      ref_first = ref_first, no_mean_imputation = no_mean_imputation,
      frequency_file = frequency_file, error_on_freq_calc = error_on_freq_calc,
      list_variants = list_variants,
      score_name = if (length(score_names) == 1L) score_names[[1L]] else score_names,
      score_names = score_names, wide = wide,
      plink2_args = plink2_args
    )
  )
  class(out) <- c("ukb_prs_plan", "ukb_plink2_plan", "list")
  out
}

#' Read and Combine PLINK2 PRS Scores
#'
#' @param files One or more PLINK2 `.sscore` files.
#' @param score_col Optional score-sum column name(s). By default, all columns
#'   ending in `_SUM` are used.
#' @param score_name Output score name(s). When omitted for a multi-score file,
#'   names are derived from the `_SUM` column headers.
#' @param standardize Add a z-score column.
#'
#' @return A data.table with `FID`, `IID`, score sums, allele denominator, and
#'   number of combined files. `PRS_AVG` is retained for a single score.
#' @export
read_prs_scores <- function(files,
                            score_col = NULL,
                            score_name = "PRS",
                            standardize = FALSE) {
  score_name_missing <- missing(score_name)
  files <- .gpt_paths(files, "files")
  .gpp_assert_flag(standardize, "standardize")
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("PRS score file(s) do not exist: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  input_score_columns <- NULL
  parts <- lapply(seq_along(files), function(i) {
    x <- data.table::fread(files[[i]], data.table = TRUE, check.names = FALSE)
    iid <- intersect(c("IID", "#IID"), names(x))
    if (length(iid) != 1L) {
      stop("Each .sscore file must contain one IID column.", call. = FALSE)
    }
    candidates <- if (is.null(score_col)) grep("_SUM$", names(x), value = TRUE) else score_col
    if (length(candidates) == 0L || any(!candidates %in% names(x))) {
      stop("No valid score-sum columns were found; use `score_col` to select them.", call. = FALSE)
    }
    if (is.null(input_score_columns)) {
      input_score_columns <<- candidates
    } else if (!identical(candidates, input_score_columns)) {
      stop("All .sscore files must contain the same score-sum columns in the same order.", call. = FALSE)
    }
    fid <- intersect(c("FID", "#FID"), names(x))
    out <- data.table::data.table(
      FID = if (length(fid) == 1L) as.character(x[[fid]]) else as.character(x[[iid]]),
      IID = as.character(x[[iid]]),
      allele_ct = if ("ALLELE_CT" %in% names(x)) as.numeric(x$ALLELE_CT) else NA_real_,
      denom = if ("DENOM" %in% names(x)) as.numeric(x$DENOM) else NA_real_,
      source_file = i
    )
    for (j in seq_along(candidates)) {
      out[[paste0("score_", j)]] <- suppressWarnings(as.numeric(x[[candidates[[j]]]]))
    }
    if (anyDuplicated(out[, .(FID, IID)])) {
      stop("A .sscore file contains duplicated FID/IID rows.", call. = FALSE)
    }
    out
  })
  if (score_name_missing) {
    score_names <- if (length(input_score_columns) == 1L) {
      "PRS"
    } else {
      sub("_SUM$", "", input_score_columns)
    }
  } else {
    score_names <- score_name
  }
  score_names <- .prs_score_names(score_names, length(input_score_columns))
  internal <- paste0("score_", seq_along(input_score_columns))
  long <- data.table::rbindlist(parts, use.names = TRUE)
  score_values <- long[, lapply(.SD, .prs_sum_or_na),
                       by = .(FID, IID), .SDcols = internal]
  metrics <- long[, .(
    ALLELE_CT = .prs_sum_or_na(allele_ct),
    DENOM = .prs_sum_or_na(denom),
    N_FILES = data.table::uniqueN(source_file)
  ), by = .(FID, IID)]
  out <- merge(score_values, metrics, by = c("FID", "IID"), sort = FALSE)
  if (any(out$N_FILES != length(files))) {
    warning("Some participants are absent from one or more .sscore files.", call. = FALSE)
  }
  data.table::setnames(out, internal, score_names)
  if (length(score_names) == 1L) {
    out$PRS_AVG <- out[[score_names]] / out$DENOM
  }
  if (isTRUE(standardize)) {
    for (name in score_names) {
      out[[paste0(name, "_Z")]] <- .prs_zscore(out[[name]])
    }
  }
  attr(out, "source_files") <- files
  attr(out, "source_score_columns") <- input_score_columns
  attr(out, "score_names") <- score_names
  class(out) <- unique(c("ukb_prs_scores", class(out)))
  out[]
}

#' Run a Planned PRS Calculation
#'
#' @param plan A `ukb_prs_plan`.
#' @param execute Execute PLINK2. Defaults to `FALSE`.
#' @param executable PLINK2 command or path.
#' @param require_rap Require a RAP-like environment before execution.
#' @param timeout Per-command timeout in seconds.
#' @param read_scores Read and combine successful `.sscore` outputs.
#' @param standardize Add a combined PRS z-score.
#' @param overwrite Permit PLINK2 to replace existing expected outputs.
#' @param workers Number of chromosome commands to execute concurrently.
#' @param resume Reuse complete existing chromosome outputs and execute only
#'   unfinished commands.
#'
#' @return The unchanged plan in dry-run mode; otherwise a
#'   `ukb_prs_run` containing PLINK2 logs and combined scores.
#' @export
ukb_run_prs <- function(plan,
                        execute = FALSE,
                        executable = "plink2",
                        require_rap = TRUE,
                        timeout = 86400,
                        read_scores = TRUE,
                        standardize = FALSE,
                        overwrite = FALSE,
                        workers = 1L,
                        resume = FALSE) {
  if (!inherits(plan, "ukb_prs_plan")) {
    stop("`plan` must be a `ukb_prs_plan`.", call. = FALSE)
  }
  .gpp_assert_flag(read_scores, "read_scores")
  .gpp_assert_flag(standardize, "standardize")
  .gpp_assert_flag(overwrite, "overwrite")
  .gpp_assert_flag(resume, "resume")
  .gpt_positive_integer(workers, "workers")
  if (!isTRUE(execute)) {
    return(plan)
  }
  inputs <- unique(c(
    plan$settings$weights, unlist(plan$settings$genotype_validation, use.names = FALSE),
    plan$settings$sample,
    plan$settings$keep_file, plan$settings$frequency_file
  ))
  missing <- inputs[!file.exists(inputs)]
  if (length(missing) > 0L) {
    stop("PRS input file(s) do not exist: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  labels <- names(plan$commands)
  expected_by_command <- stats::setNames(lapply(labels, function(label) {
    variants <- if (label %in% names(plan$expected_variants)) {
      plan$expected_variants[[label]]
    } else {
      character()
    }
    c(plan$expected_sscore[[label]], variants)
  }), labels)
  complete <- vapply(expected_by_command, function(paths) {
    all(file.exists(paths) & file.info(paths)$size > 0L)
  }, logical(1))
  expected <- unlist(expected_by_command, use.names = FALSE)
  existing <- expected[file.exists(expected)]
  if (length(existing) > 0L && !isTRUE(overwrite) && !isTRUE(resume)) {
    stop(
      "Expected PRS output file(s) already exist; use `overwrite = TRUE`: ",
      paste(existing, collapse = ", "),
      call. = FALSE
    )
  }
  partial <- vapply(expected_by_command, function(paths) {
    any(file.exists(paths)) && !all(file.exists(paths) & file.info(paths)$size > 0L)
  }, logical(1))
  if (isTRUE(resume) && any(partial) && !isTRUE(overwrite)) {
    stop(
      "Incomplete PRS output(s) exist; use `overwrite = TRUE` to rerun: ",
      paste(names(partial)[partial], collapse = ", "),
      call. = FALSE
    )
  }
  output_dirs <- unique(dirname(expected))
  if (any(!dir.exists(output_dirs))) {
    stop("The parent directory of a PRS output does not exist.", call. = FALSE)
  }
  if (isTRUE(require_rap)) {
    .ukb_assert_rap_env("ukb_run_prs()")
  }
  executable_path <- unname(Sys.which(executable))
  if (!nzchar(executable_path)) {
    stop("PLINK2 executable was not found on PATH.", call. = FALSE)
  }
  .gpt_positive_number(timeout, "timeout")
  pending <- if (isTRUE(resume)) labels[!complete] else labels
  run_one <- function(label) {
    command <- plan$commands[[label]]
    log <- system2(
      executable_path,
      args = vapply(command$args, shQuote, character(1)),
      stdout = TRUE,
      stderr = TRUE,
      timeout = timeout
    )
    status <- attr(log, "status")
    if (is.null(status)) status <- 0L
    list(
      label = label, status = as.integer(status),
      success = identical(as.integer(status), 0L), skipped = FALSE,
      command = .gpt_display_command(executable, command$args), log = log
    )
  }
  if (length(pending) > 0L && workers > 1L && .Platform$OS.type != "windows") {
    executed <- parallel::mclapply(
      pending, run_one, mc.cores = min(as.integer(workers), length(pending)),
      mc.preschedule = FALSE
    )
  } else {
    if (length(pending) > 0L && workers > 1L && .Platform$OS.type == "windows") {
      warning("Parallel PRS execution is unavailable on Windows; using one worker.", call. = FALSE)
    }
    executed <- lapply(pending, run_one)
  }
  names(executed) <- pending
  skipped <- labels[isTRUE(resume) & complete]
  skipped_results <- stats::setNames(lapply(skipped, function(label) {
    command <- plan$commands[[label]]
    list(
      label = label, status = 0L, success = TRUE, skipped = TRUE,
      command = .gpt_display_command(executable, command$args),
      log = "Skipped: complete output already exists."
    )
  }), skipped)
  results <- c(executed, skipped_results)[labels]
  failed <- labels[!vapply(results, function(x) isTRUE(x$success), logical(1))]
  if (length(failed) > 0L) {
    warning("PLINK2 command(s) failed: ", paste(failed, collapse = ", "), call. = FALSE)
  }
  execution <- list(tool = "PLINK2", results = results, plan = plan)
  class(execution) <- c("ukb_plink2_run", "list")
  success <- all(vapply(execution$results, function(x) isTRUE(x$success), logical(1)))
  scores <- if (isTRUE(read_scores) && success) {
    read_prs_scores(
      plan$expected_sscore,
      score_name = plan$settings$score_names,
      standardize = standardize
    )
  } else {
    NULL
  }
  out <- list(
    scores = scores, execution = execution, plan = plan,
    workers = as.integer(workers), resumed = isTRUE(resume),
    skipped = skipped
  )
  class(out) <- c("ukb_prs_run", "list")
  out
}

#' @export
print.ukb_prs_data <- function(x, ...) {
  cat("UKB published PRS data\n")
  cat("  Participants: ", nrow(x), "\n", sep = "")
  cat("  PRS columns: ", sum(grepl("^prs_", names(x))), "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_prs_weights <- function(x, ...) {
  cat("Prepared PRS weights\n")
  cat("  Variants: ", nrow(x), "\n", sep = "")
  cat("  File: ", if (is.null(attr(x, "path"))) "not written" else attr(x, "path"), "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_prs_plan <- function(x, ...) {
  cat("UKB PLINK2 PRS plan\n")
  cat("  Commands: ", length(x$commands), "\n", sep = "")
  cat("  Score(s): ", paste(x$settings$score_names, collapse = ", "), "\n", sep = "")
  cat("  Genotype format: ", toupper(x$settings$genotype_format), "\n", sep = "")
  cat("  Expected .sscore files: ", length(x$expected_sscore), "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_prs_scores <- function(x, ...) {
  cat("PLINK2 PRS scores\n")
  cat("  Participants: ", nrow(x), "\n", sep = "")
  n_files <- if (nrow(x) == 0L) 0L else max(x$N_FILES, na.rm = TRUE)
  cat("  Chromosome/files combined: ", n_files, "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_prs_run <- function(x, ...) {
  cat("UKB PLINK2 PRS run\n")
  cat("  Commands: ", length(x$execution$results), "\n", sep = "")
  cat("  Scores loaded: ", !is.null(x$scores), "\n", sep = "")
  invisible(x)
}

.prs_find_id <- function(columns, id_col) {
  candidates <- unique(c(id_col, paste0("participant.", id_col)))
  hit <- intersect(candidates, columns)
  if (length(hit) != 1L) {
    stop(sprintf("Could not identify one participant ID column `%s`.", id_col), call. = FALSE)
  }
  hit
}

.prs_find_field <- function(columns, field_id, required = TRUE) {
  pattern <- paste0("(^|\\.)p", field_id, "(?:_i[0-9]+)?(?:_a[0-9]+)?$")
  hit <- columns[grepl(pattern, columns, perl = TRUE)]
  if (length(hit) == 0L && isTRUE(required)) {
    stop(sprintf("UKB field %s was not found in `data`.", field_id), call. = FALSE)
  }
  if (length(hit) > 1L) {
    stop(sprintf("UKB field %s matched multiple columns; select one explicitly first.", field_id), call. = FALSE)
  }
  if (length(hit) == 0L) NULL else hit
}

.prs_resolve_column <- function(x, data) {
  if (is.numeric(x) && length(x) == 1L && !is.na(x) && x == as.integer(x) &&
      x >= 1L && x <= ncol(data)) {
    return(names(data)[[as.integer(x)]])
  }
  if (is.character(x) && length(x) == 1L && !is.na(x) && x %in% names(data)) {
    return(x)
  }
  stop("PRS column must be one valid name or one-based position.", call. = FALSE)
}

.prs_optional_column <- function(x, data) {
  if (is.null(x)) {
    return(NULL)
  }
  .prs_resolve_column(x, data)
}

.prs_positive_integer <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1L || x != as.integer(x)) {
    stop("PRS score columns must be positive integers.", call. = FALSE)
  }
  as.integer(x)
}

.prs_score_name <- function(x) {
  reserved <- c("FID", "IID", "ALLELE_CT", "DENOM", "N_FILES", "PRS_AVG")
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[A-Za-z][A-Za-z0-9_.-]*$", x) || x %in% reserved) {
    stop("`score_name` must be one simple, non-reserved name.", call. = FALSE)
  }
  x
}

.prs_zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

.prs_sum_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

.prs_testing_eligible <- function(x) {
  !is.na(x) & tolower(trimws(as.character(x))) %in% c("1", "true")
}

.prs_genome_build <- function(x, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop("Genome build must be `GRCh37` or `GRCh38`.", call. = FALSE)
  }
  key <- gsub("[^a-z0-9]", "", tolower(x))
  builds <- c(grch37 = "GRCh37", hg19 = "GRCh37", b37 = "GRCh37",
              grch38 = "GRCh38", hg38 = "GRCh38", b38 = "GRCh38")
  if (!key %in% names(builds)) {
    stop("Genome build must be `GRCh37` or `GRCh38`.", call. = FALSE)
  }
  unname(builds[[key]])
}

.prs_normalize_chr <- function(x) {
  x <- toupper(sub("^CHR", "", trimws(as.character(x)), ignore.case = TRUE))
  if (any(!grepl("^([1-9]|1[0-9]|2[0-2]|X|Y|XY|MT)$", x))) {
    stop("Chromosomes must use 1-22, X, Y, XY, or MT.", call. = FALSE)
  }
  x
}

.prs_column_range <- function(x) {
  x <- as.integer(x)
  if (length(x) == 0L || anyNA(x) || any(diff(x) != 1L)) {
    stop("Multi-PRS score columns must form one contiguous range.", call. = FALSE)
  }
  if (length(x) == 1L) as.character(x) else paste0(x[[1L]], "-", x[[length(x)]])
}

.prs_genotype_spec <- function(paths, format, sample, ref_first) {
  paths <- .gpt_paths(paths, "genotype")
  inferred <- if (identical(format, "auto")) {
    vapply(paths, .prs_infer_genotype_format, character(1))
  } else {
    rep(format, length(paths))
  }
  if (length(unique(inferred)) != 1L) {
    stop("All genotype inputs in one PRS plan must use the same format.", call. = FALSE)
  }
  format <- inferred[[1L]]
  labels <- .gpt_step_labels(paths)
  if (identical(format, "bgen")) {
    sample <- .gpt_recycle_paths(sample, length(paths), "sample")
    args <- lapply(seq_along(paths), function(i) {
      c("--bgen", paths[[i]], if (isTRUE(ref_first)) "ref-first" else character(),
        "--sample", sample[[i]])
    })
    validation <- lapply(seq_along(paths), function(i) c(paths[[i]], sample[[i]]))
    return(list(
      paths = paths, format = format, sample = sample, labels = labels,
      args = args, validation = validation
    ))
  }
  if (!is.null(sample)) {
    stop("`sample` must be NULL for PGEN or BED input.", call. = FALSE)
  }
  if (identical(format, "pgen")) {
    prefix <- sub("\\.(pgen|pvar(?:\\.zst)?|psam)$", "", paths, ignore.case = TRUE, perl = TRUE)
    compressed <- grepl("\\.pvar\\.zst$", paths, ignore.case = TRUE) |
      file.exists(paste0(prefix, ".pvar.zst"))
    args <- lapply(seq_along(prefix), function(i) {
      c("--pfile", prefix[[i]], if (compressed[[i]]) "vzs" else character())
    })
    validation <- lapply(seq_along(prefix), function(i) {
      c(paste0(prefix[[i]], ".pgen"),
        paste0(prefix[[i]], if (compressed[[i]]) ".pvar.zst" else ".pvar"),
        paste0(prefix[[i]], ".psam"))
    })
  } else {
    prefix <- sub("\\.(bed|bim|fam)$", "", paths, ignore.case = TRUE)
    args <- lapply(prefix, function(x) c("--bfile", x))
    validation <- lapply(prefix, function(x) paste0(x, c(".bed", ".bim", ".fam")))
  }
  list(
    paths = paths, format = format, sample = NULL, labels = labels,
    args = args, validation = validation
  )
}

.prs_infer_genotype_format <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.bgen$", lower)) return("bgen")
  if (grepl("\\.(pgen|pvar|pvar\\.zst|psam)$", lower)) return("pgen")
  if (grepl("\\.(bed|bim|fam)$", lower)) return("bed")
  if (file.exists(paste0(path, ".pgen"))) return("pgen")
  if (file.exists(paste0(path, ".bed"))) return("bed")
  stop("Could not infer genotype format; set `genotype_format` explicitly.", call. = FALSE)
}
