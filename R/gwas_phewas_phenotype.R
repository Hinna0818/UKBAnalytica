#' UKB Genetic Sample QC Profiles
#'
#' @description
#' Create a declarative sample-QC profile for GWAS phenotype preparation.
#' The default profile does not remove participants. The \code{"rap_example"}
#' profile reproduces the sample filters used in the DNAnexus UKB RAP
#' end-to-end GWAS/PheWAS tutorial: concordant reported/genetic sex, genetic
#' ethnic grouping value 1, no reported sex-chromosome aneuploidy, and
#' inclusion in the genetic principal-component calculation set.
#'
#' @param profile QC profile. \code{"none"} performs no built-in filtering;
#'   \code{"rap_example"} reproduces the official tutorial profile.
#' @param ancestry_values Values accepted in UKB field 22006 when
#'   \code{profile = "rap_example"}. Defaults to 1.
#' @param require_pca_sample Logical. Require field 22020 to equal 1 in the
#'   tutorial profile.
#'
#' @return A list of class \code{ukb_gwas_qc_profile}.
#' @export
ukb_gwas_qc_profile <- function(profile = c("none", "rap_example"),
                                ancestry_values = 1,
                                require_pca_sample = TRUE) {
  profile <- match.arg(profile)
  .gpp_assert_flag(require_pca_sample, "require_pca_sample")

  if (identical(profile, "none")) {
    out <- list(
      name = "none",
      sex_concordance = NULL,
      keep_values = list(),
      require_missing = character()
    )
  } else {
    if (length(ancestry_values) == 0L || anyNA(ancestry_values)) {
      stop("`ancestry_values` must contain at least one non-missing value.", call. = FALSE)
    }
    keep_values <- list(p22006 = ancestry_values)
    if (isTRUE(require_pca_sample)) {
      keep_values$p22020 <- 1
    }
    out <- list(
      name = "rap_example",
      sex_concordance = c("p31", "p22001"),
      keep_values = keep_values,
      require_missing = "p22019"
    )
  }

  class(out) <- c("ukb_gwas_qc_profile", class(out))
  out
}

#' Build a GWAS and PheWAS Phenotype Object
#'
#' @description
#' Prepare aligned, in-memory phenotype inputs for a UK Biobank GWAS followed
#' by PheWAS. The function produces a REGENIE-compatible wide table with
#' \code{FID}/\code{IID}, one or more phenotypes, and covariates, plus the
#' four-column long diagnosis table used by the PheWAS R package
#' (\code{id}, \code{vocabulary_id}, \code{code}, \code{count}).
#'
#' No files are written and no participant data are cached by the package.
#' Use \code{\link{ukb_write_gwas_phewas_phenotype}} explicitly inside RAP to
#' materialize analysis files in RAP-controlled storage.
#'
#' @param data Participant-level data.frame or data.table with one row per
#'   participant.
#' @param phenotype_cols One or more phenotype columns. Binary phenotypes must
#'   use 0/1/NA coding; quantitative phenotypes must be numeric.
#' @param diagnoses Long diagnosis table. It should normally be the result of
#'   \code{\link{parse_icd10_diagnoses}} or an equivalent table containing a
#'   participant ID and diagnosis code.
#' @param covariates Covariate columns to retain in the GWAS/PheWAS inputs.
#' @param categorical_covariates Subset of \code{covariates} that should be
#'   treated as categorical by downstream GWAS tools.
#' @param id_col Participant identifier column in \code{data}.
#' @param family_id_col Optional family identifier column. If \code{NULL},
#'   \code{FID} is set equal to \code{IID}, as in the RAP tutorial.
#' @param diagnosis_id_col Participant identifier column in \code{diagnoses}.
#' @param diagnosis_code_col Diagnosis-code column in \code{diagnoses}.
#' @param diagnosis_date_col Optional diagnosis-date column. When supplied,
#'   distinct non-missing dates are counted per participant-code pair.
#' @param diagnosis_count_col Optional existing numeric count column. It cannot
#'   be used together with \code{diagnosis_date_col}.
#' @param vocabulary_id Vocabulary label expected by PheWAS. For UKB field
#'   41270 use \code{"ICD10"}.
#' @param genotype_ids Optional vector of participant IDs present in the
#'   genotype files. When supplied, phenotype participants are intersected
#'   with this vector.
#' @param sample_qc A profile returned by
#'   \code{\link{ukb_gwas_qc_profile}}, or a compatible custom list.
#' @param trait_types Optional named character vector with one value per
#'   phenotype: \code{"binary"} or \code{"quantitative"}. If \code{NULL}, types
#'   are inferred conservatively.
#' @param covariate_missing How to handle missing covariates:
#'   \code{"keep"} preserves missing values, \code{"complete_case"} removes
#'   affected participants, and \code{"mean_mode"} performs explicit
#'   mean/mode imputation.
#' @param sex_col Optional sex column used to construct PheWAS sex
#'   restrictions.
#' @param sex_map Named character vector mapping values in \code{sex_col} to
#'   \code{"F"} or \code{"M"}. The default supports UKB field 31 coding and
#'   already formatted values.
#' @param min_code_count Minimum retained diagnosis count per
#'   participant-code pair. This filters the long code table only; the
#'   downstream PheWAS case threshold remains independently configurable.
#'
#' @return An object of class \code{ukb_gwas_phewas_phenotype}.
#' @export
build_gwas_phewas_phenotype <- function(
    data,
    phenotype_cols,
    diagnoses,
    covariates = character(),
    categorical_covariates = character(),
    id_col = "eid",
    family_id_col = NULL,
    diagnosis_id_col = id_col,
    diagnosis_code_col = "icd10_code",
    diagnosis_date_col = NULL,
    diagnosis_count_col = NULL,
    vocabulary_id = "ICD10",
    genotype_ids = NULL,
    sample_qc = ukb_gwas_qc_profile("none"),
    trait_types = NULL,
    covariate_missing = c("keep", "complete_case", "mean_mode"),
    sex_col = NULL,
    sex_map = c("0" = "F", "1" = "M", "F" = "F", "M" = "M"),
    min_code_count = 1L) {
  covariate_missing <- match.arg(covariate_missing)
  phenotype_cols <- .gpp_character_vector(phenotype_cols, "phenotype_cols")
  covariates <- .gpp_character_vector(covariates, "covariates", allow_empty = TRUE)
  categorical_covariates <- .gpp_character_vector(
    categorical_covariates,
    "categorical_covariates",
    allow_empty = TRUE
  )
  if (!all(categorical_covariates %in% covariates)) {
    stop("`categorical_covariates` must be a subset of `covariates`.", call. = FALSE)
  }
  if (length(intersect(phenotype_cols, covariates)) > 0L) {
    stop("Phenotype columns and covariates must be distinct.", call. = FALSE)
  }
  if (!is.numeric(min_code_count) || length(min_code_count) != 1L ||
      is.na(min_code_count) || min_code_count < 1 ||
      min_code_count != as.integer(min_code_count)) {
    stop("`min_code_count` must be one positive integer.", call. = FALSE)
  }
  min_code_count <- as.integer(min_code_count)
  if (!is.null(diagnosis_date_col) && !is.null(diagnosis_count_col)) {
    stop("Use only one of `diagnosis_date_col` and `diagnosis_count_col`.", call. = FALSE)
  }
  if (!is.character(vocabulary_id) || length(vocabulary_id) != 1L ||
      is.na(vocabulary_id) || !nzchar(vocabulary_id)) {
    stop("`vocabulary_id` must be one non-empty string.", call. = FALSE)
  }

  participant <- data.table::as.data.table(data.table::copy(data))
  required <- unique(c(
    id_col, family_id_col, phenotype_cols, covariates, sex_col,
    .gpp_qc_columns(sample_qc)
  ))
  required <- required[!is.na(required) & nzchar(required)]
  .gpp_require_columns(participant, required, "data")
  if (anyDuplicated(participant[[id_col]])) {
    stop("`data` must contain one row per participant ID.", call. = FALSE)
  }
  if (anyNA(participant[[id_col]])) {
    stop("Participant IDs in `data` cannot be missing.", call. = FALSE)
  }

  qc <- .gpp_apply_sample_qc(
    data = participant,
    id_col = id_col,
    sample_qc = sample_qc,
    genotype_ids = genotype_ids
  )
  participant <- qc$data

  missing_result <- .gpp_handle_covariate_missing(
    data = participant,
    covariates = covariates,
    method = covariate_missing
  )
  participant <- missing_result$data
  if (nrow(participant) == 0L) {
    stop("No participants remain after sample QC and missing-data handling.", call. = FALSE)
  }

  trait_types <- .gpp_trait_types(participant, phenotype_cols, trait_types)
  participant <- .gpp_normalize_traits(participant, phenotype_cols, trait_types)

  iid <- participant[[id_col]]
  fid <- if (is.null(family_id_col)) iid else participant[[family_id_col]]
  if (anyNA(fid)) {
    stop("Family IDs cannot be missing.", call. = FALSE)
  }

  gwas <- data.table::data.table(FID = fid, IID = iid)
  gwas <- cbind(
    gwas,
    participant[, c(phenotype_cols, covariates), with = FALSE]
  )
  data.table::setDF(gwas)

  covariate_table <- data.table::data.table(FID = fid, IID = iid)
  if (length(covariates) > 0L) {
    covariate_table <- cbind(
      covariate_table,
      participant[, covariates, with = FALSE]
    )
  }
  data.table::setDF(covariate_table)

  sample_ids <- data.frame(FID = fid, IID = iid, stringsAsFactors = FALSE)
  sex <- .gpp_prepare_sex(participant, id_col, sex_col, sex_map)
  phewas_long <- .gpp_prepare_phewas_long(
    diagnoses = diagnoses,
    included_ids = iid,
    diagnosis_id_col = diagnosis_id_col,
    diagnosis_code_col = diagnosis_code_col,
    diagnosis_date_col = diagnosis_date_col,
    diagnosis_count_col = diagnosis_count_col,
    vocabulary_id = vocabulary_id,
    min_code_count = min_code_count
  )
  trait_summary <- .gpp_trait_summary(participant, phenotype_cols, trait_types)
  unsupported_binary <- trait_summary$trait[
    trait_summary$trait_type == "binary" &
      (trait_summary$n_case == 0L | trait_summary$n_control == 0L)
  ]
  if (length(unsupported_binary) > 0L) {
    warning(
      sprintf(
        paste(
          "Binary phenotype(s) without both cases and controls cannot be",
          "tested: %s"
        ),
        paste(unsupported_binary, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- list(
    gwas = gwas,
    covariates = covariate_table,
    phewas_long = phewas_long,
    sample_ids = sample_ids,
    sex = sex,
    trait_summary = trait_summary,
    qc_flow = qc$flow,
    covariate_missing = missing_result$summary,
    settings = list(
      id_col = id_col,
      family_id_col = family_id_col,
      phenotype_cols = phenotype_cols,
      trait_types = trait_types,
      covariates = covariates,
      categorical_covariates = categorical_covariates,
      vocabulary_id = vocabulary_id,
      min_code_count = min_code_count,
      covariate_missing = covariate_missing,
      sample_qc = .gpp_or(sample_qc$name, "custom")
    )
  )
  class(out) <- c("ukb_gwas_phewas_phenotype", "list")
  out
}

#' Write a GWAS/PheWAS Phenotype Object
#'
#' @description
#' Materialize phenotype inputs in an explicit output directory. By default the
#' function can only write from a RAP-like session to a path under
#' \code{/mnt/project}. This keeps participant-level analysis files inside
#' RAP-controlled storage.
#'
#' @param x Object returned by \code{\link{build_gwas_phewas_phenotype}}.
#' @param output_dir Existing or new local directory. With
#'   \code{require_rap = TRUE}, it must be under \code{/mnt/project}.
#' @param prefix File prefix.
#' @param require_rap Require a RAP-like environment and RAP-scoped path.
#' @param overwrite Overwrite existing files.
#'
#' @return A list of class \code{ukb_gwas_phewas_files} containing generated
#'   file paths and a non-participant-level manifest.
#' @export
ukb_write_gwas_phewas_phenotype <- function(x,
                                            output_dir,
                                            prefix = "ukba_gwas_phewas",
                                            require_rap = TRUE,
                                            overwrite = FALSE) {
  if (!inherits(x, "ukb_gwas_phewas_phenotype")) {
    stop("`x` must be a `ukb_gwas_phewas_phenotype` object.", call. = FALSE)
  }
  .gpp_assert_flag(require_rap, "require_rap")
  .gpp_assert_flag(overwrite, "overwrite")
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("`output_dir` must be one non-empty path.", call. = FALSE)
  }
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", prefix)) {
    stop("`prefix` may contain letters, numbers, dot, underscore, and hyphen.", call. = FALSE)
  }

  output_dir <- path.expand(output_dir)
  if (isTRUE(require_rap)) {
    .ukb_assert_rap_env("ukb_write_gwas_phewas_phenotype()")
    if (!.ukb_output_path_is_rap_scoped(output_dir) ||
        !startsWith(normalizePath(output_dir, mustWork = FALSE), "/mnt/project")) {
      stop("Participant-level phenotype files must be written under `/mnt/project`.", call. = FALSE)
    }
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create `output_dir`.", call. = FALSE)
  }

  paths <- c(
    gwas = file.path(output_dir, paste0(prefix, ".gwas.phe")),
    covariates = file.path(output_dir, paste0(prefix, ".covariates.tsv")),
    phewas_long = file.path(output_dir, paste0(prefix, ".phewas_long.csv")),
    sample_keep = file.path(output_dir, paste0(prefix, ".samples.keep")),
    trait_summary = file.path(output_dir, paste0(prefix, ".trait_summary.tsv")),
    qc_flow = file.path(output_dir, paste0(prefix, ".qc_flow.tsv")),
    covariate_missing = file.path(output_dir, paste0(prefix, ".covariate_missing.tsv")),
    manifest = file.path(output_dir, paste0(prefix, ".manifest.tsv"))
  )
  if (!is.null(x$sex)) {
    paths <- c(paths, sex = file.path(output_dir, paste0(prefix, ".sex.tsv")))
  }
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L && !isTRUE(overwrite)) {
    stop(
      sprintf("Output file(s) already exist: %s", paste(basename(existing), collapse = ", ")),
      call. = FALSE
    )
  }

  data.table::fwrite(x$gwas, paths[["gwas"]], sep = "\t", na = "NA", quote = FALSE)
  data.table::fwrite(x$covariates, paths[["covariates"]], sep = "\t", na = "NA", quote = FALSE)
  data.table::fwrite(x$phewas_long, paths[["phewas_long"]], sep = ",", na = "NA", quote = TRUE)
  data.table::fwrite(
    x$sample_ids,
    paths[["sample_keep"]],
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
  )
  data.table::fwrite(x$trait_summary, paths[["trait_summary"]], sep = "\t", na = "NA")
  data.table::fwrite(x$qc_flow, paths[["qc_flow"]], sep = "\t", na = "NA")
  data.table::fwrite(
    x$covariate_missing,
    paths[["covariate_missing"]],
    sep = "\t",
    na = "NA"
  )
  if (!is.null(x$sex)) {
    data.table::fwrite(x$sex, paths[["sex"]], sep = "\t", na = "NA")
  }

  component_rows <- c(
    gwas = nrow(x$gwas),
    covariates = nrow(x$covariates),
    phewas_long = nrow(x$phewas_long),
    sample_keep = nrow(x$sample_ids),
    trait_summary = nrow(x$trait_summary),
    qc_flow = nrow(x$qc_flow),
    covariate_missing = nrow(x$covariate_missing)
  )
  component_cols <- c(
    gwas = ncol(x$gwas),
    covariates = ncol(x$covariates),
    phewas_long = ncol(x$phewas_long),
    sample_keep = ncol(x$sample_ids),
    trait_summary = ncol(x$trait_summary),
    qc_flow = ncol(x$qc_flow),
    covariate_missing = ncol(x$covariate_missing)
  )
  if (!is.null(x$sex)) {
    component_rows <- c(component_rows, sex = nrow(x$sex))
    component_cols <- c(component_cols, sex = ncol(x$sex))
  }
  manifest_names <- setdiff(names(paths), "manifest")
  manifest <- data.frame(
    component = manifest_names,
    file = unname(paths[manifest_names]),
    n_rows = unname(component_rows[manifest_names]),
    n_columns = unname(component_cols[manifest_names]),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(manifest, paths[["manifest"]], sep = "\t", quote = FALSE)

  out <- list(paths = paths, manifest = manifest, settings = x$settings)
  class(out) <- c("ukb_gwas_phewas_files", "list")
  out
}

#' @export
print.ukb_gwas_phewas_phenotype <- function(x, ...) {
  cat("UKB GWAS + PheWAS phenotype\n")
  cat("  Participants: ", nrow(x$sample_ids), "\n", sep = "")
  cat("  GWAS traits: ", paste(x$settings$phenotype_cols, collapse = ", "), "\n", sep = "")
  cat("  Covariates: ", length(x$settings$covariates), "\n", sep = "")
  cat("  PheWAS code rows: ", nrow(x$phewas_long), "\n", sep = "")
  cat("  QC profile: ", x$settings$sample_qc, "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_gwas_phewas_files <- function(x, ...) {
  cat("UKB GWAS + PheWAS files\n")
  cat("  Output files: ", length(x$paths), "\n", sep = "")
  cat("  Directory: ", dirname(unname(x$paths[[1L]])), "\n", sep = "")
  invisible(x)
}

.gpp_assert_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.gpp_or <- function(x, y) {
  if (is.null(x)) y else x
}

.gpp_character_vector <- function(x, name, allow_empty = FALSE) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) ||
      (!isTRUE(allow_empty) && length(x) == 0L)) {
    requirement <- if (isTRUE(allow_empty)) {
      "a character vector of non-empty names"
    } else {
      "a non-empty character vector of non-empty names"
    }
    stop(
      sprintf("`%s` must be %s.", name, requirement),
      call. = FALSE
    )
  }
  unique(x)
}

.gpp_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(
      sprintf("Missing column(s) in `%s`: %s", label, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.gpp_qc_columns <- function(sample_qc) {
  if (is.null(sample_qc)) {
    return(character())
  }
  if (!is.list(sample_qc)) {
    stop("`sample_qc` must be a QC profile or compatible list.", call. = FALSE)
  }
  c(
    .gpp_or(sample_qc$sex_concordance, character()),
    names(.gpp_or(sample_qc$keep_values, list())),
    .gpp_or(sample_qc$require_missing, character())
  )
}

.gpp_apply_sample_qc <- function(data, id_col, sample_qc, genotype_ids) {
  if (is.null(sample_qc)) {
    sample_qc <- ukb_gwas_qc_profile("none")
  }
  .gpp_qc_columns(sample_qc)
  out <- data.table::copy(data)
  flow <- data.table::data.table(
    step = "input",
    n_before = nrow(out),
    n_after = nrow(out),
    n_removed = 0L
  )

  apply_keep <- function(keep, label) {
    before <- nrow(out)
    keep[is.na(keep)] <- FALSE
    out <<- out[keep]
    flow <<- data.table::rbindlist(list(
      flow,
      data.table::data.table(
        step = label,
        n_before = before,
        n_after = nrow(out),
        n_removed = before - nrow(out)
      )
    ))
  }

  sex_columns <- .gpp_or(sample_qc$sex_concordance, character())
  if (length(sex_columns) > 0L) {
    if (length(sex_columns) != 2L) {
      stop("`sample_qc$sex_concordance` must contain two column names.", call. = FALSE)
    }
    .gpp_require_columns(out, sex_columns, "data")
    apply_keep(
      !is.na(out[[sex_columns[[1L]]]]) &
        !is.na(out[[sex_columns[[2L]]]]) &
        out[[sex_columns[[1L]]]] == out[[sex_columns[[2L]]]],
      "reported/genetic sex concordance"
    )
  }

  keep_values <- .gpp_or(sample_qc$keep_values, list())
  if (length(keep_values) > 0L && (is.null(names(keep_values)) || any(!nzchar(names(keep_values))))) {
    stop("`sample_qc$keep_values` must be a named list.", call. = FALSE)
  }
  for (column in names(keep_values)) {
    .gpp_require_columns(out, column, "data")
    values <- keep_values[[column]]
    if (length(values) == 0L || anyNA(values)) {
      stop(sprintf("QC keep values for `%s` cannot be empty or missing.", column), call. = FALSE)
    }
    apply_keep(
      !is.na(out[[column]]) & out[[column]] %in% values,
      paste0(column, " in {", paste(values, collapse = ","), "}")
    )
  }

  require_missing <- .gpp_or(sample_qc$require_missing, character())
  for (column in require_missing) {
    .gpp_require_columns(out, column, "data")
    apply_keep(is.na(out[[column]]), paste0(column, " is missing"))
  }

  if (!is.null(genotype_ids)) {
    genotype_ids <- unique(as.character(genotype_ids[!is.na(genotype_ids)]))
    apply_keep(
      as.character(out[[id_col]]) %in% genotype_ids,
      "present in genotype data"
    )
  }

  list(data = out, flow = as.data.frame(flow))
}

.gpp_handle_covariate_missing <- function(data, covariates, method) {
  out <- data.table::copy(data)
  if (length(covariates) == 0L) {
    return(list(
      data = out,
      summary = data.frame(
        covariate = character(),
        n_missing_before = integer(),
        action = character(),
        fill_value = character(),
        stringsAsFactors = FALSE
      )
    ))
  }

  summary <- data.frame(
    covariate = covariates,
    n_missing_before = vapply(covariates, function(x) sum(is.na(out[[x]])), integer(1)),
    action = method,
    fill_value = NA_character_,
    stringsAsFactors = FALSE
  )
  if (identical(method, "complete_case")) {
    out <- out[stats::complete.cases(out[, covariates, with = FALSE])]
  } else if (identical(method, "mean_mode")) {
    for (i in seq_along(covariates)) {
      column <- covariates[[i]]
      values <- out[[column]]
      missing <- is.na(values)
      if (!any(missing)) {
        summary$fill_value[[i]] <- ""
        next
      }
      observed <- values[!missing]
      if (length(observed) == 0L) {
        warning(sprintf("Covariate `%s` is entirely missing and was not imputed.", column), call. = FALSE)
        next
      }
      fill <- if (is.numeric(values)) {
        mean(observed)
      } else {
        tab <- sort(table(as.character(observed)), decreasing = TRUE)
        names(tab)[[1L]]
      }
      if (is.factor(values)) {
        replacement <- as.character(values)
        replacement[missing] <- as.character(fill)
        out[[column]] <- factor(replacement, levels = union(levels(values), as.character(fill)))
      } else {
        replacement <- if (is.numeric(values)) as.numeric(values) else values
        replacement[missing] <- fill
        data.table::set(out, j = column, value = replacement)
      }
      summary$fill_value[[i]] <- as.character(fill)
    }
  }
  list(data = out, summary = summary)
}

.gpp_trait_types <- function(data, phenotype_cols, trait_types) {
  if (is.null(trait_types)) {
    inferred <- vapply(phenotype_cols, function(column) {
      values <- data[[column]]
      observed <- values[!is.na(values)]
      if (length(observed) == 0L) {
        stop(sprintf("Phenotype `%s` is entirely missing.", column), call. = FALSE)
      }
      if ((is.numeric(observed) || is.logical(observed)) &&
          all(as.numeric(observed) %in% c(0, 1))) {
        "binary"
      } else if (is.numeric(observed)) {
        "quantitative"
      } else {
        stop(
          sprintf(
            "Could not infer GWAS trait type for `%s`; provide `trait_types`.",
            column
          ),
          call. = FALSE
        )
      }
    }, character(1))
    names(inferred) <- phenotype_cols
    return(inferred)
  }

  if (!is.character(trait_types) || anyNA(trait_types)) {
    stop("`trait_types` must be a character vector.", call. = FALSE)
  }
  if (is.null(names(trait_types))) {
    if (length(trait_types) != length(phenotype_cols)) {
      stop("Unnamed `trait_types` must match `phenotype_cols` in length.", call. = FALSE)
    }
    names(trait_types) <- phenotype_cols
  }
  if (!all(phenotype_cols %in% names(trait_types))) {
    stop("Named `trait_types` must include every phenotype column.", call. = FALSE)
  }
  trait_types <- trait_types[phenotype_cols]
  if (!all(trait_types %in% c("binary", "quantitative"))) {
    stop("Trait types must be `binary` or `quantitative`.", call. = FALSE)
  }
  trait_types
}

.gpp_normalize_traits <- function(data, phenotype_cols, trait_types) {
  out <- data.table::copy(data)
  for (column in phenotype_cols) {
    values <- out[[column]]
    if (identical(unname(trait_types[[column]]), "binary")) {
      observed <- values[!is.na(values)]
      if (!(is.numeric(values) || is.logical(values)) ||
          !all(as.numeric(observed) %in% c(0, 1))) {
        stop(sprintf("Binary phenotype `%s` must use 0/1/NA coding.", column), call. = FALSE)
      }
      data.table::set(out, j = column, value = as.integer(values))
    } else if (!is.numeric(values)) {
      stop(sprintf("Quantitative phenotype `%s` must be numeric.", column), call. = FALSE)
    }
  }
  out
}

.gpp_trait_summary <- function(data, phenotype_cols, trait_types) {
  rows <- lapply(phenotype_cols, function(column) {
    values <- data[[column]]
    observed <- values[!is.na(values)]
    binary <- identical(unname(trait_types[[column]]), "binary")
    data.frame(
      trait = column,
      trait_type = unname(trait_types[[column]]),
      n_total = length(values),
      n_nonmissing = length(observed),
      n_missing = sum(is.na(values)),
      n_case = if (binary) sum(observed == 1) else NA_integer_,
      n_control = if (binary) sum(observed == 0) else NA_integer_,
      prevalence = if (binary && length(observed) > 0L) mean(observed == 1) else NA_real_,
      mean = if (!binary && length(observed) > 0L) mean(observed) else NA_real_,
      sd = if (!binary && length(observed) > 1L) stats::sd(observed) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.gpp_prepare_sex <- function(data, id_col, sex_col, sex_map) {
  if (is.null(sex_col)) {
    return(NULL)
  }
  if (!is.character(sex_map) || is.null(names(sex_map)) || anyNA(sex_map) ||
      any(!sex_map %in% c("F", "M"))) {
    stop("`sex_map` must be a named character vector mapping to `F` or `M`.", call. = FALSE)
  }
  raw <- as.character(data[[sex_col]])
  mapped <- unname(sex_map[raw])
  data.frame(id = data[[id_col]], sex = mapped, stringsAsFactors = FALSE)
}

.gpp_prepare_phewas_long <- function(diagnoses,
                                     included_ids,
                                     diagnosis_id_col,
                                     diagnosis_code_col,
                                     diagnosis_date_col,
                                     diagnosis_count_col,
                                     vocabulary_id,
                                     min_code_count) {
  diagnosis <- data.table::as.data.table(data.table::copy(diagnoses))
  required <- c(diagnosis_id_col, diagnosis_code_col)
  if (!is.null(diagnosis_date_col)) {
    required <- c(required, diagnosis_date_col)
  }
  if (!is.null(diagnosis_count_col)) {
    required <- c(required, diagnosis_count_col)
  }
  .gpp_require_columns(diagnosis, unique(required), "diagnoses")

  diagnosis <- diagnosis[
    as.character(get(diagnosis_id_col)) %in% as.character(included_ids)
  ]
  diagnosis[
    ,
    ".gpp_code" := toupper(trimws(as.character(get(diagnosis_code_col))))
  ]
  diagnosis <- diagnosis[
    !is.na(get(".gpp_code")) & nzchar(get(".gpp_code"))
  ]

  if (nrow(diagnosis) == 0L) {
    return(data.frame(
      id = included_ids[0],
      vocabulary_id = character(),
      code = character(),
      count = integer(),
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(diagnosis_count_col)) {
    counts <- suppressWarnings(as.numeric(diagnosis[[diagnosis_count_col]]))
    if (any(is.na(counts) & !is.na(diagnosis[[diagnosis_count_col]])) ||
        any(counts < 0, na.rm = TRUE)) {
      stop("`diagnosis_count_col` must contain non-negative numeric counts.", call. = FALSE)
    }
    diagnosis[, ".gpp_count" := counts]
    long <- diagnosis[
      ,
      .(count = sum(get(".gpp_count"), na.rm = TRUE)),
      by = c(diagnosis_id_col, ".gpp_code")
    ]
  } else if (!is.null(diagnosis_date_col)) {
    long <- diagnosis[
      ,
      .(count = {
        values <- get(diagnosis_date_col)
        nonmissing <- values[!is.na(values)]
        if (length(nonmissing) == 0L) .N else data.table::uniqueN(nonmissing)
      }),
      by = c(diagnosis_id_col, ".gpp_code")
    ]
  } else {
    long <- diagnosis[
      ,
      .(count = .N),
      by = c(diagnosis_id_col, ".gpp_code")
    ]
  }

  long <- long[get("count") >= min_code_count]
  data.table::setnames(long, c(diagnosis_id_col, ".gpp_code"), c("id", "code"))
  vocab_value <- vocabulary_id
  long[, "vocabulary_id" := vocab_value]
  data.table::setcolorder(long, c("id", "vocabulary_id", "code", "count"))
  data.table::setorderv(long, c("id", "code"))
  as.data.frame(long)
}
