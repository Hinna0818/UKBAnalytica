#' Curated GWAS and PheWAS Tool Set for UKB RAP
#'
#' @description
#' Return the tool set considered by UKBAnalytica for UK Biobank GWAS/PheWAS
#' workflows. The table distinguishes the first supported path from optional
#' alternatives. It is metadata only and does not install or execute software.
#'
#' @return A data.frame describing tool roles and support status.
#' @export
ukb_gwas_phewas_tools <- function() {
  data.frame(
    tool = c("REGENIE", "PLINK2", "PheWAS", "SAIGE", "PHESANT", "DeepPheWAS", "BOLT-LMM"),
    stage = c(
      "GWAS association",
      "genotype QC/export/LD",
      "ICD-to-phecode PheWAS",
      "GWAS association",
      "general UKB phenome scan",
      "deep phenotype PheWAS",
      "GWAS association"
    ),
    support = c(
      "primary",
      "primary",
      "primary_optional",
      "planned_alternative",
      "interoperable_input",
      "interoperable_input",
      "not_planned"
    ),
    rap_mode = c(
      "RAP REGENIE app or Swiss Army Knife",
      "Swiss Army Knife or RAP shell",
      "R in RAP",
      "custom RAP environment",
      "R in RAP",
      "R plus PLINK2",
      "custom environment"
    ),
    rationale = c(
      "Official RAP GWAS path; scalable multi-trait two-step analysis.",
      "Standard UKB genotype QC, conversion, extraction, and LD utility.",
      "Used by the official RAP GWAS/PheWAS tutorial for ICD10 phecodes.",
      "Strong option for severe case-control imbalance and related samples.",
      "Designed for heterogeneous UKB data fields rather than ICD-only PheWAS.",
      "Rich phenotype generation but less established and adds another framework.",
      "Highly cited, especially for quantitative traits, but not the RAP tutorial path."
    ),
    stringsAsFactors = FALSE
  )
}

#' Plan an Arbitrary REGENIE Command
#'
#' @description
#' Create an auditable, non-executing REGENIE command plan from raw
#' command-line arguments. Arguments are preserved as individual tokens, so
#' this escape hatch can represent official REGENIE options that are not
#' exposed by the structured \code{\link{ukb_plan_regenie}} wrapper.
#'
#' UKBAnalytica validates token boundaries, but it does not attempt to validate
#' every option or option combination. The installed REGENIE version remains
#' authoritative.
#'
#' @param args Character vector containing REGENIE command-line argument
#'   tokens in execution order.
#' @param expected_outputs Optional character vector of expected output paths.
#'   This is metadata only; REGENIE remains responsible for creating them.
#' @param label Command label used in execution results.
#'
#' @return A list of class \code{ukb_regenie_command_plan} that also inherits
#'   from \code{ukb_regenie_plan}.
#' @export
ukb_plan_regenie_command <- function(args,
                                     expected_outputs = character(),
                                     label = "regenie") {
  args <- .gpt_cli_args(args, "args")
  expected_outputs <- .gpt_cli_args(
    expected_outputs,
    "expected_outputs",
    allow_empty = TRUE
  )
  label <- .gpt_command_label(label)

  commands <- list(.gpt_command("regenie", args))
  names(commands) <- label
  out <- list(
    tool = "REGENIE",
    commands = commands,
    expected_outputs = expected_outputs,
    settings = list(mode = "raw_args")
  )
  class(out) <- c("ukb_regenie_command_plan", "ukb_regenie_plan", "list")
  out
}

#' Plan a REGENIE GWAS
#'
#' @description
#' Build auditable REGENIE Step 1 and Step 2 command arguments from files
#' produced by \code{\link{ukb_write_gwas_phewas_phenotype}}. This function
#' does not submit a RAP job and does not execute REGENIE.
#'
#' The generated commands are intended for the REGENIE binary available in the
#' RAP Swiss Army Knife app or another controlled RAP execution environment.
#' For the managed RAP REGENIE orchestration app, use the same phenotype,
#' covariate, sample, and variant-list files through that app's versioned input
#' interface.
#'
#' @param files A \code{ukb_gwas_phewas_files} object or path to a
#'   REGENIE-compatible phenotype file.
#' @param bed_prefix Step 1 PLINK BED/BIM/FAM prefix, without extension.
#' @param step1_extract Step 1 variant inclusion file.
#' @param step2_bgen One or more Step 2 BGEN paths.
#' @param step2_sample One sample path or one path per BGEN.
#' @param output_prefix Output prefix in RAP-controlled storage.
#' @param phenotypes Phenotype columns. Automatically obtained when
#'   \code{files} is a UKBAnalytica files object.
#' @param covariate_file Optional covariate file. Defaults to the generated
#'   REGENIE-compatible covariate file.
#' @param covariates Covariate columns.
#' @param categorical_covariates Subset of covariates treated as categorical.
#' @param trait_type \code{"binary"} or \code{"quantitative"}. Automatically
#'   obtained when all selected traits have the same type.
#' @param keep_file Optional two-column FID/IID keep file. Defaults to the
#'   generated sample file.
#' @param step2_extract Optional Step 2 variant inclusion file. Supply one path
#'   or one path per BGEN.
#' @param prediction_file Step 1 prediction-list path expected by Step 2. If
#'   \code{NULL}, uses REGENIE's conventional
#'   \code{<output_prefix>_step1_pred.list} path.
#' @param ref_first Add \code{--ref-first} to Step 2.
#' @param firth Use Firth fallback for binary traits.
#' @param spa Use SPA instead of Firth fallback for binary traits. Only one of
#'   \code{firth} and \code{spa} may be true.
#' @param p_threshold P-value threshold for approximate Firth/SPA fallback.
#' @param bsize_step1,bsize_step2 REGENIE block sizes.
#' @param threads Optional positive integer thread count.
#' @param lowmem Use REGENIE Step 1 low-memory mode.
#' @param step1_args,step2_args Additional REGENIE argument tokens appended to
#'   the corresponding step. Flags managed by this wrapper are rejected; use
#'   \code{\link{ukb_plan_regenie_command}} for a fully custom command.
#'
#' @return A list of class \code{ukb_regenie_plan}.
#' @export
ukb_plan_regenie <- function(
    files,
    bed_prefix,
    step1_extract,
    step2_bgen,
    step2_sample,
    output_prefix,
    phenotypes = NULL,
    covariate_file = NULL,
    covariates = NULL,
    categorical_covariates = NULL,
    trait_type = NULL,
    keep_file = NULL,
    step2_extract = NULL,
    prediction_file = NULL,
    ref_first = TRUE,
    firth = TRUE,
    spa = FALSE,
    p_threshold = 0.01,
    bsize_step1 = 1000L,
    bsize_step2 = 400L,
    threads = NULL,
    lowmem = TRUE,
    step1_args = character(),
    step2_args = character()) {
  .gpp_assert_flag(ref_first, "ref_first")
  .gpp_assert_flag(firth, "firth")
  .gpp_assert_flag(spa, "spa")
  .gpp_assert_flag(lowmem, "lowmem")
  if (isTRUE(firth) && isTRUE(spa)) {
    stop("Use only one of `firth` and `spa`.", call. = FALSE)
  }
  .gpt_positive_number(p_threshold, "p_threshold", allow_zero = TRUE)
  .gpt_positive_integer(bsize_step1, "bsize_step1")
  .gpt_positive_integer(bsize_step2, "bsize_step2")
  if (!is.null(threads)) {
    .gpt_positive_integer(threads, "threads")
  }
  step1_extra_args <- .gpt_cli_args(
    step1_args,
    "step1_args",
    allow_empty = TRUE
  )
  step2_extra_args <- .gpt_cli_args(
    step2_args,
    "step2_args",
    allow_empty = TRUE
  )
  .gpt_reject_managed_flags(
    step1_extra_args,
    c(
      "--step", "--bed", "--bgen", "--pgen",
      "--phenoFile", "--phenoCol", "--phenoColList",
      "--covarFile", "--covarCol", "--covarColList", "--catCovarList",
      "--keep", "--extract",
      "--bt", "--qt", "--t2e", "--force-qt", "-1", "--cc12",
      "--bsize", "--lowmem", "--lowmem-prefix", "--threads", "--out"
    ),
    "ukb_plan_regenie(step1_args=)",
    tool = "REGENIE",
    generic = "ukb_plan_regenie_command()"
  )
  .gpt_reject_managed_flags(
    step2_extra_args,
    c(
      "--step", "--bed", "--bgen", "--pgen", "--sample",
      "--phenoFile", "--phenoCol", "--phenoColList",
      "--covarFile", "--covarCol", "--covarColList", "--catCovarList",
      "--keep", "--extract", "--ref-first", "--pred",
      "--bt", "--qt", "--t2e", "--force-qt", "-1", "--cc12",
      "--bsize", "--firth", "--approx", "--spa", "--pThresh",
      "--threads", "--out"
    ),
    "ukb_plan_regenie(step2_args=)",
    tool = "REGENIE",
    generic = "ukb_plan_regenie_command()"
  )

  resolved <- .gpt_resolve_files(
    files = files,
    phenotypes = phenotypes,
    covariate_file = covariate_file,
    covariates = covariates,
    categorical_covariates = categorical_covariates,
    trait_type = trait_type,
    keep_file = keep_file
  )
  phenotype_file <- resolved$phenotype_file
  phenotypes <- .gpt_names(resolved$phenotypes, "phenotypes")
  covariates <- .gpt_names(resolved$covariates, "covariates", allow_empty = TRUE)
  categorical_covariates <- .gpt_names(
    resolved$categorical_covariates,
    "categorical_covariates",
    allow_empty = TRUE
  )
  if (!all(categorical_covariates %in% covariates)) {
    stop("`categorical_covariates` must be a subset of `covariates`.", call. = FALSE)
  }
  trait_type <- match.arg(resolved$trait_type, c("binary", "quantitative"))

  bed_prefix <- .gpt_path(bed_prefix, "bed_prefix")
  step1_extract <- .gpt_path(step1_extract, "step1_extract")
  output_prefix <- .gpt_path(output_prefix, "output_prefix")
  step2_bgen <- .gpt_paths(step2_bgen, "step2_bgen")
  step2_sample <- .gpt_recycle_paths(step2_sample, length(step2_bgen), "step2_sample")
  if (!is.null(step2_extract)) {
    step2_extract <- .gpt_recycle_paths(step2_extract, length(step2_bgen), "step2_extract")
  }
  prediction_file <- if (is.null(prediction_file)) {
    paste0(output_prefix, "_step1_pred.list")
  } else {
    .gpt_path(prediction_file, "prediction_file")
  }

  common <- c(
    "--phenoFile", phenotype_file,
    "--phenoColList", paste(phenotypes, collapse = ",")
  )
  if (length(covariates) > 0L) {
    common <- c(
      common,
      "--covarFile", resolved$covariate_file,
      "--covarColList", paste(covariates, collapse = ",")
    )
  }
  if (length(categorical_covariates) > 0L) {
    common <- c(common, "--catCovarList", paste(categorical_covariates, collapse = ","))
  }
  trait_flag <- if (identical(trait_type, "binary")) "--bt" else "--qt"

  step1_command_args <- c(
    "--step", "1",
    "--bed", bed_prefix,
    "--extract", step1_extract
  )
  if (!is.null(resolved$keep_file)) {
    step1_command_args <- c(
      step1_command_args,
      "--keep", resolved$keep_file
    )
  }
  step1_command_args <- c(
    step1_command_args,
    common,
    trait_flag,
    "--bsize", as.character(as.integer(bsize_step1))
  )
  if (isTRUE(lowmem)) {
    step1_command_args <- c(
      step1_command_args,
      "--lowmem",
      "--lowmem-prefix", paste0(output_prefix, "_step1_tmp")
    )
  }
  if (!is.null(threads)) {
    step1_command_args <- c(
      step1_command_args,
      "--threads", as.character(as.integer(threads))
    )
  }
  step1_command_args <- c(
    step1_command_args,
    step1_extra_args,
    "--out", paste0(output_prefix, "_step1")
  )

  labels <- .gpt_step_labels(step2_bgen)
  step2 <- lapply(seq_along(step2_bgen), function(i) {
    args <- c(
      "--step", "2",
      "--bgen", step2_bgen[[i]],
      "--sample", step2_sample[[i]]
    )
    if (isTRUE(ref_first)) {
      args <- c(args, "--ref-first")
    }
    if (!is.null(step2_extract)) {
      args <- c(args, "--extract", step2_extract[[i]])
    }
    if (!is.null(resolved$keep_file)) {
      args <- c(args, "--keep", resolved$keep_file)
    }
    args <- c(
      args,
      common,
      trait_flag,
      "--pred", prediction_file,
      "--bsize", as.character(as.integer(bsize_step2))
    )
    if (identical(trait_type, "binary")) {
      if (isTRUE(firth)) {
        args <- c(args, "--firth", "--approx", "--pThresh", as.character(p_threshold))
      } else if (isTRUE(spa)) {
        args <- c(args, "--spa", "--pThresh", as.character(p_threshold))
      }
    }
    if (!is.null(threads)) {
      args <- c(args, "--threads", as.character(as.integer(threads)))
    }
    args <- c(
      args,
      step2_extra_args,
      "--out", paste0(output_prefix, "_step2_", labels[[i]])
    )
    .gpt_command("regenie", args)
  })
  names(step2) <- labels

  out <- list(
    tool = "REGENIE",
    phenotype_file = phenotype_file,
    covariate_file = resolved$covariate_file,
    phenotype_cols = phenotypes,
    covariates = covariates,
    categorical_covariates = categorical_covariates,
    trait_type = trait_type,
    step1 = .gpt_command("regenie", step1_command_args),
    step2 = step2,
    expected_prediction_file = prediction_file,
    settings = list(
      ref_first = ref_first,
      firth = firth,
      spa = spa,
      p_threshold = p_threshold,
      bsize_step1 = as.integer(bsize_step1),
      bsize_step2 = as.integer(bsize_step2),
      threads = threads,
      lowmem = lowmem,
      step1_args = step1_extra_args,
      step2_args = step2_extra_args
    )
  )
  class(out) <- c("ukb_regenie_plan", "list")
  out
}

#' Run a Planned REGENIE Analysis
#'
#' @description
#' Execute direct REGENIE commands from \code{\link{ukb_plan_regenie}} or
#' \code{\link{ukb_plan_regenie_command}}.
#' Execution is opt-in and requires a REGENIE binary. On RAP this is intended
#' for an existing Swiss Army Knife or controlled custom execution
#' environment; it does not itself submit a DNAnexus job.
#'
#' @param plan A plan inheriting from \code{ukb_regenie_plan}.
#' @param execute Logical. Defaults to \code{FALSE}; when false, the plan is
#'   returned without running anything.
#' @param steps Commands to run. For a structured plan, use
#'   \code{"step1"}, \code{"step2"}, or both. For a raw command plan, use its
#'   command label. \code{NULL} runs every command in the plan.
#' @param executable REGENIE executable path or command name.
#' @param require_rap Require a RAP-like environment before execution.
#' @param timeout Per-command timeout in seconds.
#'
#' @return The unchanged plan when \code{execute = FALSE}; otherwise a list of
#'   command statuses and captured logs.
#' @export
ukb_run_regenie <- function(plan,
                            execute = FALSE,
                            steps = c("step1", "step2"),
                            executable = "regenie",
                            require_rap = TRUE,
                            timeout = 86400) {
  if (!inherits(plan, "ukb_regenie_plan")) {
    stop("`plan` must be a `ukb_regenie_plan`.", call. = FALSE)
  }
  .gpp_assert_flag(execute, "execute")
  .gpp_assert_flag(require_rap, "require_rap")
  selected_steps <- if (missing(steps)) NULL else steps
  commands <- .gpt_regenie_commands(plan, selected_steps)
  if (!isTRUE(execute)) {
    return(plan)
  }
  if (isTRUE(require_rap)) {
    .ukb_assert_rap_env("ukb_run_regenie()")
  }
  executable_path <- unname(Sys.which(executable))
  if (!nzchar(executable_path)) {
    stop("REGENIE executable was not found on PATH.", call. = FALSE)
  }
  .gpt_positive_number(timeout, "timeout")

  results <- lapply(names(commands), function(label) {
    command <- commands[[label]]
    log <- system2(
      executable_path,
      args = vapply(command$args, shQuote, character(1)),
      stdout = TRUE,
      stderr = TRUE,
      timeout = timeout
    )
    status <- attr(log, "status")
    if (is.null(status)) {
      status <- 0L
    }
    list(
      label = label,
      status = as.integer(status),
      success = identical(as.integer(status), 0L),
      command = .gpt_display_command(executable, command$args),
      log = log
    )
  })
  names(results) <- names(commands)

  failed <- names(results)[!vapply(results, function(x) isTRUE(x$success), logical(1))]
  if (length(failed) > 0L) {
    warning(
      sprintf("REGENIE command(s) failed: %s", paste(failed, collapse = ", ")),
      call. = FALSE
    )
  }
  out <- list(tool = "REGENIE", results = results, plan = plan)
  class(out) <- c("ukb_regenie_run", "list")
  out
}

#' Plan an Arbitrary PLINK2 Command
#'
#' @description
#' Create an auditable, non-executing PLINK2 command plan from raw command-line
#' arguments. This is the escape hatch for PLINK2 operations that do not have
#' a dedicated UKBAnalytica wrapper. Arguments are passed as individual tokens
#' and quoted before execution.
#'
#' Dedicated wrappers such as \code{\link{convert_gwas_datatype}} and
#' \code{\link{ukb_plan_plink2_dosage}} remain preferable for their validated
#' inputs and output metadata.
#'
#' @param args Character vector containing PLINK2 command-line argument tokens
#'   in execution order.
#' @param expected_outputs Optional character vector of expected output paths.
#'   This is metadata only; PLINK2 remains responsible for creating them.
#' @param label Command label used in execution results.
#'
#' @return A list of class \code{ukb_plink2_plan}.
#' @export
ukb_plan_plink2 <- function(args,
                            expected_outputs = character(),
                            label = "plink2") {
  args <- .gpt_cli_args(args, "args")
  expected_outputs <- .gpt_cli_args(
    expected_outputs,
    "expected_outputs",
    allow_empty = TRUE
  )
  label <- .gpt_command_label(label)

  commands <- list(.gpt_command("plink2", args))
  names(commands) <- label
  out <- list(
    tool = "PLINK2",
    commands = commands,
    expected_outputs = expected_outputs,
    settings = list(mode = "raw_args")
  )
  class(out) <- c("ukb_plink2_plan", "list")
  out
}

#' Plan PLINK2 Dosage Extraction for PheWAS
#'
#' @description
#' Build auditable PLINK2 commands that extract selected variants from one or
#' more UKB BGEN files and export additive dosages. The function only creates a
#' plan; it does not execute PLINK2 or submit a RAP job.
#'
#' @param bgen One or more BGEN file paths.
#' @param sample One sample file path or one path per BGEN.
#' @param variant_file One variant inclusion file or one path per BGEN.
#' @param output_prefix Output prefix in RAP-controlled storage.
#' @param keep_file Optional two-column FID/IID sample inclusion file.
#' @param ref_first Pass the PLINK2 BGEN \code{ref-first} modifier.
#' @param plink2_args Additional PLINK2 argument tokens. Input, output,
#'   extraction, and sample-selection flags managed by this wrapper are
#'   rejected; use \code{\link{ukb_plan_plink2}} for fully custom commands.
#'
#' @return A list of class \code{ukb_plink2_dosage_plan}.
#' @export
ukb_plan_plink2_dosage <- function(bgen,
                                   sample,
                                   variant_file,
                                   output_prefix,
                                   keep_file = NULL,
                                   ref_first = TRUE,
                                   plink2_args = character()) {
  .gpp_assert_flag(ref_first, "ref_first")
  bgen <- .gpt_paths(bgen, "bgen")
  sample <- .gpt_recycle_paths(sample, length(bgen), "sample")
  variant_file <- .gpt_recycle_paths(
    variant_file,
    length(bgen),
    "variant_file"
  )
  output_prefix <- .gpt_path(output_prefix, "output_prefix")
  if (!is.null(keep_file)) {
    keep_file <- .gpt_path(keep_file, "keep_file")
  }
  plink2_args <- .gpt_cli_args(plink2_args, "plink2_args", allow_empty = TRUE)
  .gpt_reject_managed_flags(
    plink2_args,
    c(
      "--bgen", "--sample", "--extract", "--keep", "--export", "--out",
      "--bfile", "--pfile", "--vcf", "--bcf", "--make-bed", "--make-pgen"
    ),
    "ukb_plan_plink2_dosage()"
  )

  labels <- .gpt_step_labels(bgen)
  outputs <- if (length(bgen) == 1L) {
    output_prefix
  } else {
    paste0(output_prefix, "_", labels)
  }
  commands <- lapply(seq_along(bgen), function(i) {
    args <- c("--bgen", bgen[[i]])
    if (isTRUE(ref_first)) {
      args <- c(args, "ref-first")
    }
    args <- c(
      args,
      "--sample", sample[[i]],
      "--extract", variant_file[[i]]
    )
    if (!is.null(keep_file)) {
      args <- c(args, "--keep", keep_file)
    }
    args <- c(args, plink2_args, "--export", "A", "--out", outputs[[i]])
    .gpt_command("plink2", args)
  })
  names(commands) <- labels

  out <- list(
    tool = "PLINK2",
    commands = commands,
    expected_raw = setNames(paste0(outputs, ".raw"), labels),
    settings = list(
      bgen = bgen,
      sample = sample,
      variant_file = variant_file,
      keep_file = keep_file,
      ref_first = ref_first,
      plink2_args = plink2_args
    )
  )
  class(out) <- c("ukb_plink2_dosage_plan", "ukb_plink2_plan", "list")
  out
}

#' Run a Planned PLINK2 Dosage Extraction
#'
#' @description
#' Execute commands from \code{\link{ukb_plan_plink2_dosage}}. Execution is
#' opt-in, requires a PLINK2 binary, and does not itself submit a DNAnexus job.
#'
#' @param plan A plan inheriting from \code{ukb_plink2_plan}, including objects
#'   returned by \code{\link{ukb_plan_plink2}},
#'   \code{\link{ukb_plan_plink2_dosage}}, and
#'   \code{\link{convert_gwas_datatype}}.
#' @param execute Logical. Defaults to \code{FALSE}.
#' @param executable PLINK2 executable path or command name.
#' @param require_rap Require a RAP-like environment before execution.
#' @param timeout Per-command timeout in seconds.
#'
#' @return The unchanged plan when \code{execute = FALSE}; otherwise a list of
#'   command statuses and captured logs.
#' @export
ukb_run_plink2 <- function(plan,
                           execute = FALSE,
                           executable = "plink2",
                           require_rap = TRUE,
                           timeout = 86400) {
  if (!inherits(plan, "ukb_plink2_plan")) {
    stop("`plan` must inherit from `ukb_plink2_plan`.", call. = FALSE)
  }
  .gpp_assert_flag(execute, "execute")
  .gpp_assert_flag(require_rap, "require_rap")
  if (!isTRUE(execute)) {
    return(plan)
  }
  if (isTRUE(require_rap)) {
    .ukb_assert_rap_env("ukb_run_plink2()")
  }
  executable_path <- unname(Sys.which(executable))
  if (!nzchar(executable_path)) {
    stop("PLINK2 executable was not found on PATH.", call. = FALSE)
  }
  .gpt_positive_number(timeout, "timeout")

  results <- lapply(names(plan$commands), function(label) {
    command <- plan$commands[[label]]
    log <- system2(
      executable_path,
      args = vapply(command$args, shQuote, character(1)),
      stdout = TRUE,
      stderr = TRUE,
      timeout = timeout
    )
    status <- attr(log, "status")
    if (is.null(status)) {
      status <- 0L
    }
    list(
      label = label,
      status = as.integer(status),
      success = identical(as.integer(status), 0L),
      command = .gpt_display_command(executable, command$args),
      log = log
    )
  })
  names(results) <- names(plan$commands)

  failed <- names(results)[!vapply(results, function(x) isTRUE(x$success), logical(1))]
  if (length(failed) > 0L) {
    warning(
      sprintf("PLINK2 command(s) failed: %s", paste(failed, collapse = ", ")),
      call. = FALSE
    )
  }
  out <- list(tool = "PLINK2", results = results, plan = plan)
  class(out) <- c("ukb_plink2_run", "list")
  out
}

#' Run PheWAS Using the PheWAS R Package
#'
#' @description
#' Convert the long ICD table in a
#' \code{ukb_gwas_phewas_phenotype} object to phecodes and run the external
#' PheWAS R package. The dependency is optional and is never installed
#' automatically.
#'
#' @param x A \code{ukb_gwas_phewas_phenotype}.
#' @param genotypes Data frame containing participant IDs and one or more
#'   additive genotype/dosage columns.
#' @param genotype_id_col Participant ID column in \code{genotypes}.
#' @param genotype_cols Genotype columns. Defaults to every non-ID column.
#' @param covariate_cols Covariates from \code{x} used in every model.
#' @param min_code_count Minimum code count used by
#'   \code{PheWAS::createPhenotypes}.
#' @param min_records Minimum case/control records used by
#'   \code{PheWAS::phewas}.
#' @param cores Number of R worker processes.
#' @param additive_genotypes Treat genotype inputs as additive dosages.
#' @param add_phecode_exclusions Apply phecode control exclusions.
#' @param sex_restrictions Apply sex-specific phenotype restrictions when a
#'   sex table was prepared in \code{x}.
#' @param vocabulary_map PheWAS vocabulary map: \code{"icd10"} matches the
#'   official RAP tutorial, while \code{"default"} supports the package's
#'   combined map.
#' @param significance_threshold Threshold columns requested from PheWAS.
#' @param alpha Family-wise/FDR alpha.
#'
#' @return A list of class \code{ukb_phewas_result}.
#' @export
ukb_run_phewas <- function(
    x,
    genotypes,
    genotype_id_col = "IID",
    genotype_cols = NULL,
    covariate_cols = x$settings$covariates,
    min_code_count = 1L,
    min_records = 20L,
    cores = 1L,
    additive_genotypes = TRUE,
    add_phecode_exclusions = TRUE,
    sex_restrictions = TRUE,
    vocabulary_map = c("icd10", "default"),
    significance_threshold = c("p-value", "bonferroni", "fdr"),
    alpha = 0.05) {
  if (!inherits(x, "ukb_gwas_phewas_phenotype")) {
    stop("`x` must be a `ukb_gwas_phewas_phenotype` object.", call. = FALSE)
  }
  if (!requireNamespace("PheWAS", quietly = TRUE)) {
    stop(
      paste(
        "The optional PheWAS package is required.",
        "Install the reviewed upstream package with",
        "`remotes::install_github(\"PheWAS/PheWAS\")` inside RAP."
      ),
      call. = FALSE
    )
  }
  .gpp_assert_flag(additive_genotypes, "additive_genotypes")
  .gpp_assert_flag(add_phecode_exclusions, "add_phecode_exclusions")
  .gpp_assert_flag(sex_restrictions, "sex_restrictions")
  .gpt_positive_integer(min_code_count, "min_code_count")
  .gpt_positive_integer(min_records, "min_records")
  .gpt_positive_integer(cores, "cores")
  .gpt_positive_number(alpha, "alpha")
  if (alpha >= 1) {
    stop("`alpha` must be less than 1.", call. = FALSE)
  }
  vocabulary_map <- match.arg(vocabulary_map)
  significance_threshold <- match.arg(
    significance_threshold,
    c("p-value", "bonferroni", "fdr"),
    several.ok = TRUE
  )

  if (nrow(x$phewas_long) == 0L) {
    stop("The phenotype object contains no PheWAS diagnosis codes.", call. = FALSE)
  }
  genotype <- as.data.table(copy(genotypes))
  .gpp_require_columns(genotype, genotype_id_col, "genotypes")
  if (anyDuplicated(genotype[[genotype_id_col]])) {
    stop("`genotypes` must contain one row per participant ID.", call. = FALSE)
  }
  if (is.null(genotype_cols)) {
    genotype_cols <- setdiff(names(genotype), genotype_id_col)
  }
  genotype_cols <- .gpt_names(genotype_cols, "genotype_cols")
  .gpp_require_columns(genotype, genotype_cols, "genotypes")
  if (isTRUE(additive_genotypes)) {
    invalid <- vapply(genotype_cols, function(column) {
      values <- genotype[[column]]
      observed <- values[!is.na(values)]
      !is.numeric(values) || any(observed < 0 | observed > 2)
    }, logical(1))
    if (any(invalid)) {
      stop(
        sprintf(
          "Additive genotype columns must be numeric dosages in [0,2]: %s",
          paste(genotype_cols[invalid], collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  genotype <- genotype[, c(genotype_id_col, genotype_cols), with = FALSE]
  setnames(genotype, genotype_id_col, "id")
  set(genotype, j = "id", value = as.character(genotype[["id"]]))
  genotype <- as.data.frame(genotype)

  covariate_cols <- .gpt_names(covariate_cols, "covariate_cols", allow_empty = TRUE)
  if (!all(covariate_cols %in% x$settings$covariates)) {
    stop("`covariate_cols` must be present in the phenotype object.", call. = FALSE)
  }
  covariate <- if (length(covariate_cols) > 0L) {
    value <- x$covariates[, c("IID", covariate_cols), drop = FALSE]
    names(value)[[1L]] <- "id"
    value$id <- as.character(value$id)
    factor_cols <- intersect(
      covariate_cols,
      x$settings$categorical_covariates
    )
    for (column in factor_cols) {
      value[[column]] <- factor(value[[column]])
    }
    value
  } else {
    c(NA)
  }

  code_input <- x$phewas_long[, c("id", "vocabulary_id", "code", "count"), drop = FALSE]
  code_input$id <- as.character(code_input$id)
  full_population_ids <- as.character(x$sample_ids$IID)
  shared_ids <- intersect(full_population_ids, genotype$id)
  if (length(shared_ids) == 0L) {
    stop(
      "There are no shared participant IDs between the phenotype and genotype tables.",
      call. = FALSE
    )
  }
  create_args <- list(
    id.vocab.code.index = code_input,
    min.code.count = as.integer(min_code_count),
    add.phecode.exclusions = add_phecode_exclusions,
    translate = TRUE,
    full.population.ids = full_population_ids,
    aggregate.fun = sum,
    vocabulary.map = if (identical(vocabulary_map, "icd10")) {
      getExportedValue("PheWAS", "phecode_map_icd10")
    } else {
      getExportedValue("PheWAS", "phecode_map")
    }
  )
  if (isTRUE(sex_restrictions) && !is.null(x$sex)) {
    id_sex <- x$sex
    id_sex$id <- as.character(id_sex$id)
    create_args$id.sex <- id_sex
  }
  phenotype_matrix <- do.call(
    getExportedValue("PheWAS", "createPhenotypes"),
    create_args
  )

  phewas_args <- list(
    phenotypes = phenotype_matrix,
    genotypes = genotype,
    covariates = covariate,
    cores = as.integer(cores),
    additive.genotypes = additive_genotypes,
    significance.threshold = significance_threshold,
    alpha = alpha,
    min.records = as.integer(min_records)
  )
  results <- do.call(getExportedValue("PheWAS", "phewas"), phewas_args)
  results <- .gpt_add_phecode_info(results)

  out <- list(
    results = results,
    n_participants = length(shared_ids),
    n_phenotype_participants = nrow(phenotype_matrix),
    n_phenotypes = max(0L, ncol(phenotype_matrix) - 1L),
    genotype_cols = genotype_cols,
    covariate_cols = covariate_cols,
    settings = list(
      min_code_count = as.integer(min_code_count),
      min_records = as.integer(min_records),
      cores = as.integer(cores),
      additive_genotypes = additive_genotypes,
      add_phecode_exclusions = add_phecode_exclusions,
      sex_restrictions = sex_restrictions,
      vocabulary_map = vocabulary_map,
      significance_threshold = significance_threshold,
      alpha = alpha
    )
  )
  class(out) <- c("ukb_phewas_result", "list")
  out
}

#' @export
print.ukb_regenie_command_plan <- function(x, ...) {
  cat("UKB REGENIE command plan\n")
  cat("  Commands: ", length(x$commands), "\n", sep = "")
  if (length(x$commands) > 0L) {
    cat("  First command: ", x$commands[[1L]]$display, "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.ukb_regenie_plan <- function(x, ...) {
  cat("UKB REGENIE plan\n")
  cat("  Trait type: ", x$trait_type, "\n", sep = "")
  cat("  Phenotypes: ", paste(x$phenotype_cols, collapse = ", "), "\n", sep = "")
  cat("  Step 2 jobs: ", length(x$step2), "\n", sep = "")
  cat("  Step 1: ", x$step1$display, "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_plink2_dosage_plan <- function(x, ...) {
  cat("UKB PLINK2 dosage-extraction plan\n")
  cat("  Commands: ", length(x$commands), "\n", sep = "")
  cat("  Expected .raw files: ", length(x$expected_raw), "\n", sep = "")
  invisible(x)
}

#' @export
print.ukb_plink2_plan <- function(x, ...) {
  cat("UKB PLINK2 command plan\n")
  cat("  Commands: ", length(x$commands), "\n", sep = "")
  if (length(x$commands) > 0L) {
    cat("  First command: ", x$commands[[1L]]$display, "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.ukb_phewas_result <- function(x, ...) {
  cat("UKB PheWAS result\n")
  cat("  Participants: ", x$n_participants, "\n", sep = "")
  cat("  Phenotypes tested: ", x$n_phenotypes, "\n", sep = "")
  cat("  Genotype predictors: ", length(x$genotype_cols), "\n", sep = "")
  invisible(x)
}

.gpt_resolve_files <- function(files,
                               phenotypes,
                               covariate_file,
                               covariates,
                               categorical_covariates,
                               trait_type,
                               keep_file) {
  if (inherits(files, "ukb_gwas_phewas_files")) {
    settings <- files$settings
    selected <- if (is.null(phenotypes)) settings$phenotype_cols else phenotypes
    available_types <- settings$trait_types[selected]
    if (is.null(trait_type)) {
      if (anyNA(available_types) || length(unique(unname(available_types))) != 1L) {
        stop("Selected traits must share one type, or provide `trait_type`.", call. = FALSE)
      }
      trait_type <- unique(unname(available_types))
    }
    return(list(
      phenotype_file = files$paths[["gwas"]],
      covariate_file = if (is.null(covariate_file)) files$paths[["covariates"]] else covariate_file,
      phenotypes = selected,
      covariates = if (is.null(covariates)) settings$covariates else covariates,
      categorical_covariates = if (is.null(categorical_covariates)) {
        settings$categorical_covariates
      } else {
        categorical_covariates
      },
      trait_type = trait_type,
      keep_file = if (is.null(keep_file)) files$paths[["sample_keep"]] else keep_file
    ))
  }
  list(
    phenotype_file = .gpt_path(files, "files"),
    covariate_file = if (is.null(covariate_file)) .gpt_path(files, "files") else .gpt_path(covariate_file, "covariate_file"),
    phenotypes = phenotypes,
    covariates = if (is.null(covariates)) character() else covariates,
    categorical_covariates = if (is.null(categorical_covariates)) character() else categorical_covariates,
    trait_type = trait_type,
    keep_file = if (is.null(keep_file)) NULL else .gpt_path(keep_file, "keep_file")
  )
}

.gpt_command <- function(executable, args) {
  list(
    executable = executable,
    args = args,
    display = .gpt_display_command(executable, args)
  )
}

.gpt_display_command <- function(executable, args) {
  paste(c(executable, vapply(args, shQuote, character(1))), collapse = " ")
}

.gpt_cli_args <- function(x, name, allow_empty = FALSE) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) ||
      (!isTRUE(allow_empty) && length(x) == 0L)) {
    stop(
      sprintf(
        "`%s` must be %scharacter vector of non-empty command tokens.",
        name,
        if (isTRUE(allow_empty)) "a " else "a non-empty "
      ),
      call. = FALSE
    )
  }
  x
}

.gpt_command_label <- function(label) {
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("`label` must be one non-empty string.", call. = FALSE)
  }
  label
}

.gpt_reject_managed_flags <- function(args,
                                      forbidden,
                                      wrapper,
                                      tool = "PLINK2",
                                      generic = "ukb_plan_plink2()") {
  if (length(args) == 0L) {
    return(invisible(TRUE))
  }
  flags <- sub("=.*$", "", args)
  conflicts <- unique(flags[flags %in% forbidden])
  if (length(conflicts) > 0L) {
    stop(
      sprintf(
        "%s manages these %s flag(s): %s. Use `%s` for a fully custom command.",
        wrapper,
        tool,
        paste(conflicts, collapse = ", "),
        generic
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.gpt_regenie_commands <- function(plan, steps) {
  if (inherits(plan, "ukb_regenie_command_plan")) {
    commands <- plan$commands
    if (is.null(steps)) {
      return(commands)
    }
    steps <- unique(.gpt_cli_args(steps, "steps"))
    unknown <- setdiff(steps, names(commands))
    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "Unknown REGENIE command label(s): %s.",
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(commands[steps])
  }

  steps <- if (is.null(steps)) {
    c("step1", "step2")
  } else {
    unique(match.arg(steps, c("step1", "step2"), several.ok = TRUE))
  }
  commands <- list()
  if ("step1" %in% steps) {
    commands$step1 <- plan$step1
  }
  if ("step2" %in% steps) {
    commands <- c(commands, plan$step2)
  }
  commands
}

.gpt_add_phecode_info <- function(results) {
  info <- tryCatch(
    getExportedValue("PheWAS", "pheinfo"),
    error = function(e) NULL
  )
  if (is.null(info)) {
    warning(
      "PheWAS completed, but upstream phecode annotations were unavailable.",
      call. = FALSE
    )
    return(results)
  }
  candidates <- grep(
    "pheno|phewas|phecode",
    names(results),
    ignore.case = TRUE,
    value = TRUE
  )
  phenotype_col <- if (length(candidates) > 0L) candidates[[1L]] else names(results)[[1L]]
  index <- match(as.character(results[[phenotype_col]]), as.character(info[["phecode"]]))
  for (column in intersect(c("description", "group"), names(info))) {
    if (!column %in% names(results)) {
      results[[column]] <- info[[column]][index]
    }
  }
  results
}

.gpt_names <- function(x, name, allow_empty = FALSE) {
  if (is.null(x)) {
    if (isTRUE(allow_empty)) {
      return(character())
    }
    stop(sprintf("`%s` is required.", name), call. = FALSE)
  }
  .gpp_character_vector(as.character(x), name, allow_empty = allow_empty)
}

.gpt_path <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be one non-empty path.", name), call. = FALSE)
  }
  x
}

.gpt_paths <- function(x, name) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) || any(!nzchar(x))) {
    stop(sprintf("`%s` must contain one or more non-empty paths.", name), call. = FALSE)
  }
  x
}

.gpt_recycle_paths <- function(x, n, name) {
  x <- .gpt_paths(x, name)
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(sprintf("`%s` must have length 1 or match `step2_bgen`.", name), call. = FALSE)
  }
  x
}

.gpt_step_labels <- function(paths) {
  base <- basename(paths)
  labels <- sub(
    ".*(?:chr|_c)([0-9]+|X|Y|XY|MT).*",
    "chr\\1",
    base,
    ignore.case = TRUE,
    perl = TRUE
  )
  unchanged <- labels == base
  labels[unchanged] <- paste0("part", seq_along(paths)[unchanged])
  make.unique(labels, sep = "_")
}

.gpt_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.gpt_positive_number <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(
      sprintf("`%s` must be one %s number.", name, if (allow_zero) "non-negative" else "positive"),
      call. = FALSE
    )
  }
  lower_ok <- if (isTRUE(allow_zero)) x >= 0 else x > 0
  if (!lower_ok) {
    stop(
      sprintf("`%s` must be one %s number.", name, if (allow_zero) "non-negative" else "positive"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
