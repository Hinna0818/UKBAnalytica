#' Convert Common GWAS Genotype File Formats
#'
#' @description
#' Build or execute a PLINK2 command that converts between common GWAS
#' genotype formats. Supported high-level formats are PLINK 1 binary
#' (\code{bed/bim/fam}), PLINK 2 binary (\code{pgen/pvar/psam}), BGEN,
#' VCF, and BCF.
#'
#' The function is a dry run by default. It does not install PLINK2, submit a
#' DNAnexus job, or copy participant data into the package. With
#' \code{execute = TRUE}, it calls an existing PLINK2 binary through
#' \code{\link{ukb_run_plink2}}.
#'
#' @param input Input path. For BED and PGEN filesets this may be the fileset
#'   prefix or any primary component path; for BGEN, VCF, and BCF it is the
#'   data-file path.
#' @param to Output format: \code{"bed"}, \code{"pgen"}, \code{"bgen"},
#'   \code{"vcf"}, or \code{"bcf"}.
#' @param output_prefix Output prefix without a genotype-format extension.
#' @param from Input format or \code{"auto"}. Automatic detection uses the
#'   supplied extension and, for a prefix, companion files visible in the
#'   current filesystem.
#' @param sample_file Optional BGEN \code{.sample} file or VCF/BCF
#'   \code{.psam} file.
#' @param bgen_ref Reference-allele interpretation for BGEN input.
#' @param vcf_dosage Optional VCF/BCF dosage field name, for example
#'   \code{"DS"}.
#' @param bgen_version BGEN export version.
#' @param bgen_bits Integer probability precision for BGEN export. PLINK2
#'   defaults to 16 bits.
#' @param export_ref_first Put REF first when exporting BGEN.
#' @param vcf_bgz Block-gzip VCF output.
#' @param vcf_dosage_export Optional VCF/BCF dosage export mode: \code{"GP"},
#'   \code{"DS"}, \code{"DS-force"}, \code{"DS-only"}, \code{"HDS"}, or
#'   \code{"HDS-force"}.
#' @param plink2_args Additional PLINK2 argument tokens, such as
#'   \code{c("--maf", "0.01", "--geno", "0.05", "--threads", "16")}.
#'   Input/output flags managed by this function are rejected. Use
#'   \code{\link{ukb_plan_plink2}} for a fully custom command.
#' @param execute Execute the conversion. Defaults to \code{FALSE}.
#' @param executable PLINK2 executable path or command name.
#' @param require_rap Require a RAP-like environment before execution.
#' @param overwrite Permit PLINK2 to replace an existing expected output.
#' @param timeout Command timeout in seconds.
#'
#' @return With \code{execute = FALSE}, a
#'   \code{ukb_gwas_conversion_plan}; otherwise a \code{ukb_plink2_run}.
#' @export
convert_gwas_datatype <- function(
    input,
    to,
    output_prefix,
    from = c("auto", "bed", "pgen", "bgen", "vcf", "bcf"),
    sample_file = NULL,
    bgen_ref = c("ref-first", "ref-last", "ref-unknown"),
    vcf_dosage = NULL,
    bgen_version = c("1.2", "1.3", "1.1"),
    bgen_bits = 16L,
    export_ref_first = TRUE,
    vcf_bgz = TRUE,
    vcf_dosage_export = NULL,
    plink2_args = character(),
    execute = FALSE,
    executable = "plink2",
    require_rap = TRUE,
    overwrite = FALSE,
    timeout = 86400) {
  input <- .gpt_path(input, "input")
  output_prefix <- .gpt_path(output_prefix, "output_prefix")
  from <- match.arg(from)
  to <- match.arg(to, c("bed", "pgen", "bgen", "vcf", "bcf"))
  bgen_ref <- match.arg(bgen_ref)
  bgen_version <- match.arg(bgen_version)
  .gpp_assert_flag(export_ref_first, "export_ref_first")
  .gpp_assert_flag(vcf_bgz, "vcf_bgz")
  .gpp_assert_flag(execute, "execute")
  .gpp_assert_flag(require_rap, "require_rap")
  .gpp_assert_flag(overwrite, "overwrite")
  .gpt_positive_integer(bgen_bits, "bgen_bits")
  if (bgen_bits > 32L) {
    stop("`bgen_bits` must be between 1 and 32.", call. = FALSE)
  }
  .gpt_positive_number(timeout, "timeout")

  if (!is.null(sample_file)) {
    sample_file <- .gpt_path(sample_file, "sample_file")
  }
  if (!is.null(vcf_dosage)) {
    vcf_dosage <- .gcd_single_token(vcf_dosage, "vcf_dosage")
  }
  if (!is.null(vcf_dosage_export)) {
    vcf_dosage_export <- match.arg(
      vcf_dosage_export,
      c("GP", "DS", "DS-force", "DS-only", "HDS", "HDS-force")
    )
  }
  plink2_args <- .gpt_cli_args(plink2_args, "plink2_args", allow_empty = TRUE)
  .gpt_reject_managed_flags(
    plink2_args,
    c(
      "--bfile", "--bed", "--bim", "--fam",
      "--pfile", "--pgen", "--pvar", "--psam",
      "--bgen", "--sample", "--vcf", "--bcf",
      "--make-bed", "--make-pgen", "--make-bpgen", "--export", "--out"
    ),
    "convert_gwas_datatype()"
  )

  if (.gcd_has_format_extension(output_prefix)) {
    stop(
      "`output_prefix` must omit the genotype-format extension.",
      call. = FALSE
    )
  }
  from <- .gcd_resolve_format(input, from)
  normalized_input <- .gcd_normalize_input(input, from)
  pgen_vzs <- identical(from, "pgen") && (
    grepl("\\.pvar\\.zst$", input, ignore.case = TRUE) ||
      (!file.exists(paste0(normalized_input, ".pvar")) &&
         file.exists(paste0(normalized_input, ".pvar.zst")))
  )
  source_prefix <- .gcd_source_prefix(normalized_input, from)
  if (identical(
    normalizePath(source_prefix, winslash = "/", mustWork = FALSE),
    normalizePath(path.expand(output_prefix), winslash = "/", mustWork = FALSE)
  )) {
    stop(
      "`output_prefix` must differ from the input prefix to protect source files.",
      call. = FALSE
    )
  }

  if (!is.null(sample_file) && from %in% c("bed", "pgen")) {
    stop(
      "`sample_file` is only used with BGEN, VCF, or BCF input.",
      call. = FALSE
    )
  }
  if (!is.null(vcf_dosage) && !from %in% c("vcf", "bcf")) {
    stop("`vcf_dosage` is only valid for VCF/BCF input.", call. = FALSE)
  }
  if (!is.null(vcf_dosage_export) && !to %in% c("vcf", "bcf")) {
    stop(
      "`vcf_dosage_export` is only valid for VCF/BCF output.",
      call. = FALSE
    )
  }

  input_args <- .gcd_input_args(
    input = normalized_input,
    from = from,
    sample_file = sample_file,
    bgen_ref = bgen_ref,
    vcf_dosage = vcf_dosage,
    pgen_vzs = pgen_vzs
  )
  output_args <- .gcd_output_args(
    to = to,
    bgen_version = bgen_version,
    bgen_bits = as.integer(bgen_bits),
    export_ref_first = export_ref_first,
    vcf_bgz = vcf_bgz,
    vcf_dosage_export = vcf_dosage_export
  )
  expected_outputs <- .gcd_expected_outputs(
    output_prefix = output_prefix,
    to = to,
    vcf_bgz = vcf_bgz
  )
  args <- c(input_args, plink2_args, output_args, "--out", output_prefix)

  command <- .gpt_command("plink2", args)
  out <- list(
    tool = "PLINK2",
    commands = list(convert = command),
    expected_outputs = expected_outputs,
    input_format = from,
    output_format = to,
    settings = list(
      input = normalized_input,
      output_prefix = output_prefix,
      sample_file = sample_file,
      pgen_vzs = pgen_vzs,
      bgen_ref = bgen_ref,
      vcf_dosage = vcf_dosage,
      bgen_version = bgen_version,
      bgen_bits = as.integer(bgen_bits),
      export_ref_first = export_ref_first,
      vcf_bgz = vcf_bgz,
      vcf_dosage_export = vcf_dosage_export,
      plink2_args = plink2_args
    )
  )
  class(out) <- c("ukb_gwas_conversion_plan", "ukb_plink2_plan", "list")

  if (!isTRUE(execute)) {
    return(out)
  }
  .gcd_validate_execution_inputs(out)
  existing <- expected_outputs[file.exists(expected_outputs)]
  if (length(existing) > 0L && !isTRUE(overwrite)) {
    stop(
      sprintf(
        "Expected output file(s) already exist: %s",
        paste(basename(existing), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  output_dir <- dirname(path.expand(output_prefix))
  if (!dir.exists(output_dir)) {
    stop("The parent directory of `output_prefix` does not exist.", call. = FALSE)
  }
  ukb_run_plink2(
    out,
    execute = TRUE,
    executable = executable,
    require_rap = require_rap,
    timeout = timeout
  )
}

#' @export
print.ukb_gwas_conversion_plan <- function(x, ...) {
  cat("UKB GWAS datatype conversion plan\n")
  cat("  Conversion: ", x$input_format, " -> ", x$output_format, "\n", sep = "")
  cat("  Command: ", x$commands$convert$display, "\n", sep = "")
  cat("  Expected files: ", paste(x$expected_outputs, collapse = ", "), "\n", sep = "")
  invisible(x)
}

.gcd_single_token <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      grepl("[[:space:]]", x)) {
    stop(sprintf("`%s` must be one non-empty token.", name), call. = FALSE)
  }
  x
}

.gcd_has_format_extension <- function(path) {
  grepl(
    "\\.(bed|bim|fam|pgen|pvar|pvar\\.zst|psam|bgen|vcf|vcf\\.gz|vcf\\.bgz|bcf)$",
    path,
    ignore.case = TRUE
  )
}

.gcd_resolve_format <- function(input, from) {
  if (!identical(from, "auto")) {
    return(from)
  }
  lower <- tolower(input)
  if (grepl("\\.(bed|bim|fam)$", lower)) {
    return("bed")
  }
  if (grepl("\\.(pgen|pvar|pvar\\.zst|psam)$", lower)) {
    return("pgen")
  }
  if (grepl("\\.bgen$", lower)) {
    return("bgen")
  }
  if (grepl("\\.vcf(\\.gz|\\.bgz)?$", lower)) {
    return("vcf")
  }
  if (grepl("\\.bcf$", lower)) {
    return("bcf")
  }

  expanded <- path.expand(input)
  candidates <- c(
    bed = all(file.exists(paste0(expanded, c(".bed", ".bim", ".fam")))),
    pgen = file.exists(paste0(expanded, ".pgen")) &&
      (file.exists(paste0(expanded, ".pvar")) ||
         file.exists(paste0(expanded, ".pvar.zst"))) &&
      file.exists(paste0(expanded, ".psam"))
  )
  detected <- names(candidates)[candidates]
  if (length(detected) == 1L) {
    return(detected)
  }
  if (length(detected) > 1L) {
    stop(
      "Both BED and PGEN companion files were detected; specify `from`.",
      call. = FALSE
    )
  }
  stop(
    "Could not infer the input format; specify `from` explicitly.",
    call. = FALSE
  )
}

.gcd_normalize_input <- function(input, from) {
  expanded <- path.expand(input)
  if (identical(from, "bed")) {
    return(sub("\\.(bed|bim|fam)$", "", expanded, ignore.case = TRUE))
  }
  if (identical(from, "pgen")) {
    return(sub("\\.(pgen|pvar|pvar\\.zst|psam)$", "", expanded, ignore.case = TRUE))
  }
  expanded
}

.gcd_source_prefix <- function(input, from) {
  if (from %in% c("bed", "pgen")) {
    return(input)
  }
  sub(
    "\\.(bgen|vcf|vcf\\.gz|vcf\\.bgz|bcf)$",
    "",
    input,
    ignore.case = TRUE
  )
}

.gcd_input_args <- function(input,
                            from,
                            sample_file,
                            bgen_ref,
                            vcf_dosage,
                            pgen_vzs) {
  if (identical(from, "bed")) {
    return(c("--bfile", input))
  }
  if (identical(from, "pgen")) {
    args <- c("--pfile", input)
    if (isTRUE(pgen_vzs)) {
      args <- c(args, "vzs")
    }
    return(args)
  }
  if (identical(from, "bgen")) {
    args <- c("--bgen", input, bgen_ref)
    if (!is.null(sample_file)) {
      args <- c(args, "--sample", sample_file)
    }
    return(args)
  }
  args <- c(if (identical(from, "vcf")) "--vcf" else "--bcf", input)
  if (!is.null(vcf_dosage)) {
    args <- c(args, paste0("dosage=", vcf_dosage))
  }
  if (!is.null(sample_file)) {
    args <- c(args, "--psam", sample_file)
  }
  args
}

.gcd_output_args <- function(to,
                             bgen_version,
                             bgen_bits,
                             export_ref_first,
                             vcf_bgz,
                             vcf_dosage_export) {
  if (identical(to, "bed")) {
    return("--make-bed")
  }
  if (identical(to, "pgen")) {
    return("--make-pgen")
  }
  if (identical(to, "bgen")) {
    args <- c(
      "--export",
      paste0("bgen-", bgen_version),
      paste0("bits=", bgen_bits)
    )
    if (isTRUE(export_ref_first)) {
      args <- c(args, "ref-first")
    }
    return(args)
  }

  args <- c("--export", to)
  if (identical(to, "vcf") && isTRUE(vcf_bgz)) {
    args <- c(args, "bgz")
  }
  if (!is.null(vcf_dosage_export)) {
    args <- c(args, paste0("vcf-dosage=", vcf_dosage_export))
  }
  args
}

.gcd_expected_outputs <- function(output_prefix, to, vcf_bgz) {
  switch(
    to,
    bed = paste0(output_prefix, c(".bed", ".bim", ".fam")),
    pgen = paste0(output_prefix, c(".pgen", ".pvar", ".psam")),
    bgen = paste0(output_prefix, c(".bgen", ".sample")),
    vcf = paste0(output_prefix, if (isTRUE(vcf_bgz)) ".vcf.gz" else ".vcf"),
    bcf = paste0(output_prefix, ".bcf")
  )
}

.gcd_validate_execution_inputs <- function(plan) {
  input <- plan$settings$input
  from <- plan$input_format
  paths <- if (identical(from, "bed")) {
    paste0(input, c(".bed", ".bim", ".fam"))
  } else if (identical(from, "pgen")) {
    pvar <- if (file.exists(paste0(input, ".pvar"))) {
      paste0(input, ".pvar")
    } else {
      paste0(input, ".pvar.zst")
    }
    c(paste0(input, ".pgen"), pvar, paste0(input, ".psam"))
  } else {
    input
  }
  if (!is.null(plan$settings$sample_file)) {
    paths <- c(paths, plan$settings$sample_file)
  }
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "Input file(s) do not exist: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
