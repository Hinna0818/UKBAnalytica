#' @title Create Disease Definition Object
#'
#' @description
#' Helper function to create a standardized disease definition object
#' containing ICD-10/ICD-9 patterns and self-report codes.
#'
#' @param name Full disease name (e.g., "Aortic Aneurysm").
#' @param icd10_pattern Regular expression pattern for ICD-10 codes (optional).
#' @param icd9_pattern Regular expression pattern for ICD-9 codes (optional).
#' @param sr_codes Integer vector of UKB self-report illness codes (optional).
#'
#' @return A list containing the disease definition parameters.
#'
#' @examples
#' \dontrun{
#' aa_def <- create_disease_definition(
#'   name = "Aortic Aneurysm",
#'   icd10_pattern = "^I71",
#'   icd9_pattern = "^441"
#' )
#' }
#'
#' @export
create_disease_definition <- function(name,
                                       icd10_pattern = NULL,
                                       icd9_pattern = NULL,
                                       sr_codes = NULL) {
  list(
    name = name,
    icd10_pattern = icd10_pattern,
    icd9_pattern = icd9_pattern,
    sr_codes = sr_codes
  )
}


#' @title Get Predefined Disease Definitions
#'
#' @description
#' Returns a list of commonly used cardiovascular and metabolic disease
#' definitions with validated ICD-10, ICD-9, and self-report code mappings.
#'
#' @return A named list of disease definition objects.
#'
#' @details
#' Included diseases:
#' \describe{
#'   \item{AA}{Aortic Aneurysm (I71, 441)}
#'   \item{TAA}{Thoracic Aortic Aneurysm}
#'   \item{AAA}{Abdominal Aortic Aneurysm}
#'   \item{CVD}{Cardiovascular Disease}
#'   \item{MI}{Myocardial Infarction}
#'   \item{HF}{Heart Failure}
#'   \item{Stroke}{Stroke (ischemic and hemorrhagic)}
#'   \item{Hypertension}{Essential and secondary hypertension}
#'   \item{Diabetes}{Diabetes Mellitus (all types)}
#'   \item{Vascular_Disease}{Peripheral vascular disease}
#' }
#'
#' @examples
#' \dontrun{
#' diseases <- get_predefined_diseases()
#' my_diseases <- diseases[c("AA", "Hypertension", "Diabetes")]
#' }
#'
#' @export
get_predefined_diseases <- function() {
  list(
    # Aortic diseases
    AA = create_disease_definition(
      name = "Aortic Aneurysm",
      icd10_pattern = "^I71",
      icd9_pattern = "^441"
    ),
    TAA = create_disease_definition(
      name = "Thoracic Aortic Aneurysm",
      icd10_pattern = "^(I710|I711|I712|I715|I716)",
      icd9_pattern = "^(4410|4411|4412)"
    ),
    AAA = create_disease_definition(
      name = "Abdominal Aortic Aneurysm",
      icd10_pattern = "^(I713|I714)",
      icd9_pattern = "^(4413|4414)"
    ),

    # Cardiovascular diseases
    CVD = create_disease_definition(
      name = "Cardiovascular Disease",
      icd10_pattern = "^(I21|I22|I23|I24|I25)",
      icd9_pattern = "^(410|411|412|413|414)",
      sr_codes = c(1066, 1067)
    ),
    MI = create_disease_definition(
      name = "Myocardial Infarction",
      icd10_pattern = "^(I21|I22)",
      icd9_pattern = "^410",
      sr_codes = c(1066)
    ),
    HF = create_disease_definition(
      name = "Heart Failure",
      icd10_pattern = "^(I50|I420|I426|I427|I429|I110)",
      icd9_pattern = "^(428|4254)",
      sr_codes = c(1076)
    ),
    Stroke = create_disease_definition(
      name = "Stroke",
      icd10_pattern = "^(I60|I61|I62|I63|I64)",
      icd9_pattern = "^(430|431|432|433|434|436)",
      sr_codes = c(1068)
    ),

    # Metabolic diseases
    Hypertension = create_disease_definition(
      name = "Hypertension",
      icd10_pattern = "^(I10|I11|I12|I13|I14|I15)",
      icd9_pattern = "^(401|402|403|404|405)",
      sr_codes = c(1065)
    ),
    Diabetes = create_disease_definition(
      name = "Diabetes Mellitus",
      icd10_pattern = "^(E10|E11|E12|E13|E14)",
      icd9_pattern = "^250",
      sr_codes = c(1220, 1221, 1222, 1223)
    ),

    # Vascular diseases
    Vascular_Disease = create_disease_definition(
      name = "Vascular Disease",
      icd10_pattern = "^(I71|I72|I73|I77|I78|I79)",
      icd9_pattern = "^(441|442|443|447)"
    )
  )
}


#' @title Combine Multiple Disease Definitions
#'
#' @description
#' Merges multiple disease definitions into a single composite endpoint definition.
#' Useful for creating MACE (Major Adverse Cardiovascular Events) or similar
#' composite outcomes.
#'
#' @param ... Disease definition objects to combine.
#' @param name Name for the composite outcome.
#'
#' @return A combined disease definition object.
#'
#' @examples
#' \dontrun{
#' diseases <- get_predefined_diseases()
#' mace <- combine_disease_definitions(
#'   diseases$MI, diseases$Stroke, diseases$HF,
#'   name = "MACE"
#' )
#' }
#'
#' @export
combine_disease_definitions <- function(..., name = "Combined") {
  defs <- list(...)

  # Combine ICD-10 patterns
  icd10_patterns <- sapply(defs, function(x) x$icd10_pattern)
  icd10_patterns <- icd10_patterns[!sapply(icd10_patterns, is.null)]
  icd10_combined <- if (length(icd10_patterns) > 0) {
    paste0("(", paste(icd10_patterns, collapse = "|"), ")")
  } else NULL

  # Combine ICD-9 patterns
  icd9_patterns <- sapply(defs, function(x) x$icd9_pattern)
  icd9_patterns <- icd9_patterns[!sapply(icd9_patterns, is.null)]
  icd9_combined <- if (length(icd9_patterns) > 0) {
    paste0("(", paste(icd9_patterns, collapse = "|"), ")")
  } else NULL

  # Combine self-report codes
  sr_codes <- unlist(lapply(defs, function(x) x$sr_codes))
  sr_codes <- unique(sr_codes[!is.na(sr_codes)])
  if (length(sr_codes) == 0) sr_codes <- NULL

  create_disease_definition(
    name = name,
    icd10_pattern = icd10_combined,
    icd9_pattern = icd9_combined,
    sr_codes = sr_codes
  )
}
