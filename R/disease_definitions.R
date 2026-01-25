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
    ),

    # Coronary and rhythm disorders
    Angina = create_disease_definition(
      name = "Angina Pectoris",
      icd10_pattern = "^I20",
      icd9_pattern = "^413",
      sr_codes = c(1074)
    ),
    Atrial_Fibrillation = create_disease_definition(
      name = "Atrial Fibrillation/Flutter",
      icd10_pattern = "^I48",
      icd9_pattern = "^(4273|4274)",
      sr_codes = c(1471, 1483, 1485)
    ),

    # Respiratory diseases
    Asthma = create_disease_definition(
      name = "Asthma",
      icd10_pattern = "^(J45|J46)",
      icd9_pattern = "^493",
      sr_codes = c(1111)
    ),
    COPD = create_disease_definition(
      name = "Chronic Obstructive Pulmonary Disease",
      icd10_pattern = "^(J41|J42|J43|J44)",
      icd9_pattern = "^(491|492|496)",
      sr_codes = c(1112, 1113, 1472)
    ),

    # Renal and metabolic
    CKD = create_disease_definition(
      name = "Chronic Kidney Disease",
      icd10_pattern = "^(N18|N19)",
      icd9_pattern = "^(585|586)",
      sr_codes = c(1192, 1193, 1194, 1405, 1582, 1675)
    ),
    Hyperlipidemia = create_disease_definition(
      name = "Hyperlipidemia/High Cholesterol",
      icd10_pattern = "^E78",
      icd9_pattern = "^272",
      sr_codes = c(1473)
    ),

    # Peripheral vascular and thromboembolic
    PAD = create_disease_definition(
      name = "Peripheral Arterial Disease",
      icd10_pattern = "^(I702|I703|I704|I708|I709|I738|I739|I74)",
      icd9_pattern = "^(4402|4403|4409|4439)",
      sr_codes = c(1067, 1104, 1105)
    ),
    VTE = create_disease_definition(
      name = "Venous Thromboembolism",
      icd10_pattern = "^(I26|I80|I81|I82)",
      icd9_pattern = "^(4151|451|453)",
      sr_codes = c(1068, 1093, 1094)
    ),

    # Musculoskeletal and rheumatologic
    Osteoarthritis = create_disease_definition(
      name = "Osteoarthritis",
      icd10_pattern = "^(M15|M16|M17|M18|M19)",
      icd9_pattern = "^715",
      sr_codes = c(1465)
    ),
    Rheumatoid_Arthritis = create_disease_definition(
      name = "Rheumatoid Arthritis",
      icd10_pattern = "^(M05|M06)",
      icd9_pattern = "^714",
      sr_codes = c(1464, 1477)
    ),

    # Neurologic and psychiatric
    Parkinsons = create_disease_definition(
      name = "Parkinson's Disease",
      icd10_pattern = "^G20",
      icd9_pattern = "^3320",
      sr_codes = c(1262)
    ),
    Dementia = create_disease_definition(
      name = "Dementia/Alzheimer's Disease",
      icd10_pattern = "^(F00|F01|F02|F03|G30)",
      icd9_pattern = "^(290|3310)",
      sr_codes = c(1263)
    ),
    Epilepsy = create_disease_definition(
      name = "Epilepsy",
      icd10_pattern = "^(G40|G41)",
      icd9_pattern = "^345",
      sr_codes = c(1264)
    ),
    Depression = create_disease_definition(
      name = "Depressive Disorders",
      icd10_pattern = "^(F32|F33)",
      icd9_pattern = "^(2962|2963|311)",
      sr_codes = c(1286, 1531, 1682)
    ),
    Anxiety = create_disease_definition(
      name = "Anxiety Disorders",
      icd10_pattern = "^(F40|F41)",
      icd9_pattern = "^(3000|3002|3003)",
      sr_codes = c(1287)
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
