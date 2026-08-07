#' @title Create Disease Definition Object
#'
#' @description
#' Helper function to create a standardized disease definition object
#' containing ICD-10/ICD-9 patterns, self-report codes, UK Biobank First
#' Occurrence fields, and optionally a UK Biobank algorithmically-defined
#' outcome date field.
#'
#' @param name Full disease name (e.g., "Aortic Aneurysm"). If NULL,
#'   defaults to "Custom disease".
#' @param icd10_pattern Regular expression pattern for ICD-10 codes (optional).
#' @param icd9_pattern Regular expression pattern for ICD-9 codes (optional).
#' @param sr_codes Integer vector of UKB self-report illness codes (optional).
#' @param death_icd10 Optional regular expression pattern (or code vector) for
#'   death-cause ICD-10 matching. If NULL, defaults to \code{icd10_pattern}.
#' @param opcs4_pattern Optional regular expression pattern (or code vector) for
#'   OPCS4 operative procedure matching. If NULL, operative procedures are not
#'   used in case ascertainment.
#' @param first_occurrence_fields Optional integer vector of UK Biobank First
#'   Occurrence date field IDs. These fields are generated for 3-character ICD-10
#'   codes in Category 1712, e.g. 131298 for I21 (acute myocardial infarction)
#'   and 130708 for E11 (type 2 diabetes). The source field is normally the next
#'   field ID and is inferred automatically.
#' @param first_occurrence_source_fields Optional integer vector of First
#'   Occurrence source field IDs. If NULL, uses \code{first_occurrence_fields + 1}.
#' @param cancer_icd10_pattern Optional regular expression pattern for UKB
#'   cancer registry ICD-10 codes (Field 40006).
#' @param cancer_histology Optional integer vector of tumour histology codes
#'   (Field 40011) to retain.
#' @param cancer_behaviour Optional integer vector of tumour behaviour codes
#'   (Field 40012) to retain. Use \code{3L} for malignant tumours.
#' @param algo_date_field Integer vector. UKB field ID(s) for
#'   algorithmically-defined outcome dates (Category 42). For example, 42016
#'   for COPD and 42014 for Asthma.
#'   The corresponding data column can be \code{p{field}_i0} or \code{p{field}}.
#'   Records with date \code{1900-01-01} are treated as unknown and excluded.
#' @param algo_source_field Optional integer vector paired with
#'   \code{algo_date_field}. For example, 42017 for COPD source and 42015 for
#'   Asthma source. Use \code{NA} when a date field has no source field.
#' @param icd10 Deprecated alias of \code{icd10_pattern}.
#' @param icd9 Deprecated alias of \code{icd9_pattern}.
#' @param self_report Deprecated alias of \code{sr_codes}.
#'
#' @return A list containing the disease definition parameters.
#'
#' @export
create_disease_definition <- function(name = NULL,
                                      icd10_pattern = NULL,
                                      icd9_pattern = NULL,
                                      sr_codes = NULL,
                                      death_icd10 = NULL,
                                      opcs4_pattern = NULL,
                                      first_occurrence_fields = NULL,
                                      first_occurrence_source_fields = NULL,
                                      cancer_icd10_pattern = NULL,
                                      cancer_histology = NULL,
                                      cancer_behaviour = NULL,
                                      algo_date_field = NULL,
                                      algo_source_field = NULL,
                                      icd10 = NULL,
                                      icd9 = NULL,
                                      self_report = NULL) {

  .normalize_pattern <- function(x, field_name) {
    if (is.null(x)) return(NULL)
    if (!is.character(x)) {
      stop(sprintf("'%s' must be character (regex or code vector)", field_name))
    }
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) return(NULL)
    if (length(x) == 1) return(x)

    # Treat vector inputs as code prefixes unless already anchored.
    x <- ifelse(grepl("^\\^", x), x, paste0("^", x))
    paste0("(", paste(x, collapse = "|"), ")")
  }

  if (!is.null(icd10)) {
    warning("'icd10' is deprecated; use 'icd10_pattern'", call. = FALSE)
    if (is.null(icd10_pattern)) icd10_pattern <- icd10
  }

  if (!is.null(icd9)) {
    warning("'icd9' is deprecated; use 'icd9_pattern'", call. = FALSE)
    if (is.null(icd9_pattern)) icd9_pattern <- icd9
  }

  if (!is.null(self_report)) {
    warning("'self_report' is deprecated; use 'sr_codes'", call. = FALSE)
    if (is.null(sr_codes)) sr_codes <- self_report
  }

  if (is.null(name) || !is.character(name) || length(name) != 1 || !nzchar(name)) {
    name <- "Custom disease"
  }

  icd10_pattern <- .normalize_pattern(icd10_pattern, "icd10_pattern")
  icd9_pattern <- .normalize_pattern(icd9_pattern, "icd9_pattern")
  death_icd10 <- .normalize_pattern(death_icd10, "death_icd10")
  opcs4_pattern <- .normalize_pattern(opcs4_pattern, "opcs4_pattern")
  cancer_icd10_pattern <- .normalize_pattern(cancer_icd10_pattern, "cancer_icd10_pattern")
  first_occurrence_fields <- .normalize_first_occurrence_fields(first_occurrence_fields)
  if (!is.null(first_occurrence_source_fields)) {
    first_occurrence_source_fields <- .normalize_first_occurrence_fields(
      first_occurrence_source_fields,
      "first_occurrence_source_fields"
    )
    if (length(first_occurrence_source_fields) != length(first_occurrence_fields)) {
      stop(
        "'first_occurrence_source_fields' must have the same length as 'first_occurrence_fields'",
        call. = FALSE
      )
    }
  }

  # Backward-compatible default: death matching follows ICD-10 pattern unless specified.
  if (is.null(death_icd10)) {
    death_icd10 <- icd10_pattern
  }
  algorithm_fields <- .normalize_algorithm_field_pairs(
    algo_date_field,
    algo_source_field
  )

  list(
    name = name,
    icd10_pattern = icd10_pattern,
    icd9_pattern = icd9_pattern,
    sr_codes = sr_codes,
    death_icd10 = death_icd10,
    opcs4_pattern = opcs4_pattern,
    first_occurrence_fields = first_occurrence_fields,
    first_occurrence_source_fields = first_occurrence_source_fields,
    cancer_icd10_pattern = cancer_icd10_pattern,
    cancer_histology = .extract_integer_code(cancer_histology),
    cancer_behaviour = .extract_integer_code(cancer_behaviour),
    algo_date_field = algorithm_fields$date,
    algo_source_field = algorithm_fields$source
  )
}

.normalize_algorithm_field_pairs <- function(date_field, source_field = NULL) {
  if (is.null(date_field) || length(date_field) == 0L) {
    if (!is.null(source_field) && length(source_field) > 0L) {
      stop(
        "'algo_source_field' cannot be supplied without 'algo_date_field'.",
        call. = FALSE
      )
    }
    return(list(date = NULL, source = NULL))
  }

  date_numeric <- suppressWarnings(as.numeric(as.character(date_field)))
  date_integer <- suppressWarnings(as.integer(date_numeric))
  if (
    anyNA(date_integer) ||
      any(!is.finite(date_numeric)) ||
      any(date_numeric != date_integer) ||
      any(date_integer <= 0L)
  ) {
    stop("'algo_date_field' must contain positive integer field IDs.", call. = FALSE)
  }

  if (is.null(source_field) || length(source_field) == 0L) {
    source_integer <- rep(NA_integer_, length(date_integer))
  } else {
    if (length(source_field) != length(date_integer)) {
      stop(
        "'algo_source_field' must have the same length as 'algo_date_field'.",
        call. = FALSE
      )
    }
    source_numeric <- suppressWarnings(as.numeric(as.character(source_field)))
    source_integer <- suppressWarnings(as.integer(source_numeric))
    invalid_source <- !is.na(source_field) & (
      is.na(source_integer) |
        !is.finite(source_numeric) |
        source_numeric != source_integer |
        source_integer <= 0L
    )
    if (any(invalid_source)) {
      stop(
        "'algo_source_field' must contain positive integer field IDs or NA.",
        call. = FALSE
      )
    }
  }

  pair_key <- paste(date_integer, ifelse(is.na(source_integer), "NA", source_integer))
  keep <- !duplicated(pair_key)
  date_integer <- date_integer[keep]
  source_integer <- source_integer[keep]

  list(
    date = date_integer,
    source = if (all(is.na(source_integer))) NULL else source_integer
  )
}


#' @title Get Predefined Disease Definitions
#'
#' @description
#' Returns a list of commonly used cardiovascular and metabolic disease
#' definitions with validated ICD-10, ICD-9, and self-report code mappings.
#'
#' @param source Definition source. `"curated"` returns the original manually
#'   curated UKBAnalytica definitions. `"pomegranate"` returns definitions
#'   converted from the built-in Pomegranate-derived disease catalog. `"both"`
#'   returns diseases that can be matched between both sources, with standardized
#'   curated names and either intersected or unioned source definitions depending
#'   on `merge_type`.
#' @param merge_type Merge strategy for `source = "both"`. `"intersection"`
#'   keeps codes supported by both curated and Pomegranate definitions. `"union"`
#'   combines codes from both definitions.
#' @param disease Optional disease key or name pattern used to subset the
#'   returned definition list.
#' @param supported_only Logical. For Pomegranate-derived definitions, keep only
#'   code systems currently supported by UKBAnalytica parsers.
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
#'   \item{T1DM}{Type 1 Diabetes Mellitus}
#'   \item{T2DM}{Type 2 Diabetes Mellitus}
#'   \item{Vascular_Disease}{Peripheral vascular disease}
#'   \item{Arrhythmia}{Broad cardiac arrhythmia endpoint including OPCS4 procedures}
#'   \item{Atrial_Fibrillation}{Atrial arrhythmia / atrial fibrillation-flutter}
#'   \item{Ventricular_Arrhythmia}{Ventricular arrhythmia endpoint}
#'   \item{AV_Block}{Atrioventricular conduction block}
#'   \item{Intraventricular_Block}{Intraventricular conduction block}
#'   \item{SVT}{Supraventricular tachycardia}
#'   \item{Lung_Cancer}{Lung cancer using ICD-10/death and cancer registry}
#'   \item{Additional chronic diseases}{Common respiratory, renal,
#'   gastrointestinal, neurologic, psychiatric, eye, skin, musculoskeletal, and
#'   cancer endpoints used in UKB epidemiology workflows}
#' }
#'
#' @export
get_predefined_diseases <- function(source = c("curated", "pomegranate", "both"),
                                    merge_type = c("intersection", "union"),
                                    disease = NULL,
                                    supported_only = TRUE) {
  source <- match.arg(source)
  merge_type <- match.arg(merge_type)

  curated <- list(
    # Aortic diseases
    AA = create_disease_definition(
      name = "Aortic Aneurysm",
      icd10_pattern = "^I71",
      icd9_pattern = "^441",
      first_occurrence_fields = 131382
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
      sr_codes = c(
        1066, 1067, 1074, 1075, 1076, 1077, 1078, 1079,
        1081, 1082, 1086, 1087, 1425, 1426, 1471, 1479,
        1483, 1484, 1485, 1486, 1487, 1488, 1489, 1490,
        1491, 1492, 1583, 1584, 1585, 1586, 1587, 1588,
        1589, 1590, 1591
      ),
      first_occurrence_fields = c(131298, 131300, 131302, 131304, 131306)
    ),
    MI = create_disease_definition(
      name = "Myocardial Infarction",
      icd10_pattern = "^(I21|I22)",
      icd9_pattern = "^410",
      sr_codes = c(1066, 1075),
      first_occurrence_fields = c(131298, 131300),
      algo_date_field = 42000,
      algo_source_field = 42001
    ),
    STEMI = create_disease_definition(
      name = "ST-Elevation Myocardial Infarction",
      algo_date_field = 42002,
      algo_source_field = 42003
    ),
    NSTEMI = create_disease_definition(
      name = "Non-ST-Elevation Myocardial Infarction",
      algo_date_field = 42004,
      algo_source_field = 42005
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
      sr_codes = c(1068),
      first_occurrence_fields = c(131360, 131362, 131364, 131366, 131368),
      algo_date_field = 42006,
      algo_source_field = 42007
    ),
    Ischaemic_Stroke = create_disease_definition(
      name = "Ischaemic Stroke",
      algo_date_field = 42008,
      algo_source_field = 42009
    ),
    Intracerebral_Haemorrhage = create_disease_definition(
      name = "Intracerebral Haemorrhage",
      algo_date_field = 42010,
      algo_source_field = 42011
    ),
    Subarachnoid_Haemorrhage = create_disease_definition(
      name = "Subarachnoid Haemorrhage",
      algo_date_field = 42012,
      algo_source_field = 42013
    ),

    # Metabolic diseases
    Hypertension = create_disease_definition(
      name = "Hypertension",
      icd10_pattern = "^(I10|I11|I12|I13|I14|I15)",
      icd9_pattern = "^(401|402|403|404|405)",
      sr_codes = c(1065, 1072),
      first_occurrence_fields = c(131286, 131288, 131290, 131292, 131294)
    ),
    Diabetes = create_disease_definition(
      name = "Diabetes Mellitus",
      icd10_pattern = "^(E10|E11|E12|E13|E14|G590|G632|H280|H360|M142|N083|O24|P702|T383|Y423)",
      icd9_pattern = "^(249|250|3572|3620|6480|7751|9623)",
      sr_codes = c(1220, 1221, 1222, 1223, 1276, 1468, 1607),
      first_occurrence_fields = c(130706, 130708, 130710, 130712, 130714)
    ),
    T1DM = create_disease_definition(
      name = "Type 1 Diabetes Mellitus",
      icd10_pattern = "^E10",
      sr_codes = c(1222),
      first_occurrence_fields = 130706
    ),
    T2DM = create_disease_definition(
      name = "Type 2 Diabetes Mellitus",
      icd10_pattern = "^E11",
      sr_codes = c(1223),
      first_occurrence_fields = 130708
    ),
    Hyperglycaemia = create_disease_definition(
      name = "Hyperglycaemia",
      icd10_pattern = "^R73",
      icd9_pattern = "^7902"
    ),
    PCOS = create_disease_definition(
      name = "Polycystic Ovary Syndrome",
      sr_codes = c(1350)
    ),

    # Vascular diseases
    Vascular_Disease = create_disease_definition(
      name = "Vascular Disease",
      icd10_pattern = "^(I71|I72|I73|I77|I78|I79)",
      icd9_pattern = "^(441|442|443|447)",
      first_occurrence_fields = c(131382, 131384, 131386, 131390, 131392, 131394)
    ),

    # Coronary and rhythm disorders
    Arrhythmia = create_disease_definition(
      name = "Cardiac Arrhythmia",
      icd10_pattern = "^(I44|I45|I46|I47|I48|I49)",
      opcs4_pattern = "^(K576|K59|K60|K61|K62|K641|K72|K73|K74)",
      first_occurrence_fields = c(131342, 131344, 131346, 131348, 131350, 131352)
    ),
    Angina = create_disease_definition(
      name = "Angina Pectoris",
      icd10_pattern = "^I20",
      icd9_pattern = "^413",
      sr_codes = c(1074),
      first_occurrence_fields = 131296
    ),
    Atrial_Fibrillation = create_disease_definition(
      name = "Atrial Fibrillation/Flutter",
      icd10_pattern = "^I48",
      icd9_pattern = "^(4273|4274)",
      sr_codes = c(1471, 1483, 1485),
      opcs4_pattern = "^K62",
      first_occurrence_fields = 131350
    ),
    Ventricular_Arrhythmia = create_disease_definition(
      name = "Ventricular Arrhythmia",
      icd10_pattern = "^(I470|I472|I490|I493)",
      opcs4_pattern = "^(K576|K641)"
    ),
    AV_Block = create_disease_definition(
      name = "Atrioventricular Block",
      icd10_pattern = "^(I440|I441|I442|I443|I458)"
    ),
    Intraventricular_Block = create_disease_definition(
      name = "Intraventricular Conduction Block",
      icd10_pattern = "^(I444|I445|I446|I447|I450|I451|I452|I453|I454)"
    ),
    SVT = create_disease_definition(
      name = "Supraventricular Tachycardia",
      icd10_pattern = "^I471"
    ),

    # Respiratory diseases
    Asthma = create_disease_definition(
      name = "Asthma",
      icd10_pattern = "^(J45|J46)",
      icd9_pattern = "^493",
      sr_codes = c(1111),
      first_occurrence_fields = c(131494, 131496),
      algo_date_field = 42014,
      algo_source_field = 42015
    ),
    COPD = create_disease_definition(
      name = "Chronic Obstructive Pulmonary Disease",
      icd10_pattern = "^(J40|J41|J42|J43|J44)",
      icd9_pattern = "^(491|492|4932|496)",
      sr_codes = c(1112, 1113, 1472),
      first_occurrence_fields = c(131484, 131486, 131488, 131490, 131492),
      algo_date_field = 42016,
      algo_source_field = 42017
    ),
    Lung_Cancer = create_disease_definition(
      name = "Lung Cancer",
      icd10_pattern = "^C34",
      death_icd10 = "^C34",
      cancer_icd10_pattern = "^C34",
      cancer_behaviour = 3L
    ),

    # Renal and metabolic
    CKD = create_disease_definition(
      name = "Chronic Kidney Disease",
      icd10_pattern = "^(N18|N19)",
      icd9_pattern = "^(585|586)",
      sr_codes = c(1192, 1193, 1194, 1405, 1582, 1675),
      first_occurrence_fields = c(132032, 132034)
    ),
    ESRD = create_disease_definition(
      name = "End Stage Renal Disease",
      algo_date_field = 42026,
      algo_source_field = 42027
    ),
    Hyperlipidemia = create_disease_definition(
      name = "Hyperlipidemia/High Cholesterol",
      icd10_pattern = "^E78",
      icd9_pattern = "^272",
      sr_codes = c(1473),
      first_occurrence_fields = 130814
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
      sr_codes = c(1068, 1093, 1094),
      first_occurrence_fields = c(131308, 131396, 131398, 131400)
    ),

    # Musculoskeletal and rheumatologic
    Osteoarthritis = create_disease_definition(
      name = "Osteoarthritis",
      icd10_pattern = "^(M15|M16|M17|M18|M19)",
      icd9_pattern = "^715",
      sr_codes = c(1465),
      first_occurrence_fields = c(131868, 131870, 131872, 131874, 131876)
    ),
    Rheumatoid_Arthritis = create_disease_definition(
      name = "Rheumatoid Arthritis",
      icd10_pattern = "^(M05|M06)",
      icd9_pattern = "^714",
      sr_codes = c(1464, 1477),
      first_occurrence_fields = c(131848, 131850)
    ),

    # Neurologic and psychiatric
    Parkinsons = create_disease_definition(
      name = "Parkinson's Disease",
      icd10_pattern = "^G20",
      icd9_pattern = "^3320",
      sr_codes = c(1262),
      first_occurrence_fields = 131022,
      algo_date_field = 42032,
      algo_source_field = 42033
    ),
    Parkinsonism = create_disease_definition(
      name = "All-Cause Parkinsonism",
      algo_date_field = 42030,
      algo_source_field = 42031
    ),
    Progressive_Supranuclear_Palsy = create_disease_definition(
      name = "Progressive Supranuclear Palsy",
      algo_date_field = 42034,
      algo_source_field = 42035
    ),
    Multiple_System_Atrophy = create_disease_definition(
      name = "Multiple System Atrophy",
      algo_date_field = 42036,
      algo_source_field = 42037
    ),
    Dementia = create_disease_definition(
      name = "Dementia/Alzheimer's Disease",
      icd10_pattern = "^(F00|F01|F02|F03|G30)",
      icd9_pattern = "^(290|3310)",
      sr_codes = c(1263),
      first_occurrence_fields = c(130836, 130838, 130840, 130842, 131036),
      algo_date_field = 42018,
      algo_source_field = 42019
    ),
    Alzheimers_Disease = create_disease_definition(
      name = "Alzheimer's Disease",
      algo_date_field = 42020,
      algo_source_field = 42021
    ),
    Vascular_Dementia = create_disease_definition(
      name = "Vascular Dementia",
      algo_date_field = 42022,
      algo_source_field = 42023
    ),
    Frontotemporal_Dementia = create_disease_definition(
      name = "Frontotemporal Dementia",
      algo_date_field = 42024,
      algo_source_field = 42025
    ),
    Motor_Neurone_Disease = create_disease_definition(
      name = "Motor Neurone Disease",
      algo_date_field = 42028,
      algo_source_field = 42029
    ),
    Epilepsy = create_disease_definition(
      name = "Epilepsy",
      icd10_pattern = "^(G40|G41)",
      icd9_pattern = "^345",
      sr_codes = c(1264),
      first_occurrence_fields = c(131048, 131050)
    ),
    Depression = create_disease_definition(
      name = "Depressive Disorders",
      icd10_pattern = "^(F32|F33)",
      icd9_pattern = "^(2962|2963|311)",
      sr_codes = c(1286, 1531, 1682),
      first_occurrence_fields = c(130894, 130896)
    ),
    Anxiety = create_disease_definition(
      name = "Anxiety Disorders",
      icd10_pattern = "^(F40|F41)",
      icd9_pattern = "^(3000|3002|3003)",
      sr_codes = c(1287),
      first_occurrence_fields = c(130904, 130906)
    ),

    # Additional common chronic diseases
    Stroke_TIA = create_disease_definition(
      name = "Stroke or Transient Ischaemic Attack",
      icd10_pattern = "^(I60|I61|I62|I63|I64|I65|I66|G45)",
      icd9_pattern = "^(430|431|432|433|434|435|436|437)",
      sr_codes = c(1081, 1082, 1083, 1086, 1583)
    ),
    Thyroid_Disorders = create_disease_definition(
      name = "Thyroid Disorders",
      icd10_pattern = "^(E00|E01|E02|E03|E04|E05|E06|E07)",
      icd9_pattern = "^(240|241|242|243|244|245|246)"
    ),
    Dyspepsia = create_disease_definition(
      name = "Dyspepsia and Upper Gastrointestinal Disorders",
      icd10_pattern = "^(K21|K22|K25|K26|K27|K28|K29|K30)",
      icd9_pattern = "^(530|531|532|533|534|535|536)"
    ),
    Irritable_Bowel_Syndrome = create_disease_definition(
      name = "Irritable Bowel Syndrome",
      icd10_pattern = "^K58",
      icd9_pattern = "^5641",
      sr_codes = c(1154)
    ),
    Inflammatory_Bowel_Disease = create_disease_definition(
      name = "Inflammatory Bowel Disease",
      icd10_pattern = "^(K50|K51)",
      icd9_pattern = "^(555|556)",
      sr_codes = c(1461, 1462, 1463)
    ),
    Diverticular_Disease = create_disease_definition(
      name = "Diverticular Disease",
      icd10_pattern = "^K57",
      icd9_pattern = "^562",
      sr_codes = c(1458)
    ),
    Treated_Constipation = create_disease_definition(
      name = "Treated Constipation",
      icd10_pattern = "^K590",
      icd9_pattern = "^5640",
      sr_codes = c(1599)
    ),
    Chronic_Liver_Disease = create_disease_definition(
      name = "Chronic Liver Disease",
      icd10_pattern = "^(K70|K71|K72|K73|K74|K75|K76)",
      icd9_pattern = "^(571|572|573)",
      sr_codes = c(1141, 1157, 1158, 1506)
    ),
    Alcohol_Use_Disorder = create_disease_definition(
      name = "Alcohol Use Disorder",
      icd10_pattern = "^(F10|K70|T51)",
      icd9_pattern = "^(291|303|3050|5710|5711|5712|5713)"
    ),
    Substance_Use_Disorder = create_disease_definition(
      name = "Other Psychoactive Substance Use Disorder",
      icd10_pattern = "^(F11|F12|F13|F14|F15|F16|F18|F19)",
      icd9_pattern = "^(304|305[2-9])"
    ),
    Schizophrenia_Bipolar = create_disease_definition(
      name = "Schizophrenia or Bipolar Disorder",
      icd10_pattern = "^(F20|F21|F22|F23|F24|F25|F30|F31)",
      icd9_pattern = "^(295|2960|2961|2964|2965|2966|2967|2968)",
      sr_codes = c(1289, 1291)
    ),
    Migraine = create_disease_definition(
      name = "Migraine",
      icd10_pattern = "^G43",
      icd9_pattern = "^346",
      sr_codes = c(1265)
    ),
    Bronchiectasis = create_disease_definition(
      name = "Bronchiectasis",
      icd10_pattern = "^J47",
      icd9_pattern = "^494",
      sr_codes = c(1114)
    ),
    Multiple_Sclerosis = create_disease_definition(
      name = "Multiple Sclerosis",
      icd10_pattern = "^G35",
      icd9_pattern = "^340",
      sr_codes = c(1261)
    ),
    Menieres_Disease = create_disease_definition(
      name = "Meniere's Disease",
      icd10_pattern = "^H810",
      icd9_pattern = "^3860",
      sr_codes = c(1421)
    ),
    Pernicious_Anaemia = create_disease_definition(
      name = "Pernicious Anaemia",
      icd10_pattern = "^D51",
      icd9_pattern = "^2810",
      sr_codes = c(1331)
    ),
    Psoriasis_Eczema = create_disease_definition(
      name = "Psoriasis or Eczema",
      icd10_pattern = "^(L20|L21|L22|L23|L24|L25|L26|L27|L30|L40|L41)",
      icd9_pattern = "^(691|692|696)",
      sr_codes = c(1452, 1453)
    ),
    Fracture = create_disease_definition(
      name = "Major Fracture",
      icd10_pattern = "^(S32|S42|S52|S72|S82)",
      icd9_pattern = "^(800|801|802|803|804|805|806|807|808|809|810|811|812|813|814|815|816|817|818|819|820|821|822|823|824|825|826|827|828|829)",
      sr_codes = c(1647, 1648, 1650)
    ),
    Glaucoma = create_disease_definition(
      name = "Glaucoma",
      icd10_pattern = "^H40",
      icd9_pattern = "^365",
      sr_codes = c(1277)
    ),
    Cataract = create_disease_definition(
      name = "Cataract",
      icd10_pattern = "^(H25|H26|H28)",
      icd9_pattern = "^366",
      sr_codes = c(1278)
    ),
    AMD = create_disease_definition(
      name = "Age-Related Macular Degeneration",
      icd10_pattern = "^H353",
      icd9_pattern = "^3625",
      sr_codes = c(1528)
    ),
    Prostate_Disorders = create_disease_definition(
      name = "Benign Prostate Disorders",
      icd10_pattern = "^(N40|N41|N42)",
      icd9_pattern = "^600",
      sr_codes = c(1207, 1396, 1516)
    ),
    Breast_Cancer = create_disease_definition(
      name = "Breast Cancer",
      icd10_pattern = "^C50",
      death_icd10 = "^C50",
      cancer_icd10_pattern = "^C50",
      cancer_behaviour = 3L
    ),
    Prostate_Cancer = create_disease_definition(
      name = "Prostate Cancer",
      icd10_pattern = "^C61",
      death_icd10 = "^C61",
      cancer_icd10_pattern = "^C61",
      cancer_behaviour = 3L
    ),
    Colorectal_Cancer = create_disease_definition(
      name = "Colorectal Cancer",
      icd10_pattern = "^(C18|C19|C20)",
      death_icd10 = "^(C18|C19|C20)",
      cancer_icd10_pattern = "^(C18|C19|C20)",
      cancer_behaviour = 3L
    ),
    Colon_Cancer = create_disease_definition(
      name = "Colon Cancer",
      icd10_pattern = "^C18",
      death_icd10 = "^C18",
      cancer_icd10_pattern = "^C18",
      cancer_behaviour = 3L
    ),
    Rectal_Cancer = create_disease_definition(
      name = "Rectal Cancer",
      icd10_pattern = "^(C19|C20)",
      death_icd10 = "^(C19|C20)",
      cancer_icd10_pattern = "^(C19|C20)",
      cancer_behaviour = 3L
    ),
    Melanoma = create_disease_definition(
      name = "Melanoma",
      icd10_pattern = "^C43",
      death_icd10 = "^C43",
      cancer_icd10_pattern = "^C43",
      cancer_behaviour = 3L
    ),
    Non_Melanoma_Skin_Cancer = create_disease_definition(
      name = "Non-Melanoma Skin Cancer",
      icd10_pattern = "^C44",
      death_icd10 = "^C44",
      cancer_icd10_pattern = "^C44",
      cancer_behaviour = 3L
    ),
    Ovarian_Cancer = create_disease_definition(
      name = "Ovarian Cancer",
      icd10_pattern = "^C56",
      death_icd10 = "^C56",
      cancer_icd10_pattern = "^C56",
      cancer_behaviour = 3L
    ),
    Uterus_Cancer = create_disease_definition(
      name = "Uterus Cancer",
      icd10_pattern = "^C54",
      death_icd10 = "^C54",
      cancer_icd10_pattern = "^C54",
      cancer_behaviour = 3L
    ),
    Oesophageal_Cancer = create_disease_definition(
      name = "Oesophageal Cancer",
      icd10_pattern = "^C15",
      death_icd10 = "^C15",
      cancer_icd10_pattern = "^C15",
      cancer_behaviour = 3L
    ),
    Stomach_Cancer = create_disease_definition(
      name = "Stomach Cancer",
      icd10_pattern = "^C16",
      death_icd10 = "^C16",
      cancer_icd10_pattern = "^C16",
      cancer_behaviour = 3L
    ),
    Erectile_Dysfunction = create_disease_definition(
      name = "Erectile Dysfunction",
      icd10_pattern = "^(F52|N48)"
    ),
    MACE = create_disease_definition(
      name = "Major Adverse Cardiovascular Events",
      icd10_pattern = "^(G45|I21|I22|I23|I24|I25|I63|I64)"
    ),
    Renal_Disease = create_disease_definition(
      name = "Renal Disease",
      icd10_pattern = paste0(
        "^(N00|N01|N02|N03|N04|N05|N06|N07|N08|N09|",
        "N10|N11|N12|N13|N14|N15|N16|N17|N18|N19|",
        "N25|N26|N27|N28|N29)"
      )
    ),
    Severe_Mental_Illness = create_disease_definition(
      name = "Severe Mental Illness",
      icd10_pattern = "^(F20|F25|F30|F31|F32|F33|F44)"
    ),
    Systemic_Lupus_Erythematosus = create_disease_definition(
      name = "Systemic Lupus Erythematosus",
      icd10_pattern = "^M32"
    )
  )

  if (identical(source, "curated")) {
    return(.filter_disease_definition_list(curated, disease))
  }

  pomegranate <- get_pomegranate_diseases(
    disease = if (identical(source, "pomegranate")) disease else NULL,
    supported_only = supported_only
  )
  if (identical(source, "pomegranate")) {
    return(pomegranate)
  }

  out <- .merge_predefined_diseases(curated, pomegranate, merge_type = merge_type)
  .filter_disease_definition_list(out, disease)
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

  # Combine death ICD-10 patterns (fallback to general ICD-10 when absent)
  death_patterns <- sapply(defs, function(x) {
    if (!is.null(x$death_icd10)) x$death_icd10 else x$icd10_pattern
  })
  death_patterns <- death_patterns[!sapply(death_patterns, is.null)]
  death_combined <- if (length(death_patterns) > 0) {
    paste0("(", paste(death_patterns, collapse = "|"), ")")
  } else NULL

  # Combine OPCS4 patterns
  opcs4_patterns <- sapply(defs, function(x) x$opcs4_pattern)
  opcs4_patterns <- opcs4_patterns[!sapply(opcs4_patterns, is.null)]
  opcs4_combined <- if (length(opcs4_patterns) > 0) {
    paste0("(", paste(opcs4_patterns, collapse = "|"), ")")
  } else NULL

  # Combine UKB First Occurrence date/source fields as paired metadata.
  first_occurrence_pairs <- rbindlist(lapply(defs, function(x) {
    date_fields <- x$first_occurrence_fields
    if (is.null(date_fields) || length(date_fields) == 0L) {
      return(NULL)
    }
    source_fields <- x$first_occurrence_source_fields
    if (is.null(source_fields)) {
      source_fields <- date_fields + 1L
    }
    data.table(
      date_field = as.integer(date_fields),
      source_field = as.integer(source_fields)
    )
  }), use.names = TRUE, fill = TRUE)
  if (nrow(first_occurrence_pairs) > 0L) {
    first_occurrence_pairs <- unique(first_occurrence_pairs)
    first_occurrence_fields <- first_occurrence_pairs$date_field
    first_occurrence_source_fields <- first_occurrence_pairs$source_field
  } else {
    first_occurrence_fields <- NULL
    first_occurrence_source_fields <- NULL
  }

  # Combine cancer registry definitions. Histology/behaviour restrictions are
  # preserved only when every cancer definition supplies that restriction.
  cancer_patterns <- sapply(defs, function(x) x$cancer_icd10_pattern)
  cancer_patterns <- cancer_patterns[!sapply(cancer_patterns, is.null)]
  cancer_combined <- if (length(cancer_patterns) > 0) {
    paste0("(", paste(cancer_patterns, collapse = "|"), ")")
  } else NULL

  cancer_histology_list <- lapply(defs, function(x) x$cancer_histology)
  cancer_histology_nonnull <- cancer_histology_list[!vapply(cancer_histology_list, is.null, logical(1))]
  cancer_histology <- if (length(cancer_histology_nonnull) == length(cancer_patterns) && length(cancer_patterns) > 0) {
    unique(unlist(cancer_histology_nonnull))
  } else NULL

  cancer_behaviour_list <- lapply(defs, function(x) x$cancer_behaviour)
  cancer_behaviour_nonnull <- cancer_behaviour_list[!vapply(cancer_behaviour_list, is.null, logical(1))]
  cancer_behaviour <- if (length(cancer_behaviour_nonnull) == length(cancer_patterns) && length(cancer_patterns) > 0) {
    unique(unlist(cancer_behaviour_nonnull))
  } else NULL

  # Combine self-report codes
  sr_codes <- unlist(lapply(defs, function(x) x$sr_codes))
  sr_codes <- unique(sr_codes[!is.na(sr_codes)])
  if (length(sr_codes) == 0) sr_codes <- NULL

  # Preserve every algorithm date/source pair in composite definitions.
  algorithm_pairs <- rbindlist(lapply(defs, function(x) {
    date_fields <- x$algo_date_field
    if (is.null(date_fields) || length(date_fields) == 0L) {
      return(NULL)
    }
    source_fields <- x$algo_source_field
    if (is.null(source_fields)) {
      source_fields <- rep(NA_integer_, length(date_fields))
    }
    data.table(
      date_field = as.integer(date_fields),
      source_field = as.integer(source_fields)
    )
  }), use.names = TRUE, fill = TRUE)
  if (nrow(algorithm_pairs) > 0L) {
    algorithm_pairs <- unique(algorithm_pairs)
    algo_date_field <- algorithm_pairs$date_field
    algo_source_field <- algorithm_pairs$source_field
  } else {
    algo_date_field <- NULL
    algo_source_field <- NULL
  }

  create_disease_definition(
    name = name,
    icd10_pattern = icd10_combined,
    icd9_pattern = icd9_combined,
    sr_codes = sr_codes,
    death_icd10 = death_combined,
    opcs4_pattern = opcs4_combined,
    first_occurrence_fields = first_occurrence_fields,
    first_occurrence_source_fields = first_occurrence_source_fields,
    cancer_icd10_pattern = cancer_combined,
    cancer_histology = cancer_histology,
    cancer_behaviour = cancer_behaviour,
    algo_date_field = algo_date_field,
    algo_source_field = algo_source_field
  )
}
