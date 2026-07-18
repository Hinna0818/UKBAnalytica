#' Query the built-in disease code catalog
#'
#' @description
#' Returns a source-aware disease code catalog containing curated
#' UKBAnalytica disease definitions and Pomegranate-derived UK Biobank
#' phenotype coding definitions. This function returns tabular code metadata;
#' it does not change the default behavior of [get_predefined_diseases()].
#'
#' @param source Character. One of `"all"`, `"curated"`, or `"pomegranate"`.
#' @param disease Optional disease name or definition ID pattern.
#' @param code_system Optional code system filter, such as `"ICD-10"` or
#'   `"self-report illness"`.
#' @param supported_only Logical. If TRUE, keep only catalog rows currently
#'   supported by UKBAnalytica disease parsers.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' copd_codes <- get_disease_catalog(disease = "copd")
#' head(copd_codes)
get_disease_catalog <- function(source = c("all", "curated", "pomegranate"),
                                disease = NULL,
                                code_system = NULL,
                                supported_only = FALSE) {
  source <- match.arg(source)
  out <- diseases

  if (!identical(source, "all")) {
    out <- out[out[["source"]] == source, , drop = FALSE]
  }

  if (!is.null(disease)) {
    if (!is.character(disease) || length(disease) != 1 || is.na(disease) || !nzchar(disease)) {
      stop("`disease` must be NULL or a single non-empty character string.", call. = FALSE)
    }
    pattern <- disease
    keep <- grepl(pattern, out[["definition_id"]], ignore.case = TRUE) |
      grepl(pattern, out[["disease_name"]], ignore.case = TRUE)
    out <- out[keep, , drop = FALSE]
  }

  if (!is.null(code_system)) {
    if (!is.character(code_system) || anyNA(code_system) || any(!nzchar(code_system))) {
      stop("`code_system` must be NULL or a non-empty character vector.", call. = FALSE)
    }
    out <- out[out[["code_system"]] %in% code_system, , drop = FALSE]
  }

  if (isTRUE(supported_only)) {
    out <- out[isTRUE(out[["is_supported_by_current_parser"]]) |
      out[["is_supported_by_current_parser"]] %in% TRUE, , drop = FALSE]
  }

  row.names(out) <- NULL
  out
}

#' Query the built-in medication code catalog
#'
#' @description
#' Returns a medication code catalog containing UKBAnalytica curated
#' medication definitions and UK Biobank official coding 4 entries for field
#' 20003.
#'
#' @param medication Optional medication name, ID, or code pattern.
#' @param medication_class Optional medication class filter.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' metformin <- get_medication_catalog("metformin")
#' head(metformin)
get_medication_catalog <- function(medication = NULL,
                                   medication_class = NULL) {
  out <- medications

  if (!is.null(medication)) {
    if (!is.character(medication) || length(medication) != 1 || is.na(medication) || !nzchar(medication)) {
      stop("`medication` must be NULL or a single non-empty character string.", call. = FALSE)
    }
    keep <- grepl(medication, out[["medication_id"]], ignore.case = TRUE) |
      grepl(medication, out[["medication_name"]], ignore.case = TRUE) |
      grepl(medication, out[["code"]], ignore.case = TRUE)
    out <- out[keep, , drop = FALSE]
  }

  if (!is.null(medication_class)) {
    if (!is.character(medication_class) || anyNA(medication_class) || any(!nzchar(medication_class))) {
      stop("`medication_class` must be NULL or a non-empty character vector.", call. = FALSE)
    }
    out <- out[out[["medication_class"]] %in% medication_class, , drop = FALSE]
  }

  row.names(out) <- NULL
  out
}

#' Load the Pomegranate portal coding evidence table
#'
#' @description
#' Loads a long-form Pomegranate portal extraction from a user-supplied local
#' CSV or CSV.GZ file for audit and traceability. The canonical Pomegranate
#' disease catalog used by `get_disease_catalog(source = "pomegranate")` is
#' built into the package from the public GitHub YAML algorithms; the portal
#' audit table is not required for endpoint construction and is not shipped in
#' the CRAN build.
#'
#' @param path Path to a local Pomegranate portal CSV or CSV.GZ file.
#'
#' @return A data.frame.
#' @export
#'
load_pomegranate_portal_coding <- function(path = NULL) {
  if (is.null(path)) {
    stop(
      "`path` is required. The portal audit table is not shipped with the ",
      "package build; use get_disease_catalog(source = 'pomegranate') for the ",
      "built-in canonical catalog.",
      call. = FALSE
    )
  }
  if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Pomegranate portal coding file was not found: ", path, call. = FALSE)
  }
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
}

#' Get the Pomegranate source manifest
#'
#' @description
#' Returns source provenance for the built-in Pomegranate resources, including
#' the GitHub YAML commit used for the canonical disease catalog and the portal
#' CSV retained for audit.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' get_pomegranate_source_manifest()
get_pomegranate_source_manifest <- function() {
  pomegranate_source_manifest
}

#' Get Pomegranate-derived disease definitions
#'
#' @description
#' Converts the Pomegranate-derived rows in the disease catalog into
#' UKBAnalytica disease definition objects. Only code systems currently
#' supported by the package parsers are used by default; GP and medication
#' rows remain available through \code{get_disease_catalog()}.
#'
#' @param disease Optional disease name or definition ID pattern.
#' @param supported_only Logical. If TRUE, use only rows supported by current
#'   UKBAnalytica disease parsers.
#'
#' @return A named list of disease definition objects.
#' @export
#'
#' @examples
#' pom <- get_pomegranate_diseases("asthma")
#' names(pom)
get_pomegranate_diseases <- function(disease = NULL, supported_only = TRUE) {
  x <- get_disease_catalog(
    source = "pomegranate",
    disease = disease,
    supported_only = supported_only
  )
  if (nrow(x) == 0) {
    return(list())
  }

  split_x <- split(x, x[["definition_id"]])
  out <- lapply(split_x, .catalog_to_disease_definition)
  out[!vapply(out, is.null, logical(1))]
}

.catalog_to_disease_definition <- function(x) {
  disease_name <- unique(x[["disease_name"]])
  disease_name <- disease_name[!is.na(disease_name) & nzchar(disease_name)][1]
  if (is.na(disease_name) || !nzchar(disease_name)) {
    return(NULL)
  }

  get_codes <- function(system, field_ids = NULL) {
    y <- x[x[["code_system"]] %in% system, , drop = FALSE]
    if (!is.null(field_ids)) {
      y <- y[y[["field_id"]] %in% field_ids, , drop = FALSE]
    }
    codes <- unique(y[["code"]])
    codes <- codes[!is.na(codes) & nzchar(codes)]
    codes
  }

  icd10 <- get_codes("ICD-10", field_ids = c(41202L, 41204L))
  death_icd10 <- get_codes("death registry ICD-10", field_ids = c(40001L, 40002L))
  sr <- suppressWarnings(as.integer(get_codes("self-report illness", field_ids = 20002L)))
  sr <- sr[!is.na(sr)]
  opcs4 <- get_codes("OPCS4", field_ids = c(41200L, 41210L))
  cancer_icd10 <- get_codes("cancer registry ICD-10", field_ids = 40006L)

  if (
    length(icd10) == 0 &&
      length(death_icd10) == 0 &&
      length(sr) == 0 &&
      length(opcs4) == 0 &&
      length(cancer_icd10) == 0
  ) {
    return(NULL)
  }

  create_disease_definition(
    name = disease_name,
    icd10_pattern = .catalog_codes_to_pattern(icd10),
    sr_codes = if (length(sr) > 0) unique(sr) else NULL,
    death_icd10 = .catalog_codes_to_pattern(death_icd10),
    opcs4_pattern = .catalog_codes_to_pattern(opcs4),
    cancer_icd10_pattern = .catalog_codes_to_pattern(cancer_icd10),
    cancer_behaviour = if (length(cancer_icd10) > 0) 3L else NULL
  )
}

.catalog_codes_to_pattern <- function(codes) {
  codes <- unique(as.character(codes))
  codes <- codes[!is.na(codes) & nzchar(codes)]
  if (length(codes) == 0) {
    return(NULL)
  }
  paste0("^(", paste(codes, collapse = "|"), ")")
}

.merge_predefined_diseases <- function(curated,
                                       pomegranate,
                                       merge_type = c("intersection", "union")) {
  merge_type <- match.arg(merge_type)
  map <- .curated_pomegranate_definition_map()
  map <- map[names(map) %in% names(curated)]
  map <- lapply(map, function(x) x[x %in% names(pomegranate)])
  map <- map[lengths(map) > 0]

  out <- list()
  for (curated_key in names(map)) {
    cur <- curated[[curated_key]]
    pom_defs <- pomegranate[map[[curated_key]]]
    def <- if (identical(merge_type, "union")) {
      .union_disease_definition_set(cur, pom_defs)
    } else {
      .intersect_disease_definition_set(cur, pom_defs)
    }
    if (!is.null(def)) {
      out[[curated_key]] <- def
    }
  }

  if (identical(merge_type, "union")) {
    unmatched_curated <- setdiff(names(curated), names(map))
    out <- c(out, curated[unmatched_curated])

    matched_pomegranate <- unique(unlist(map, use.names = FALSE))
    unmatched_pomegranate <- setdiff(names(pomegranate), matched_pomegranate)
    out <- c(out, pomegranate[unmatched_pomegranate])
  }

  out
}

.union_disease_definition_set <- function(curated, pomegranate_defs) {
  pom <- do.call(
    combine_disease_definitions,
    c(pomegranate_defs, list(name = curated$name))
  )

  sr <- unique(c(
    as.integer(curated$sr_codes %||% integer()),
    as.integer(pom$sr_codes %||% integer())
  ))
  sr <- sr[!is.na(sr)]

  create_disease_definition(
    name = curated$name,
    icd10_pattern = .union_definition_patterns(
      curated$icd10_pattern,
      pom$icd10_pattern
    ),
    icd9_pattern = .union_definition_patterns(
      curated$icd9_pattern,
      pom$icd9_pattern
    ),
    sr_codes = if (length(sr) > 0) sr else NULL,
    death_icd10 = .union_definition_patterns(
      curated$death_icd10,
      pom$death_icd10
    ),
    opcs4_pattern = .union_definition_patterns(
      curated$opcs4_pattern,
      pom$opcs4_pattern
    ),
    first_occurrence_fields = curated$first_occurrence_fields,
    first_occurrence_source_fields = curated$first_occurrence_source_fields,
    cancer_icd10_pattern = .union_definition_patterns(
      curated$cancer_icd10_pattern,
      pom$cancer_icd10_pattern
    ),
    cancer_histology = curated$cancer_histology,
    cancer_behaviour = curated$cancer_behaviour,
    algo_date_field = curated$algo_date_field,
    algo_source_field = curated$algo_source_field
  )
}

.union_definition_patterns <- function(...) {
  patterns <- unlist(list(...), use.names = FALSE)
  patterns <- unique(as.character(patterns))
  patterns <- patterns[!is.na(patterns) & nzchar(patterns)]
  if (length(patterns) == 0L) {
    return(NULL)
  }
  if (length(patterns) == 1L) {
    return(patterns[[1]])
  }
  paste0("(", paste0("(", patterns, ")", collapse = "|"), ")")
}

.curated_pomegranate_definition_map <- function() {
  list(
    AA = c("aaa", "abdominal_aortic_aneurysm"),
    AAA = c("aaa", "abdominal_aortic_aneurysm"),
    CVD = c("chd_nos", "coronary_heart_disease_not_otherwise_specified"),
    MI = c("myocardial_infarction"),
    HF = c("hf", "heart_failure"),
    Stroke = c("stroke_nos"),
    Ischaemic_Stroke = c("isch_stroke", "ischaemic_stroke"),
    Intracerebral_Haemorrhage = c("intracereb_haem", "intracerebral_haemorrhage"),
    Subarachnoid_Haemorrhage = c("subarach", "subarachnoid_haemorrhage"),
    Hypertension = c("hypertension"),
    Diabetes = c("diabetes_nos", "diabetes", "diabetes_type_i", "diabetes_type_ii"),
    T1DM = c("diabetes_t1", "diabetes_type_i"),
    T2DM = c("diabetes_t2", "diabetes_type_ii"),
    PCOS = c("pcos", "polycystic_ovarian_syndrome"),
    Angina = c("stable_angina", "unstable_angina"),
    Atrial_Fibrillation = c("af", "atrial_fibrillation"),
    Ventricular_Arrhythmia = c("vt", "ventricular_tachycardia"),
    AV_Block = c(
      "av_block_1",
      "av_block_2",
      "av_block_3",
      "atrioventricular_block_complete",
      "atrioventricular_block_first_degree",
      "atrioventricular_block_second_degree"
    ),
    SVT = c("svt", "supraventricular_tachycardia"),
    Asthma = c("asthma"),
    COPD = c("copd", "copd_excl_bronchitis_nos", "copd_excluding_bronchitis"),
    Lung_Cancer = c("pri_lung", "primary_malignancy_of_lung_and_trachea"),
    CKD = c("ckd", "chronic_kidney_disease"),
    ESRD = c("esrd", "end_stage_renal_disease"),
    PAD = c("peripheral_arterial_disease"),
    VTE = c("pe", "vte_ex_pe", "pulmonary_embolism", "venous_thromboembolic_disease_excl_pe"),
    Osteoarthritis = c("oa", "osteoarthritis_excl_spine"),
    Rheumatoid_Arthritis = c("rha", "rheumatoid_arthritis"),
    Parkinsons = c("parkinsons", "parkinson_s_disease"),
    Dementia = c("dementia"),
    Alzheimers_Disease = c("alzheimer", "alzheimer_s_disease"),
    Vascular_Dementia = c("vascular_dementia_vad"),
    Motor_Neurone_Disease = c("mnd", "motor_neuron_disease"),
    Epilepsy = c("epilepsy"),
    Depression = c("depression"),
    Anxiety = c("anxiety", "anxiety_disorders"),
    Stroke_TIA = c("stroke_nos", "tia", "transient_ischaemic_attack"),
    Irritable_Bowel_Syndrome = c("ibs", "irritable_bowel_syndrome"),
    Inflammatory_Bowel_Disease = c("ibd", "inflammatory_bowel_disease_ibd"),
    Diverticular_Disease = c("diverticuli", "diverticular_disease_of_intestine_acute_and_chronic"),
    Alcohol_Use_Disorder = c("alc_problems", "liver_alc", "alcohol_problems", "alcoholic_liver_disease"),
    Substance_Use_Disorder = c("substance_misuse", "other_psychoactive_substance_misuse"),
    Schizophrenia_Bipolar = c("schizo", "bad", "schizophrenia", "bipolar_affective_disorder_and_mania"),
    Migraine = c("migraine"),
    Bronchiectasis = c("bronchiectasis"),
    Multiple_Sclerosis = c("ms", "multiple_sclerosis"),
    Menieres_Disease = c("meniere", "meniere_disease"),
    Psoriasis_Eczema = c("psoriasis"),
    Fracture = c("fracture_of_hip", "fracture_of_wrist"),
    Glaucoma = c("glaucoma"),
    Cataract = c("cataract"),
    AMD = c("macula_degen", "macular_degeneration"),
    Prostate_Disorders = c("bph", "hyperplasia_of_prostate"),
    Breast_Cancer = c("pri_breast", "primary_malignancy_of_breast"),
    Prostate_Cancer = c("pri_prost", "primary_malignancy_of_prostate"),
    Colorectal_Cancer = c("pri_bowel", "primary_malignancy_of_colorectal_and_anus"),
    Erectile_Dysfunction = c("ed", "erectile_dysfunction"),
    Systemic_Lupus_Erythematosus = c("sle", "lupus_erythematosus_local_and_systemic"),
    Uterus_Cancer = c("pri_uterine", "primary_malignancy_of_uterus"),
    Melanoma = c("pri_melanoma", "primary_malignancy_of_malignant_melanoma"),
    Non_Melanoma_Skin_Cancer = c("pri_skin", "primary_malignancy_of_other_skin_and_subcutaneous_tissue"),
    Ovarian_Cancer = c("pri_ovarian", "primary_malignancy_of_ovarian"),
    Oesophageal_Cancer = c("pri_oesoph", "primary_malignancy_of_oesophageal"),
    Stomach_Cancer = c("pri_stomach", "primary_malignancy_of_stomach")
  )
}

.intersect_disease_definition_set <- function(curated, pomegranate_defs) {
  pom <- do.call(
    combine_disease_definitions,
    c(pomegranate_defs, list(name = curated$name))
  )

  icd10 <- .intersect_code_vectors(
    .definition_pattern_codes(curated$icd10_pattern),
    .definition_pattern_codes(pom$icd10_pattern)
  )
  death_icd10 <- .intersect_code_vectors(
    .definition_pattern_codes(curated$death_icd10),
    .definition_pattern_codes(pom$death_icd10)
  )
  opcs4 <- .intersect_code_vectors(
    .definition_pattern_codes(curated$opcs4_pattern),
    .definition_pattern_codes(pom$opcs4_pattern)
  )
  cancer_icd10 <- .intersect_code_vectors(
    .definition_pattern_codes(curated$cancer_icd10_pattern),
    .definition_pattern_codes(pom$cancer_icd10_pattern)
  )
  sr <- intersect(
    as.integer(curated$sr_codes %||% integer()),
    as.integer(pom$sr_codes %||% integer())
  )
  sr <- sr[!is.na(sr)]

  if (
    length(icd10) == 0 &&
      length(death_icd10) == 0 &&
      length(opcs4) == 0 &&
      length(cancer_icd10) == 0 &&
      length(sr) == 0
  ) {
    return(NULL)
  }

  create_disease_definition(
    name = curated$name,
    icd10_pattern = .catalog_codes_to_pattern(icd10),
    sr_codes = if (length(sr) > 0) sr else NULL,
    death_icd10 = .catalog_codes_to_pattern(death_icd10),
    opcs4_pattern = .catalog_codes_to_pattern(opcs4),
    first_occurrence_fields = curated$first_occurrence_fields,
    first_occurrence_source_fields = curated$first_occurrence_source_fields,
    cancer_icd10_pattern = .catalog_codes_to_pattern(cancer_icd10),
    cancer_histology = curated$cancer_histology,
    cancer_behaviour = curated$cancer_behaviour,
    algo_date_field = curated$algo_date_field,
    algo_source_field = curated$algo_source_field
  )
}

.definition_pattern_codes <- function(pattern) {
  if (is.null(pattern) || length(pattern) == 0 || is.na(pattern) || !nzchar(pattern)) {
    return(character())
  }
  pattern <- gsub("\\\\", "", pattern)
  pattern <- gsub("\\^", "", pattern)
  pattern <- gsub("[()]", "", pattern)
  pattern <- gsub("\\.\\*", "", pattern)
  pattern <- gsub("\\$", "", pattern)
  out <- unlist(strsplit(pattern, "\\|"), use.names = FALSE)
  out <- toupper(gsub("[^A-Z0-9]", "", trimws(out)))
  unique(out[nzchar(out)])
}

.intersect_code_vectors <- function(x, y) {
  x <- unique(toupper(gsub("\\.", "", as.character(x))))
  y <- unique(toupper(gsub("\\.", "", as.character(y))))
  x <- x[!is.na(x) & nzchar(x)]
  y <- y[!is.na(y) & nzchar(y)]
  if (length(x) == 0 || length(y) == 0) {
    return(character())
  }

  out <- character()
  for (a in x) {
    for (b in y) {
      if (identical(a, b)) {
        out <- c(out, a)
      } else if (startsWith(a, b)) {
        out <- c(out, a)
      } else if (startsWith(b, a)) {
        out <- c(out, b)
      }
    }
  }
  unique(out)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.filter_disease_definition_list <- function(definitions, disease = NULL) {
  if (is.null(disease)) {
    return(definitions)
  }
  if (!is.character(disease) || anyNA(disease) || any(!nzchar(disease))) {
    stop("`disease` must be NULL or a non-empty character vector.", call. = FALSE)
  }

  keys <- names(definitions)
  display_names <- vapply(definitions, function(x) x$name %||% "", character(1))
  keep <- rep(FALSE, length(definitions))
  for (one in disease) {
    keep <- keep |
      grepl(one, keys, ignore.case = TRUE) |
      grepl(one, display_names, ignore.case = TRUE)
  }
  definitions[keep]
}
