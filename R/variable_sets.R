#' Curated UK Biobank variable sets for extraction
#'
#' @description
#' Returns curated UKB field groups for common analysis domains. These sets are
#' intended for field discovery and RAP extraction, not for automatic
#' preprocessing. Use [preprocess_baseline()] only for variables documented by
#' [get_variable_info()].
#'
#' @param set Optional set name, such as `"clinical_core"`,
#'   `"air_pollution"`, or `"family_history"`. If `NULL`, returns all rows.
#' @param category Optional broad category filter.
#'
#' @return A data.frame with one row per curated variable.
#' @export
#'
#' @examples
#' vars <- get_variable_sets("air_pollution")
#' unique(vars$field_id)
get_variable_sets <- function(set = NULL,
                              category = NULL) {
  out <- .get_variable_set_catalog()

  if (!is.null(set)) {
    if (!is.character(set) || length(set) == 0 || anyNA(set)) {
      stop("`set` must be a non-empty character vector.", call. = FALSE)
    }
    missing <- setdiff(set, unique(out$set))
    if (length(missing) > 0) {
      stop(
        "Unknown variable set(s): ", paste(missing, collapse = ", "),
        ". Use get_variable_sets() to inspect available sets.",
        call. = FALSE
      )
    }
    out <- out[out$set %in% set, , drop = FALSE]
  }

  if (!is.null(category)) {
    if (!is.character(category) || length(category) == 0 || anyNA(category)) {
      stop("`category` must be a non-empty character vector.", call. = FALSE)
    }
    missing <- setdiff(category, unique(out$category))
    if (length(missing) > 0) {
      stop(
        "Unknown category: ", paste(missing, collapse = ", "),
        ". Use unique(get_variable_sets()$category) to inspect categories.",
        call. = FALSE
      )
    }
    out <- out[out$category %in% category, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

#' Get one curated UK Biobank variable set
#'
#' @param set Set name.
#' @param output Output format. `"data.frame"` returns the full manifest,
#'   `"field_id"` returns unique UKB field IDs, and `"ukb_col"` returns UKB
#'   column stems such as `p31` or `p4080_i0_a0`.
#'
#' @return A data.frame or character/integer vector.
#' @export
#'
#' @examples
#' get_variable_set("clinical_core")
#' get_variable_set("air_pollution", output = "field_id")
get_variable_set <- function(set,
                             output = c("data.frame", "field_id", "ukb_col")) {
  output <- match.arg(output)
  out <- get_variable_sets(set = set)

  if (identical(output, "field_id")) {
    return(unique(out$field_id[!is.na(out$field_id)]))
  }
  if (identical(output, "ukb_col")) {
    return(unique(out$ukb_col[!is.na(out$ukb_col) & nzchar(out$ukb_col)]))
  }
  out
}

.get_variable_set_catalog <- function() {
  rows <- list()
  add <- function(set, category, variable, field_id, ukb_col = NULL,
                  label = NULL, role = "raw_field", notes = NA_character_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      set = set,
      category = category,
      variable = variable,
      field_id = as.integer(field_id),
      ukb_col = if (is.null(ukb_col)) paste0("p", field_id) else ukb_col,
      label = if (is.null(label)) variable else label,
      role = role,
      notes = notes,
      stringsAsFactors = FALSE
    )
  }

  # Core time and demographics
  add("time_core", "time", "year_of_birth", 34, "p34", "Year of birth")
  add("time_core", "time", "month_of_birth", 52, "p52", "Month of birth")
  add("time_core", "time", "baseline_date", 53, "p53_i0", "Date of attending assessment centre")
  add("time_core", "time", "age_at_recruitment", 21022, "p21022", "Age at recruitment")
  add("time_core", "time", "lost_to_followup_date", 191, "p191", "Date lost to follow-up")
  add("time_core", "time", "interview_date", 200, "p200", "Date of verbal interview")
  add("time_core", "time", "blood_sample_collection_time", 3166, "p3166", "Time of blood sample collection")
  add("time_core", "time", "death_date", 40000, "p40000", "Date of death")
  add("time_core", "time", "underlying_death_cause_icd10", 40001, "p40001", "Underlying primary cause of death: ICD10")
  add("time_core", "time", "age_at_death", 40007, "p40007", "Age at death")

  add("clinical_core", "demographics", "sex", 31, "p31", "Sex")
  add("clinical_core", "demographics", "age_at_recruitment", 21022, "p21022", "Age at recruitment")
  add("clinical_core", "demographics", "ethnic_background", 21000, "p21000_i0", "Ethnic background")
  add("clinical_core", "socioeconomic", "household_income", 738, "p738_i0", "Average total household income")
  add("clinical_core", "socioeconomic", "qualifications", 6138, "p6138_i0", "Qualifications")
  add("clinical_core", "socioeconomic", "townsend_deprivation_index", 189, "p189", "Townsend deprivation index")
  add("clinical_core", "anthropometrics", "bmi", 21001, "p21001_i0", "Body mass index")
  add("clinical_core", "anthropometrics", "height", 50, "p50_i0", "Standing height")
  add("clinical_core", "anthropometrics", "weight", 21002, "p21002_i0", "Weight")
  add("clinical_core", "blood_pressure", "sbp_auto_1", 4080, "p4080_i0_a0", "Automated systolic blood pressure reading 1")
  add("clinical_core", "blood_pressure", "sbp_auto_2", 4080, "p4080_i0_a1", "Automated systolic blood pressure reading 2")
  add("clinical_core", "blood_pressure", "dbp_auto_1", 4079, "p4079_i0_a0", "Automated diastolic blood pressure reading 1")
  add("clinical_core", "blood_pressure", "dbp_auto_2", 4079, "p4079_i0_a1", "Automated diastolic blood pressure reading 2")
  add("clinical_core", "blood_pressure", "sbp_manual_1", 93, "p93_i0_a0", "Manual systolic blood pressure reading 1")
  add("clinical_core", "blood_pressure", "sbp_manual_2", 93, "p93_i0_a1", "Manual systolic blood pressure reading 2")
  add("clinical_core", "blood_pressure", "dbp_manual_1", 94, "p94_i0_a0", "Manual diastolic blood pressure reading 1")
  add("clinical_core", "blood_pressure", "dbp_manual_2", 94, "p94_i0_a1", "Manual diastolic blood pressure reading 2")

  biomarkers <- list(
    wbc = c(30000, "White blood cell count"),
    rbc = c(30010, "Red blood cell count"),
    haemoglobin = c(30020, "Haemoglobin concentration"),
    haematocrit = c(30030, "Haematocrit percentage"),
    platelet = c(30080, "Platelet count"),
    lymphocyte_count = c(30120, "Lymphocyte count"),
    monocyte_count = c(30130, "Monocyte count"),
    neutrophil_count = c(30140, "Neutrophil count"),
    alanine_aminotransferase = c(30620, "Alanine aminotransferase"),
    albumin = c(30600, "Albumin"),
    alkaline_phosphatase = c(30610, "Alkaline phosphatase"),
    apolipoprotein_a = c(30630, "Apolipoprotein A"),
    apolipoprotein_b = c(30640, "Apolipoprotein B"),
    aspartate_aminotransferase = c(30650, "Aspartate aminotransferase"),
    c_reactive_protein = c(30710, "C-reactive protein"),
    calcium = c(30680, "Calcium"),
    cholesterol = c(30690, "Cholesterol"),
    creatinine = c(30700, "Creatinine"),
    cystatin_c = c(30720, "Cystatin C"),
    direct_bilirubin = c(30660, "Direct bilirubin"),
    gamma_glutamyltransferase = c(30730, "Gamma glutamyltransferase"),
    glucose = c(30740, "Glucose"),
    hba1c = c(30750, "Glycated haemoglobin HbA1c"),
    hdl_cholesterol = c(30760, "HDL cholesterol"),
    ldl_cholesterol = c(30780, "LDL direct"),
    lipoprotein_a = c(30790, "Lipoprotein A"),
    testosterone = c(30850, "Testosterone"),
    total_bilirubin = c(30840, "Total bilirubin"),
    triglycerides = c(30870, "Triglycerides"),
    urate = c(30880, "Urate"),
    urea = c(30670, "Urea")
  )
  for (nm in names(biomarkers)) {
    add("biomarkers_core", "biomarkers", nm, as.integer(biomarkers[[nm]][[1]]), paste0("p", biomarkers[[nm]][[1]], "_i0"), biomarkers[[nm]][[2]])
  }

  add("lifestyle_core", "lifestyle", "smoking_status", 20116, "p20116_i0", "Smoking status")
  add("lifestyle_core", "lifestyle", "alcohol_drinker_status", 20117, "p20117_i0", "Alcohol drinker status")
  add("lifestyle_core", "lifestyle", "sleep_duration", 1160, "p1160_i0", "Sleep duration")
  add("lifestyle_core", "lifestyle", "insomnia", 1200, "p1200_i0", "Sleeplessness / insomnia")
  add("lifestyle_core", "lifestyle", "snoring", 1210, "p1210_i0", "Snoring")
  add("lifestyle_core", "lifestyle", "daytime_dozing", 1220, "p1220_i0", "Daytime dozing / sleeping")
  add("lifestyle_core", "lifestyle", "ipaq_activity_group", 22032, "p22032_i0", "IPAQ activity group")

  # Self-report, medication, procedure fields
  add("self_report_diagnosis", "self_report", "doctor_diagnosed_vascular_conditions", 6150, "p6150_i0", "Vascular/diabetes doctor diagnosis")
  add("self_report_diagnosis", "self_report", "non_cancer_illness", 20002, "p20002_i0_a0", "Non-cancer illness code, self-reported")
  add("self_report_diagnosis", "self_report", "cancer_illness", 20001, "p20001_i0_a0", "Cancer code, self-reported")
  add("self_report_diagnosis", "self_report", "operation_code", 20004, "p20004_i0_a0", "Operation code")
  add("self_report_medication", "medication", "male_medication", 6177, "p6177_i0", "Medication for cholesterol, blood pressure or diabetes")
  add("self_report_medication", "medication", "female_medication", 6153, "p6153_i0", "Medication for cholesterol, blood pressure, diabetes or hormones")
  add("self_report_medication", "medication", "regular_medication_code", 20003, "p20003_i0_a0", "Treatment/medication code")
  add("self_report_medication", "medication", "regular_pain_medication", 6154, "p6154_i0", "Medication for pain relief")
  add("self_report_medication", "medication", "regular_stomach_medication", 10004, "p10004_i0", "Medication for stomach symptoms")
  add("self_report_medication", "medication", "regular_constipation_medication", 10005, "p10005_i0", "Medication for constipation")

  # Algorithmically-defined outcomes
  algo <- list(
    Dementia = 42018, Alzheimers_disease = 42020, Vascular_dementia = 42022,
    Frontotemporal_dementia = 42024, End_stage_renal_disease = 42026,
    Motor_neurone_disease = 42028, Parkinsonism = 42030,
    Parkinsons_disease = 42032, Progressive_supranuclear_palsy = 42034,
    Multiple_system_atrophy = 42036
  )
  for (nm in names(algo)) {
    add("algorithm_outcomes", "algorithm", .ukb_snake(nm), algo[[nm]], paste0("p", algo[[nm]]), nm)
  }

  # Environmental exposure fields
  env_air <- c(
    no2_2010 = 24003, nox_2010 = 24004, pm10_2010 = 24005, pm25_2010 = 24006,
    pm25_absorbance_2010 = 24007, pm25_10_2010 = 24008, traffic_nearest_road = 24009,
    inverse_distance_nearest_road = 24010, traffic_nearest_major_road = 24011,
    inverse_distance_major_road = 24012, traffic_load_major_roads = 24013,
    close_to_major_road = 24014, major_road_length_100m = 24015,
    no2_2005 = 24016, no2_2006 = 24017, no2_2007 = 24018, pm10_2007 = 24019
  )
  for (nm in names(env_air)) add("air_pollution", "environment", nm, env_air[[nm]], paste0("p", env_air[[nm]]), nm)

  water <- c(ca_concentration = 21104, caco3_concentration = 21103, mg_concentration = 21105,
             water_hardness_usgs = 21100, water_hardness_who = 21101, water_survey_year = 21102)
  for (nm in names(water)) add("water_quality", "environment", nm, water[[nm]], paste0("p", water[[nm]]), nm)

  noise <- c(noise_daytime = 24020, noise_evening = 24021, noise_night = 24022,
             noise_16h = 24023, noise_24h = 24024)
  for (nm in names(noise)) add("noise_pollution", "environment", nm, noise[[nm]], paste0("p", noise[[nm]]), nm)

  greenspace_ids <- c(24500, 24501, 24502, 24503, 24504, 24505, 24506, 24507, 24508)
  greenspace_names <- c(
    "greenspace_1000m", "domestic_garden_1000m", "water_1000m",
    "greenspace_300m", "domestic_garden_300m", "water_300m",
    "natural_environment_1000m", "natural_environment_300m", "distance_to_coast"
  )
  for (i in seq_along(greenspace_ids)) {
    add("greenspace", "environment", greenspace_names[[i]], greenspace_ids[[i]], paste0("p", greenspace_ids[[i]], "_i0"), greenspace_names[[i]])
  }

  # Family, sun, electronic device use, genetics.
  family <- c(father_alive = 1797, father_age = 2946, father_age_at_death = 1807,
              mother_alive = 1835, mother_age = 1845, mother_age_at_death = 3526,
              non_accidental_death_family = 4501, father_illnesses = 20107,
              mother_illnesses = 20110, sibling_illnesses = 20111)
  for (nm in names(family)) add("family_history", "family", nm, family[[nm]], paste0("p", family[[nm]], "_i0"), nm)

  sun <- c(time_outdoors_summer = 1050, time_outdoors_winter = 1060, skin_colour = 1717,
           ease_skin_tanning = 1727, childhood_sunburn = 1737, hair_colour = 1747,
           facial_ageing = 1757, sun_uv_protection = 2267, solarium_use = 2277)
  for (nm in names(sun)) add("sun_exposure", "lifestyle", nm, sun[[nm]], paste0("p", sun[[nm]], "_i0"), nm)

  devices <- c(mobile_phone_use_length = 1110, weekly_mobile_phone_use = 1120,
               hands_free_device_use = 1130, mobile_phone_use_change = 1140,
               usual_side_mobile_phone = 1150, computer_games = 2237)
  for (nm in names(devices)) add("electronic_device_use", "lifestyle", nm, devices[[nm]], paste0("p", devices[[nm]], "_i0"), nm)

  for (i in seq_len(20)) {
    add("genetic_pcs_20", "genetics", paste0("genetic_pc_", i), 22009, paste0("p22009_a", i), paste0("Genetic principal component ", i))
  }

  out <- do.call(rbind, rows)
  out[order(out$category, out$set, out$variable), , drop = FALSE]
}

.ukb_snake <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_+|_+$", "", x)
}
