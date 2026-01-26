#' Variable preprocessing functions for UKB baseline characteristics
#' 
#' This module provides flexible variable preprocessing with automatic field ID mapping
#' and standardized transformations for common UKB baseline variables.
#' Supports both predefined mappings and user-defined custom mappings.
#' 
#' @name variable_preprocess
#' @keywords internal
NULL

# Declare global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  "fruit_score", "vegetable_score", "fish_score", "meat_score", "diet_score",
  "p1289_i0", "p1299_i0", "p1309_i0", "p1319_i0", "p1329_i0", "p1339_i0"
))

# Variable Mapping Table

#' Get default variable to UKB field ID mapping
#' @return A named list with variable mappings
#' @keywords internal
.get_variable_mapping <- function() {
  list(
    # Demographics 
    sex = list(
      field_id = 31,
      ukb_col = "p31",
      description = "Sex (0=Female, 1=Male)"
    ),
    age = list(
      field_id = 21022,
      ukb_col = "p21022",
      description = "Age at recruitment (years)"
    ),
    ethnicity = list(
      field_id = 21000,
      ukb_col = "p21000_i0",
      description = "Ethnic background (White vs Others)"
    ),
    
    # Anthropometrics
    bmi = list(
      field_id = 21001,
      ukb_col = "p21001_i0",
      description = "Body mass index (kg/m^2)"
    ),
    height = list(
      field_id = 50,
      ukb_col = "p50_i0",
      description = "Standing height (cm)"
    ),
    weight = list(
      field_id = 21002,
      ukb_col = "p21002_i0",
      description = "Weight (kg)"
    ),
    
    # Lifestyle
    smoking = list(
      field_id = 20116,
      ukb_col = "p20116_i0",
      description = "Smoking status (0=Never, 1=Previous, 2=Current)"
    ),
    drinking = list(
      field_id = 20117,
      ukb_col = "p20117_i0",
      description = "Alcohol drinker status"
    ),
    sleep_duration = list(
      field_id = 1160,
      ukb_col = "p1160_i0",
      description = "Sleep duration (hours/day)"
    ),
    exercise_intensity = list(
      field_id = 22032,
      ukb_col = "p22032_i0",
      description = "Physical activity intensity (MET)"
    ),
    
    # Socioeconomic
    education = list(
      field_id = 6138,
      ukb_col = "p6138_i0",
      description = "Qualifications (recoded: 1=Low, 2=Medium, 3=High)"
    ),
    income = list(
      field_id = 738,
      ukb_col = "p738_i0",
      description = "Average total household income"
    ),
    
    # Blood Pressure
    sbp_auto_1 = list(field_id = 93, ukb_col = "p93_i0_a0", description = "SBP automated reading 1"),
    sbp_auto_2 = list(field_id = 93, ukb_col = "p93_i0_a1", description = "SBP automated reading 2"),
    dbp_auto_1 = list(field_id = 94, ukb_col = "p94_i0_a0", description = "DBP automated reading 1"),
    dbp_auto_2 = list(field_id = 94, ukb_col = "p94_i0_a1", description = "DBP automated reading 2"),
    sbp_manual_1 = list(field_id = 4080, ukb_col = "p4080_i0_a0", description = "SBP manual reading 1"),
    sbp_manual_2 = list(field_id = 4080, ukb_col = "p4080_i0_a1", description = "SBP manual reading 2"),
    dbp_manual_1 = list(field_id = 4079, ukb_col = "p4079_i0_a0", description = "DBP manual reading 1"),
    dbp_manual_2 = list(field_id = 4079, ukb_col = "p4079_i0_a1", description = "DBP manual reading 2"),
    
    # Medications 
    medication_male = list(
      field_id = 6177,
      ukb_col = "p6177_i0",
      description = "Medication for cholesterol, BP, diabetes (male)"
    ),
    medication_female = list(
      field_id = 6153,
      ukb_col = "p6153_i0",
      description = "Medication for cholesterol, BP, diabetes, HRT (female)"
    ),
    
    # Biomarkers
    triglycerides = list(field_id = 30870, ukb_col = "p30870_i0", description = "Triglycerides (mmol/L)"),
    ldl = list(field_id = 30780, ukb_col = "p30780_i0", description = "LDL cholesterol (mmol/L)"),
    hdl = list(field_id = 30760, ukb_col = "p30760_i0", description = "HDL cholesterol (mmol/L)"),
    hba1c = list(field_id = 30750, ukb_col = "p30750_i0", description = "HbA1c (mmol/mol)"),
    glucose = list(field_id = 30740, ukb_col = "p30740_i0", description = "Glucose (mmol/L)"),
    
    # Air Pollution
    no2_2005 = list(field_id = 24016, ukb_col = "p24016", description = "NO2 air pollution 2005"),
    no2_2006 = list(field_id = 24018, ukb_col = "p24018", description = "NO2 air pollution 2006"),
    no2_2007 = list(field_id = 24017, ukb_col = "p24017", description = "NO2 air pollution 2007"),
    no2_2010 = list(field_id = 24003, ukb_col = "p24003", description = "NO2 air pollution 2010"),
    nox = list(field_id = 24004, ukb_col = "p24004", description = "NOx air pollution"),
    pm25 = list(field_id = 24006, ukb_col = "p24006", description = "PM2.5 air pollution"),
    pm10_2007 = list(field_id = 24019, ukb_col = "p24019", description = "PM10 air pollution 2007"),
    pm10_2010 = list(field_id = 24005, ukb_col = "p24005", description = "PM10 air pollution 2010"),
    
    # Diet
    cooked_vegetable = list(field_id = 1289, ukb_col = "p1289_i0", description = "Cooked vegetable intake"),
    salad_vegetable = list(field_id = 1299, ukb_col = "p1299_i0", description = "Raw salad/vegetable intake"),
    fresh_fruit = list(field_id = 1309, ukb_col = "p1309_i0", description = "Fresh fruit intake"),
    dried_fruit = list(field_id = 1319, ukb_col = "p1319_i0", description = "Dried fruit intake"),
    oily_fish = list(field_id = 1329, ukb_col = "p1329_i0", description = "Oily fish intake"),
    non_oily_fish = list(field_id = 1339, ukb_col = "p1339_i0", description = "Non-oily fish intake"),
    processed_meat = list(field_id = 1349, ukb_col = "p1349_i0", description = "Processed meat intake"),
    poultry = list(field_id = 1359, ukb_col = "p1359_i0", description = "Poultry intake"),
    beef = list(field_id = 1369, ukb_col = "p1369_i0", description = "Beef intake"),
    lamb = list(field_id = 1379, ukb_col = "p1379_i0", description = "Lamb/mutton intake"),
    pork = list(field_id = 1389, ukb_col = "p1389_i0", description = "Pork intake"),
    milk = list(field_id = 1418, ukb_col = "p1418_i0", description = "Milk type used"),
    bread_intake = list(field_id = 1438, ukb_col = "p1438_i0", description = "Bread intake"),
    bread_type = list(field_id = 1448, ukb_col = "p1448_i0", description = "Bread type"),
    cereal_intake = list(field_id = 1458, ukb_col = "p1458_i0", description = "Cereal intake"),
    cereal_type = list(field_id = 1468, ukb_col = "p1468_i0", description = "Cereal type"),
    never_eat_eggs_dairy_wheat = list(field_id = 6144, ukb_col = "p6144_i0", description = "Never eat eggs/dairy/wheat/sugar")
  )
}

# Core Preprocessing Function

#' Preprocess UKB baseline variables
#' @description A unified function to preprocess UKB baseline characteristics with automatic
#' field mapping and standardized transformations.
#' @param df A data.table or data.frame containing UKB data from rap platform export.
#' @param variables Character vector of variable names to process. 
#'   Use `get_variable_info()` to see available variables.
#' @param custom_mapping Optional named list for user-defined variable mappings.
#'   Each element should have: `ukb_col` (required), `description` (optional).
#'   Example: `list(my_var = list(ukb_col = "p12345_i0", description = "My custom var"))`
#' @param missing_action Character. How to handle missing values:
#'   - "keep": Keep as NA (default)
#'   - "drop": Remove rows with any missing values in processed variables
#' @param invalid_codes Numeric vector of UKB codes to treat as missing.
#'   Default: c(-1, -3) which are "Prefer not to answer" and "Do not know"
#' @importFrom data.table := data.table copy is.data.table 
#' @importFrom stats complete.cases 
#' @return A data.table with original data plus processed variable columns
#' @examples
#' \dontrun{
#' # Process predefined variables
#' result <- preprocess_baseline(ukb_data, 
#'                               variables = c("sex", "age", "ethnicity", "bmi"))
#' 
#' # Use custom variable mapping
#' custom <- list(
#'   my_biomarker = list(ukb_col = "p30000_i0", description = "Custom biomarker")
#' )
#' result <- preprocess_baseline(ukb_data, 
#'                               variables = c("sex", "my_biomarker"),
#'                               custom_mapping = custom)
#' 
#' # Drop rows with missing values
#' result <- preprocess_baseline(ukb_data, 
#'                               variables = c("sex", "age"),
#'                               missing_action = "drop")
#' }
#' 
#' @export
preprocess_baseline <- function(df, 
                                variables,
                                custom_mapping = NULL,
                                missing_action = c("keep", "drop"),
                                invalid_codes = c(-1, -3)) {
  
  missing_action <- match.arg(missing_action)

  # check input
  if (!is.data.frame(df)) {
    stop("df must be a data.frame or data.table")
  }

  # use data.table for further processing (faster)
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  
  # Merge default and custom mappings
  mapping <- .get_variable_mapping()
  if (!is.null(custom_mapping)) {
    if (!is.list(custom_mapping)) {
      stop("custom_mapping must be a named list")
    }
    # Validate custom mapping structure
    for (var_name in names(custom_mapping)) {
      if (!"ukb_col" %in% names(custom_mapping[[var_name]])) {
        stop(sprintf("Custom mapping for '%s' must include 'ukb_col'", var_name))
      }
    }
    mapping <- c(mapping, custom_mapping)
  }
  
  # Validate variables
  invalid_vars <- setdiff(variables, names(mapping))
  if (length(invalid_vars) > 0) {
    stop("Unknown variables: ", paste(invalid_vars, collapse = ", "),
         "\nUse get_variable_info() to see available variables, or provide custom_mapping.")
  }
  
  result_df <- data.table::copy(df)
  processed_vars <- character(0)
  
  # Process each variable
  for (var in variables) {
    var_info <- mapping[[var]]
    ukb_col <- var_info$ukb_col
    
    # Check if column exists
    if (!ukb_col %in% names(result_df)) {
      warning(sprintf("Column '%s' not found for variable '%s'. Skipping.", ukb_col, var))
      next
    }
    
    # Get processor function
    processor <- .get_processor(var)
    
    if (!is.null(processor)) {
      # Use predefined processor
      result_df[, (var) := processor(.SD[[ukb_col]], invalid_codes)]
    } else {
      # Default: just copy and handle invalid codes for numeric
      raw_col <- result_df[[ukb_col]]
      if (is.numeric(raw_col)) {
        raw_col[raw_col %in% invalid_codes] <- NA
      }
      result_df[, (var) := raw_col]
    }
    
    processed_vars <- c(processed_vars, var)
  }
  
  # Handle missing values
  if (missing_action == "drop" && length(processed_vars) > 0) {
    n_before <- nrow(result_df)
    complete_rows <- stats::complete.cases(result_df[, processed_vars, with = FALSE])
    result_df <- result_df[complete_rows]
    n_dropped <- n_before - nrow(result_df)
    if (n_dropped > 0) {
      message(sprintf("Dropped %d rows with missing values (%.1f%%)", 
                     n_dropped, 100 * n_dropped / n_before))
    }
  }
  
  return(result_df[])
}


#' Get processor function for a variable
#' @keywords internal
.get_processor <- function(var_name) {
  processors <- list(
    # Demographics
    sex = .process_sex,
    age = .process_age,
    ethnicity = .process_ethnicity,
    
    # Anthropometrics
    bmi = .process_bmi,
    height = .process_height,
    weight = .process_weight,
    
    # Lifestyle
    smoking = .process_smoking,
    drinking = .process_drinking,
    sleep_duration = .process_sleep,
    
    # Socioeconomic
    education = .process_education,
    income = .process_income
  )
  
  processors[[var_name]]
}


# Variable-Specific Processors
#' @keywords internal
.process_sex <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.integer(x)
  x[x %in% invalid_codes] <- NA
  factor(x, levels = c(0, 1), labels = c("Female", "Male"))
}

#' @keywords internal
.process_age <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.numeric(x)
  x[x %in% invalid_codes] <- NA
  x
}

#' @keywords internal
.process_ethnicity <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.numeric(x)
  # White = 1, 1001, 1002, 1003; Others = everything else
  result <- ifelse(x %in% c(1, 1001, 1002, 1003), 1L,
                  ifelse(x %in% invalid_codes | is.na(x), NA_integer_, 2L))
  factor(result, levels = c(1, 2), labels = c("White", "Others"))
}

#' @keywords internal
.process_bmi <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.numeric(x)
  x[x %in% invalid_codes] <- NA
  #x[x < 12 | x > 60] <- NA
  x
}

#' @keywords internal
.process_height <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.numeric(x)
  x[x %in% invalid_codes] <- NA
  #x[x < 120 | x > 220] <- NA
  x
}

#' @keywords internal
.process_weight <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.numeric(x)
  x[x %in% invalid_codes] <- NA
  #x[x < 30 | x > 300] <- NA
  x
}

#' @keywords internal
.process_smoking <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.integer(x)
  x[x %in% invalid_codes] <- NA
  factor(x, levels = c(0, 1, 2), labels = c("Never", "Previous", "Current"))
}

#' @keywords internal
.process_drinking <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.integer(x)
  x[x %in% invalid_codes] <- NA
  factor(x, levels = c(0, 1, 2), labels = c("Never", "Previous", "Current"))
}

#' @keywords internal
.process_sleep <- function(x, invalid_codes = c(-1, -3), qc_range = c(1, 24)) {
  x <- as.numeric(x)
  x[x %in% invalid_codes] <- NA
  x[x < qc_range[1] | x > qc_range[2]] <- NA
  x
}

#' @keywords internal
.process_education <- function(x, invalid_codes = c(-1, -3, -7)) {
  x <- as.character(x)
  result <- rep(NA_integer_, length(x))
  
  # -7 = None of the above -> Low education
  result[x == "[-7]" | x == "-7"] <- 1L
  
  # Contains 1, 5, or 6 -> High education (College/University, NVQ, Professional)
  result[grepl("\\b[156]\\b", x)] <- 3L
  
  # Contains 2, 3, or 4 -> Medium education (A levels, O levels, CSEs)
  result[grepl("\\b[234]\\b", x) & is.na(result)] <- 2L
  
  # Invalid codes
  result[x %in% c("[-3]", "[-1]", "-3", "-1", "", NA)] <- NA_integer_
  
  factor(result, levels = c(1, 2, 3), labels = c("Low", "Medium", "High"))
}

#' @keywords internal
.process_income <- function(x, invalid_codes = c(-1, -3)) {
  x <- as.integer(x)
  x[x %in% invalid_codes] <- NA
  x
}


# Composite Variable Functions

#' Calculate blood pressure from multiple readings
#' 
#' Combines manual and automated BP readings, using manual as primary 
#' and automated as fallback. Returns mean of two readings.
#' 
#' @param df A data.table containing BP columns
#' @param type Character: "sbp" or "dbp"
#' @param prefer Character: "manual" (default) or "auto"
#' @importFrom data.table := data.table copy is.data.table as.data.table
#' 
#' @return A data.table with calculated `sbp` or `dbp` column added
#' 
#' @examples
#' \dontrun{
#' result <- calculate_blood_pressure(ukb_data, type = "sbp")
#' result <- calculate_blood_pressure(result, type = "dbp")
#' }
#' 
#' @export
calculate_blood_pressure <- function(df, type = c("sbp", "dbp"), prefer = "manual") {
  type <- match.arg(type)
  
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  result_df <- data.table::copy(df)
  
  if (type == "sbp") {
    manual_cols <- c("p4080_i0_a0", "p4080_i0_a1")
    auto_cols <- c("p93_i0_a0", "p93_i0_a1")
    out_col <- "sbp"
  } else {
    manual_cols <- c("p4079_i0_a0", "p4079_i0_a1")
    auto_cols <- c("p94_i0_a0", "p94_i0_a1")
    out_col <- "dbp"
  }
  
  # Check available columns
  has_manual <- all(manual_cols %in% names(result_df))
  has_auto <- all(auto_cols %in% names(result_df))
  
  if (!has_manual && !has_auto) {
    warning(sprintf("No %s columns found. Skipping.", toupper(type)))
    return(result_df)
  }
  
  if (has_manual && has_auto) {
    # Fill manual with auto where missing
    m1 <- result_df[[manual_cols[1]]]
    m2 <- result_df[[manual_cols[2]]]
    a1 <- result_df[[auto_cols[1]]]
    a2 <- result_df[[auto_cols[2]]]
    
    v1 <- ifelse(is.na(m1), a1, m1)
    v2 <- ifelse(is.na(m2), a2, m2)
  } else if (has_manual) {
    v1 <- result_df[[manual_cols[1]]]
    v2 <- result_df[[manual_cols[2]]]
  } else {
    v1 <- result_df[[auto_cols[1]]]
    v2 <- result_df[[auto_cols[2]]]
  }
  
  # Calculate mean
  result_df[, (out_col) := rowMeans(cbind(v1, v2), na.rm = TRUE)]
  result_df[is.nan(get(out_col)), (out_col) := NA_real_]
  
  return(result_df[])
}

#' Extract medication use from UKB drug fields
#' 
#' Processes medication fields (6177 for male, 6153 for female) to extract
#' specific medication categories.
#' 
#' @param df A data.table containing medication columns (p6177_i0, p6153_i0)
#' @param medications Character vector of medications to extract.
#'   Available: "cholesterol", "blood_pressure", "insulin"
#' @importFrom data.table := data.table copy is.data.table fifelse as.data.table
#' 
#' @return A data.table with binary medication columns added (1=Yes, 0=No, NA=Missing)
#' 
#' @examples
#' \dontrun{
#' result <- extract_medications(ukb_data, c("cholesterol", "blood_pressure"))
#' }
#' 
#' @export
extract_medications <- function(df, medications = c("cholesterol", "blood_pressure", "insulin")) {
  
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  result_df <- data.table::copy(df)
  
  # Check columns
  drug_cols <- c("p6177_i0", "p6153_i0")
  has_cols <- drug_cols %in% names(result_df)
  
  if (!any(has_cols)) {
    warning("Medication columns (p6177_i0, p6153_i0) not found.")
    return(result_df)
  }
  
  # Medication code mapping: 1=Cholesterol, 2=Blood pressure, 3=Insulin
  med_codes <- list(
    cholesterol = "1",
    blood_pressure = "2", 
    insulin = "3"
  )
  
  invalid_patterns <- c("[-1]", "[-3]")
  
  for (med in medications) {
    if (!med %in% names(med_codes)) {
      warning(sprintf("Unknown medication: %s", med))
      next
    }
    
    code <- med_codes[[med]]
    out_col <- paste0("med_", med)
    
    # Get drug values
    drug1 <- if ("p6177_i0" %in% names(result_df)) as.character(result_df$p6177_i0) else rep("", nrow(result_df))
    drug2 <- if ("p6153_i0" %in% names(result_df)) as.character(result_df$p6153_i0) else rep("", nrow(result_df))
    
    # Process: NA if invalid/missing, 1 if contains code, 0 otherwise
    result_df[, (out_col) := {
      is_invalid <- (drug1 %in% invalid_patterns) | (drug2 %in% invalid_patterns) |
                   (is.na(drug1) & is.na(drug2)) |
                   (drug1 == "" & drug2 == "")
      
      has_med <- grepl(code, drug1, fixed = TRUE) | grepl(code, drug2, fixed = TRUE)
      
      fifelse(is_invalid, NA_integer_, fifelse(has_med, 1L, 0L))
    }]
  }
  
  return(result_df[])
}

#' Calculate air pollution exposure averages
#' 
#' Computes averaged air pollution exposures from multiple time points.
#' 
#' @param df A data.table containing air pollution columns
#' @param pollutants Character vector of pollutants to calculate.
#'   Available: "NO2", "PM10", "PM2.5", "NOx"
#' @importFrom data.table := data.table copy is.data.table as.data.table fifelse
#' @return A data.table with averaged pollution columns
#' 
#' @examples
#' \dontrun{
#' result <- calculate_air_pollution(ukb_data, c("NO2", "PM2.5"))
#' }
#' 
#' @export
calculate_air_pollution <- function(df, pollutants = c("NO2", "PM10", "PM2.5", "NOx")) {
  
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  result_df <- data.table::copy(df)
  
  # Pollutant column mappings
  pollutant_cols <- list(
    NO2 = c("p24016", "p24018", "p24017", "p24003"),  # 2005, 2006, 2007, 2010
    PM10 = c("p24019", "p24005"),  # 2007, 2010
    `PM2.5` = "p24006",
    NOx = "p24004"
  )
  
  for (poll in pollutants) {
    if (!poll %in% names(pollutant_cols)) {
      warning(sprintf("Unknown pollutant: %s", poll))
      next
    }
    
    cols <- pollutant_cols[[poll]]
    available_cols <- cols[cols %in% names(result_df)]
    
    if (length(available_cols) == 0) {
      warning(sprintf("No columns found for %s", poll))
      next
    }
    
    out_col <- tolower(gsub("\\.", "", poll))  # NO2, pm10, pm25, nox
    
    if (length(available_cols) == 1) {
      result_df[, (out_col) := .SD[[available_cols]]]
    } else {
      result_df[, (out_col) := rowMeans(.SD, na.rm = TRUE), .SDcols = available_cols]
      result_df[is.nan(get(out_col)), (out_col) := NA_real_]
    }
  }
  
  return(result_df[])
}

#' Calculate diet score
#' 
#' Computes a simplified healthy diet score based on food frequency questionnaire.
#' 
#' @param df A data.table containing diet-related columns
#' @param components Character vector of diet components to include.
#'   Available: "fruit", "vegetable", "fish", "meat", "cereal", "milk"
#' @param na_handling Character: "strict" (NA if any component missing) or 
#'   "partial" (calculate from available components, NA only if insufficient data)
#' @importFrom data.table := data.table copy is.data.table as.data.table fifelse
#' @return A data.table with diet_score column (0-7 scale)
#' 
#' @export
calculate_diet_score <- function(df, 
                                 components = c("fruit", "vegetable", "fish", "meat", "cereal", "milk"),
                                 na_handling = c("strict", "partial")) {
  
  na_handling <- match.arg(na_handling)
  
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  result_df <- data.table::copy(df)
  
  # Calculate component scores
  
  # Fruit score: (fresh_fruit/1 + dried_fruit/5) >= 3
  if ("fruit" %in% components) {
    if (all(c("p1309_i0", "p1319_i0") %in% names(result_df))) {
      result_df[, fruit_score := {
        fresh <- as.numeric(p1309_i0)
        dried <- as.numeric(p1319_i0)
        fresh[fresh %in% c(-1, -3)] <- NA
        dried[dried %in% c(-1, -3)] <- NA
        fifelse(is.na(fresh) | is.na(dried), NA_integer_,
               fifelse((fresh + dried/5) >= 3, 1L, 0L))
      }]
    }
  }
  
  # Vegetable score: (cooked + salad)/3 >= 3
  if ("vegetable" %in% components) {
    if (all(c("p1289_i0", "p1299_i0") %in% names(result_df))) {
      result_df[, vegetable_score := {
        cooked <- as.numeric(p1289_i0)
        salad <- as.numeric(p1299_i0)
        cooked[cooked %in% c(-1, -3)] <- NA
        salad[salad %in% c(-1, -3)] <- NA
        fifelse(is.na(cooked) | is.na(salad), NA_integer_,
               fifelse((cooked + salad)/3 >= 3, 1L, 0L))
      }]
    }
  }
  
  # Fish score: oily fish >= 2 OR non-oily fish >= 2
  if ("fish" %in% components) {
    if (all(c("p1329_i0", "p1339_i0") %in% names(result_df))) {
      result_df[, fish_score := {
        oily <- as.numeric(p1329_i0)
        non_oily <- as.numeric(p1339_i0)
        oily[oily %in% c(-1, -3)] <- NA
        non_oily[non_oily %in% c(-1, -3)] <- NA
        fifelse(is.na(oily) & is.na(non_oily), NA_integer_,
               fifelse(oily >= 2 | non_oily >= 2, 1L, 0L))
      }]
    }
  }
  
  # Meat score: total meat <= 4 (lower is healthier)
  if ("meat" %in% components) {
    meat_cols <- c("p1349_i0", "p1359_i0", "p1369_i0", "p1379_i0", "p1389_i0")
    available_meat_cols <- meat_cols[meat_cols %in% names(result_df)]
    if (length(available_meat_cols) > 0) {
      result_df[, meat_score := {
        vals <- lapply(available_meat_cols, function(col) {
          v <- as.numeric(.SD[[col]])
          v[v %in% c(-1, -3)] <- NA
          v
        })
        total_meat <- Reduce(`+`, vals)
        fifelse(is.na(total_meat), NA_integer_,
               fifelse(total_meat <= 4, 1L, 0L))
      }, .SDcols = available_meat_cols]
    }
  }
  
  # Calculate total score
  score_cols <- grep("_score$", names(result_df), value = TRUE)
  
  if (length(score_cols) > 0) {
    if (na_handling == "strict") {
      result_df[, diet_score := {
        vals <- as.matrix(.SD)
        has_na <- rowSums(is.na(vals)) > 0
        total <- rowSums(vals, na.rm = TRUE)
        fifelse(has_na, NA_integer_, as.integer(total))
      }, .SDcols = score_cols]
    } else {
      # Partial: calculate what's available, NA only if >half missing
      result_df[, diet_score := {
        vals <- as.matrix(.SD)
        n_missing <- rowSums(is.na(vals))
        n_total <- ncol(vals)
        total <- rowSums(vals, na.rm = TRUE)
        fifelse(n_missing > n_total/2, NA_integer_, as.integer(total))
      }, .SDcols = score_cols]
    }
    
    # Clean up intermediate columns
    result_df[, (score_cols) := NULL]
  }
  
  return(result_df[])
}

# Utility Functions

#' Get information about available variables
#' 
#' Returns a data.frame describing all predefined variables available for preprocessing.
#' 
#' @param category Character. Filter by category:
#'   - "all": All variables (default)
#'   - "demographics", "anthropometrics", "lifestyle", "socioeconomic",
#'     "blood_pressure", "medications", "biomarkers", "pollution", "diet"
#' @importFrom data.table data.table
#' @return A data.frame with variable information
#' 
#' @examples
#' \dontrun{
#' # See all available variables
#' get_variable_info()
#' 
#' # Filter by category
#' get_variable_info("demographics")
#' get_variable_info("lifestyle")
#' }
#' 
#' @export
get_variable_info <- function(category = "all") {
  mapping <- .get_variable_mapping()
  
  # Category definitions
  categories <- list(
    demographics = c("sex", "age", "ethnicity"),
    anthropometrics = c("bmi", "height", "weight"),
    lifestyle = c("smoking", "drinking", "sleep_duration", "exercise_intensity"),
    socioeconomic = c("education", "income"),
    blood_pressure = c("sbp_auto_1", "sbp_auto_2", "dbp_auto_1", "dbp_auto_2",
                       "sbp_manual_1", "sbp_manual_2", "dbp_manual_1", "dbp_manual_2"),
    medications = c("medication_male", "medication_female"),
    biomarkers = c("triglycerides", "ldl", "hdl", "hba1c", "glucose"),
    pollution = c("no2_2005", "no2_2006", "no2_2007", "no2_2010", "nox", "pm25", "pm10_2007", "pm10_2010"),
    diet = c("cooked_vegetable", "salad_vegetable", "fresh_fruit", "dried_fruit",
             "oily_fish", "non_oily_fish", "processed_meat", "poultry", "beef",
             "lamb", "pork", "milk", "bread_intake", "bread_type", "cereal_intake",
             "cereal_type", "never_eat_eggs_dairy_wheat")
  )
  
  if (category == "all") {
    selected_vars <- names(mapping)
    var_categories <- sapply(names(mapping), function(v) {
      for (cat in names(categories)) {
        if (v %in% categories[[cat]]) return(cat)
      }
      return("other")
    })
  } else if (category %in% names(categories)) {
    selected_vars <- categories[[category]]
    var_categories <- rep(category, length(selected_vars))
  } else {
    stop("Invalid category. Available: all, ", paste(names(categories), collapse = ", "))
  }
  
  # Build info table
  info_list <- lapply(seq_along(selected_vars), function(i) {
    var <- selected_vars[i]
    if (var %in% names(mapping)) {
      var_info <- mapping[[var]]
      data.frame(
        variable = var,
        category = var_categories[i],
        field_id = var_info$field_id,
        ukb_column = var_info$ukb_col,
        description = var_info$description,
        stringsAsFactors = FALSE
      )
    }
  })
  
  do.call(rbind, info_list)
}
