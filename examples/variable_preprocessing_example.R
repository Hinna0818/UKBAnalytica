# Variable Preprocessing Example
# This script demonstrates how to use the variable preprocessing functions

library(UKBAnalytica)
library(data.table)

# Create sample UKB-like data for demonstration
set.seed(123)
n <- 1000

demo_ukb <- data.table(
  eid = 1:n,
  
  # Demographics  
  p31 = sample(c(0, 1), n, replace = TRUE),  # Sex: 0=Female, 1=Male
  p21022 = rnorm(n, mean = 55, sd = 10),     # Age at recruitment
  p21000_i0 = sample(c(1001, 1002, 3001, 4001, -1, -3), n, replace = TRUE, 
                     prob = c(0.7, 0.1, 0.1, 0.05, 0.03, 0.02)),  # Ethnicity
  
  # Anthropometrics
  p21001_i0 = rnorm(n, mean = 27, sd = 5),    # BMI
  p50_i0 = rnorm(n, mean = 170, sd = 10),     # Height (cm)
  p21002_i0 = rnorm(n, mean = 78, sd = 15),   # Weight (kg)
  
  # Lifestyle
  p20116_i0 = sample(c(0, 1, 2, -1, -3), n, replace = TRUE, 
                     prob = c(0.45, 0.30, 0.15, 0.05, 0.05)),  # Smoking status
  p20117_i0 = sample(c(0, 1, 2, -1), n, replace = TRUE),  # Drinking status
  p1160_i0 = rnorm(n, mean = 7, sd = 1.5),  # Sleep duration
  
  # Education  
  p6138_i0 = sample(c("[1]", "[2, 3]", "[4]", "[-7]", "[-3]", NA), n, replace = TRUE,
                    prob = c(0.2, 0.3, 0.2, 0.1, 0.1, 0.1)),
  
  # Blood pressure (manual and automated)
  p4080_i0_a0 = rnorm(n, mean = 130, sd = 15),  # SBP manual 1
  p4080_i0_a1 = rnorm(n, mean = 128, sd = 15),  # SBP manual 2
  p4079_i0_a0 = rnorm(n, mean = 80, sd = 10),   # DBP manual 1
  p4079_i0_a1 = rnorm(n, mean = 78, sd = 10),   # DBP manual 2
  p93_i0_a0 = rnorm(n, mean = 132, sd = 15),    # SBP auto 1
  p93_i0_a1 = rnorm(n, mean = 130, sd = 15),    # SBP auto 2
  p94_i0_a0 = rnorm(n, mean = 82, sd = 10),     # DBP auto 1
  p94_i0_a1 = rnorm(n, mean = 80, sd = 10),     # DBP auto 2
  
  # Medications
  p6177_i0 = sample(c("[1]", "[2]", "[1, 2]", "[3]", "[-1]", NA), n, replace = TRUE),
  p6153_i0 = sample(c("[1]", "[2]", "[1, 2, 3]", "[-3]", NA), n, replace = TRUE),
  
  # Air pollution
  p24016 = rnorm(n, mean = 25, sd = 5),   # NO2 2005
  p24018 = rnorm(n, mean = 24, sd = 5),   # NO2 2006
  p24017 = rnorm(n, mean = 23, sd = 5),   # NO2 2007
  p24003 = rnorm(n, mean = 22, sd = 5),   # NO2 2010
  p24006 = rnorm(n, mean = 12, sd = 3),   # PM2.5
  p24019 = rnorm(n, mean = 20, sd = 4),   # PM10 2007
  p24005 = rnorm(n, mean = 18, sd = 4)    # PM10 2010
)

# Add some missing values for realistic testing
demo_ukb[sample(n, 50), p31 := NA]
demo_ukb[sample(n, 100), p21022 := NA]
demo_ukb[sample(n, 80), p4080_i0_a0 := NA]  # Missing manual BP

cat("=== Variable Information ===\n")
print(get_variable_info("demographics"))
print(get_variable_info("lifestyle"))

cat("\n=== 1. Basic Preprocessing ===\n")
result1 <- preprocess_baseline(demo_ukb, 
                               variables = c("sex", "age", "ethnicity", "bmi", "smoking"))
cat("Processed columns:", paste(setdiff(names(result1), names(demo_ukb)), collapse = ", "), "\n")
cat("\nSex distribution:\n")
print(table(result1$sex, useNA = "ifany"))
cat("\nEthnicity distribution:\n")
print(table(result1$ethnicity, useNA = "ifany"))
cat("\nSmoking status:\n")
print(table(result1$smoking, useNA = "ifany"))

cat("\n=== 2. Education Processing ===\n")
result2 <- preprocess_baseline(demo_ukb, variables = c("education"))
cat("Education levels:\n")
print(table(result2$education, useNA = "ifany"))

cat("\n=== 3. Custom Variable Mapping ===\n")
# User can define their own variable to column mapping
custom_mapping <- list(
  my_biomarker = list(
    ukb_col = "p21022",  # Use age column as example
    description = "My custom biomarker"
  ),
  pollution_index = list(
    ukb_col = "p24006",
    description = "PM2.5 as pollution index"
  )
)

result3 <- preprocess_baseline(demo_ukb, 
                               variables = c("sex", "my_biomarker", "pollution_index"),
                               custom_mapping = custom_mapping)
cat("Custom variables created:", 
    paste(c("my_biomarker", "pollution_index") %in% names(result3), collapse = ", "), "\n")

cat("\n=== 4. Blood Pressure Calculation ===\n")
result4 <- calculate_blood_pressure(demo_ukb, type = "sbp")
result4 <- calculate_blood_pressure(result4, type = "dbp")
cat("SBP summary:\n")
print(summary(result4$sbp))
cat("\nDBP summary:\n")
print(summary(result4$dbp))

cat("\n=== 5. Medication Extraction ===\n")
result5 <- extract_medications(demo_ukb, c("cholesterol", "blood_pressure", "insulin"))
cat("Medication columns created:\n")
for (col in c("med_cholesterol", "med_blood_pressure", "med_insulin")) {
  cat(sprintf("  %s: Yes=%d, No=%d, NA=%d\n", col,
              sum(result5[[col]] == 1, na.rm = TRUE),
              sum(result5[[col]] == 0, na.rm = TRUE),
              sum(is.na(result5[[col]]))))
}

cat("\n=== 6. Air Pollution Calculation ===\n")
result6 <- calculate_air_pollution(demo_ukb, c("NO2", "PM2.5", "PM10"))
cat("NO2 (averaged):", round(mean(result6$no2, na.rm = TRUE), 2), "\n")
cat("PM2.5:", round(mean(result6$pm25, na.rm = TRUE), 2), "\n")
cat("PM10 (averaged):", round(mean(result6$pm10, na.rm = TRUE), 2), "\n")

cat("\n=== 7. Missing Data Handling ===\n")
result7_keep <- preprocess_baseline(demo_ukb, 
                                    variables = c("sex", "age"),
                                    missing_action = "keep")
result7_drop <- preprocess_baseline(demo_ukb, 
                                    variables = c("sex", "age"),
                                    missing_action = "drop")
cat("Rows with keep:", nrow(result7_keep), "\n")
cat("Rows with drop:", nrow(result7_drop), "\n")

cat("\n=== Complete Workflow Example ===\n")
# Process all baseline characteristics in one pipeline
final_data <- demo_ukb |>
  preprocess_baseline(variables = c("sex", "age", "ethnicity", "bmi", 
                                    "smoking", "drinking", "education")) |>
  calculate_blood_pressure(type = "sbp") |>
  calculate_blood_pressure(type = "dbp") |>
  extract_medications(c("cholesterol", "blood_pressure")) |>
  calculate_air_pollution(c("NO2", "PM2.5"))

cat("Final dataset columns:", ncol(final_data), "\n")
cat("New processed columns:\n")
new_cols <- setdiff(names(final_data), names(demo_ukb))
print(new_cols)

cat("\n=== Example completed successfully ===\n")