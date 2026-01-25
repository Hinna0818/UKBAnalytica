# Tests for dual-source functionality using real demo data
# Uses raw_demo.csv from testdata folder

test_that("build_survival_dataset with real demo data works correctly", {
  skip_if_not_installed("data.table")
  
  # Load real test data
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  n <- nrow(demo_data)
  
  # Get disease definitions
  diseases <- get_predefined_diseases()[c("Hypertension", "Diabetes")]
  
  # Test 1: Default dual-source behavior
  result1 <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # Basic structure tests
  expect_s3_class(result1, "data.table")
  expect_equal(nrow(result1), n)
  expect_true(all(c("Hypertension_history", "Hypertension_incident", 
                    "outcome_status", "outcome_surv_time") %in% names(result1)))
  expect_true(all(result1$Hypertension_history %in% c(0, 1)))
  expect_true(all(result1$Hypertension_incident %in% c(0, 1)))
  expect_true(all(result1$outcome_status %in% c(0, 1)))
  expect_true(all(result1$outcome_surv_time >= 0, na.rm = TRUE))
  
  # Test 2: Hospital records only (strictest)
  result2 <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    prevalent_sources = c("ICD10"),
    outcome_sources = c("ICD10"),
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  expect_s3_class(result2, "data.table")
  expect_equal(nrow(result2), n)
  
  # Test 3: All sources including self-report for outcome
  result3 <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    prevalent_sources = c("ICD10", "Self-report"),
    outcome_sources = c("ICD10", "Self-report"),
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  expect_s3_class(result3, "data.table")
  expect_equal(nrow(result3), n)
  
  # Test 4: Logical consistency
  # Default (with self-report for prevalent) should have same or more prevalent cases
  expect_gte(sum(result1$Hypertension_history), sum(result2$Hypertension_history))
  
  # All participants should be retained
  expect_equal(nrow(result1), n)
  expect_equal(nrow(result2), n)
  expect_equal(nrow(result3), n)
  
  # Outcome status should be binary
  expect_true(all(result1$outcome_status %in% c(0, 1)))
  expect_true(all(result2$outcome_status %in% c(0, 1)))
  expect_true(all(result3$outcome_status %in% c(0, 1)))
  
  # Test 5: Long format output
  long_result <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    primary_disease = "Hypertension",
    output = "long",
    include_all = TRUE,
    censor_date = as.Date("2023-12-31")
  ))
  
  expect_s3_class(long_result, "data.table")
  expect_true("disease" %in% names(long_result))
  expect_true(all(c("Hypertension", "Diabetes") %in% unique(long_result$disease)))
  expect_true(all(long_result$status %in% c(0, 1)))
  
  # Test 6: Parameter validation
  expect_error(
    suppressMessages(build_survival_dataset(
      dt = demo_data,
      disease_definitions = diseases,
      primary_disease = "NonexistentDisease"
    )),
    "Primary disease"
  )
  
  # Test 7: Prevalent case exclusion logic
  analysis_dt <- result1[Hypertension_history == 0]  # Exclude prevalent
  expect_lte(nrow(analysis_dt), n)  # Should be <= original size
  
  # Events in analysis cohort should be valid
  if (nrow(analysis_dt) > 0) {
    expect_true(all(analysis_dt$outcome_status %in% c(0, 1)))
    expect_true(all(analysis_dt$outcome_surv_time > 0, na.rm = TRUE))
  }
})

test_that("prevalent_sources and outcome_sources work independently", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  n <- nrow(demo_data)
  diseases <- get_predefined_diseases()["Hypertension"]
  
  # Method 1: Self-report for prevalent, ICD10 for outcome (default)
  result_default <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
    outcome_sources = c("ICD10", "ICD9", "Death"),
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # Method 2: ICD10 only for both
  result_strict <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    prevalent_sources = c("ICD10"),
    outcome_sources = c("ICD10"),
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # Both should have same number of rows (full cohort)
  expect_equal(nrow(result_default), nrow(result_strict))
  expect_equal(nrow(result_default), n)
  
  # Default (with self-report) should have >= prevalent cases
  expect_gte(
    sum(result_default$Hypertension_history),
    sum(result_strict$Hypertension_history)
  )
  
  # Both should produce valid binary outcomes
  expect_true(all(result_default$outcome_status %in% c(0, 1)))
  expect_true(all(result_strict$outcome_status %in% c(0, 1)))
})


test_that("prevalent cases are excluded from outcome_status", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  diseases <- get_predefined_diseases()["Hypertension"]
  
  result <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # Participants with history=1 should have outcome_status=0
  # (prevalent cases should not be counted as incident outcomes)
  prevalent_outcomes <- result[Hypertension_history == 1 & outcome_status == 1]
  expect_equal(nrow(prevalent_outcomes), 0)
})


test_that("self-report contributes to history but not outcome", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  diseases <- get_predefined_diseases()["Hypertension"]
  
  # With self-report for prevalent only
  result_with_sr <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    prevalent_sources = c("ICD10", "Self-report"),
    outcome_sources = c("ICD10"),
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # Without self-report
  result_no_sr <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    prevalent_sources = c("ICD10"),
    outcome_sources = c("ICD10"),
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # With self-report should have more or equal history cases
  expect_gte(
    sum(result_with_sr$Hypertension_history),
    sum(result_no_sr$Hypertension_history)
  )
})


test_that("multiple diseases can be processed together", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  
  # Test with multiple diseases
  diseases <- get_predefined_diseases()[c("Hypertension", "Diabetes", "Stroke")]
  
  result <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # Should have columns for all diseases
  for (d in names(diseases)) {
    expect_true(paste0(d, "_history") %in% names(result))
    expect_true(paste0(d, "_incident") %in% names(result))
  }
  
  # Primary disease determines outcome columns
  expect_true("outcome_status" %in% names(result))
  expect_true("outcome_surv_time" %in% names(result))
})


test_that("censor_date parameter works correctly", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  diseases <- get_predefined_diseases()["Hypertension"]
  
  # Early censor date
  result_early <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    primary_disease = "Hypertension",
    censor_date = as.Date("2015-12-31")
  ))
  
  # Late censor date
  result_late <- suppressMessages(build_survival_dataset(
    dt = demo_data,
    disease_definitions = diseases,
    primary_disease = "Hypertension",
    censor_date = as.Date("2023-12-31")
  ))
  
  # Later censor date should result in longer mean follow-up
  mean_early <- mean(result_early$outcome_surv_time, na.rm = TRUE)
  mean_late <- mean(result_late$outcome_surv_time, na.rm = TRUE)
  expect_lte(mean_early, mean_late)
})