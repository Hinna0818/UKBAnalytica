# Tests for variable preprocessing functions using raw_demo.csv

test_that("preprocess_baseline works with real demo data", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  n <- nrow(demo_data)
  
  # Test 1: Basic demographics processing
  result <- preprocess_baseline(demo_data, 
                                variables = c("sex", "age", "ethnicity"),
                                missing_action = "keep")
  
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), n)
  expect_true(all(c("sex", "age", "ethnicity") %in% names(result)))
  
  # Sex should be factor with Female/Male levels
  expect_s3_class(result$sex, "factor")
  expect_true(all(levels(result$sex) %in% c("Female", "Male")))
  
  # Age should be numeric
  expect_type(result$age, "double")
  
  # Ethnicity should be factor
  expect_s3_class(result$ethnicity, "factor")
})


test_that("preprocess_baseline handles missing_action correctly", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  
  # Test keep option
  result_keep <- preprocess_baseline(demo_data, 
                                     variables = c("sex", "age"),
                                     missing_action = "keep")
  expect_equal(nrow(result_keep), nrow(demo_data))
  
  # Test drop option (should have fewer or equal rows)
  result_drop <- suppressMessages(preprocess_baseline(demo_data, 
                                                      variables = c("sex", "age"),
                                                      missing_action = "drop"))
  expect_lte(nrow(result_drop), nrow(demo_data))
})


test_that("preprocess_baseline supports custom mapping", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  
  # Define custom mapping
  custom <- list(
    my_age_copy = list(ukb_col = "p21022", description = "Custom age variable")
  )
  
  result <- preprocess_baseline(demo_data, 
                                variables = c("sex", "my_age_copy"),
                                custom_mapping = custom)
  
  expect_true("my_age_copy" %in% names(result))
  expect_equal(result$my_age_copy, demo_data$p21022)
})


test_that("preprocess_baseline validates input", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  
  # Test invalid variable name

  expect_error(
    preprocess_baseline(demo_data, variables = c("invalid_variable")),
    "Unknown variables"
  )
  
  # Test invalid df input
  expect_error(
    preprocess_baseline("not_a_df", variables = c("sex")),
    "must be a data.frame"
  )
  
  # Test invalid custom_mapping
  expect_error(
    preprocess_baseline(demo_data, 
                        variables = c("bad_custom"),
                        custom_mapping = list(bad_custom = list(no_ukb_col = "test"))),
    "must include 'ukb_col'"
  )
})


test_that("get_variable_info returns correct structure", {
  # Test all categories
  all_info <- get_variable_info("all")
  expect_s3_class(all_info, "data.frame")
  expect_true(all(c("variable", "category", "field_id", "ukb_column", "description") %in% names(all_info)))
  
  # Test specific category
  demo_info <- get_variable_info("demographics")
  expect_equal(nrow(demo_info), 3)  # sex, age, ethnicity
  expect_true(all(demo_info$variable %in% c("sex", "age", "ethnicity")))
  
  # Test invalid category
  expect_error(get_variable_info("invalid_category"), "Invalid category")
})


test_that("calculate_blood_pressure works correctly", {
  skip_if_not_installed("data.table")
  
  # Create test data with BP columns
  test_data <- data.table::data.table(
    eid = 1:10,
    p4080_i0_a0 = c(120, 130, NA, 140, 125, 135, NA, 128, 132, 138),
    p4080_i0_a1 = c(118, 128, NA, 138, 123, 133, NA, 126, 130, 136),
    p93_i0_a0 = c(122, 132, 125, 142, 127, 137, 130, 130, 134, 140),
    p93_i0_a1 = c(120, 130, 123, 140, 125, 135, 128, 128, 132, 138)
  )
  
  result <- calculate_blood_pressure(test_data, type = "sbp")
  
  expect_true("sbp" %in% names(result))
  expect_type(result$sbp, "double")
  
  # Check that NA manual values are filled with auto
  expect_false(is.na(result$sbp[3]))  # Should be filled from auto
  expect_false(is.na(result$sbp[7]))  # Should be filled from auto
})


test_that("extract_medications works correctly", {
  skip_if_not_installed("data.table")
  
  # Create test data
  test_data <- data.table::data.table(
    eid = 1:5,
    p6177_i0 = c("[1]", "[2]", "[1, 2]", "[-1]", NA),
    p6153_i0 = c("[2]", "[3]", NA, "[1]", "[2, 3]")
  )
  
  result <- extract_medications(test_data, c("cholesterol", "blood_pressure", "insulin"))
  
  expect_true(all(c("med_cholesterol", "med_blood_pressure", "med_insulin") %in% names(result)))
  
  # Check first row: drug1=[1], drug2=[2] -> cholesterol=1, bp=1, insulin=0
  expect_equal(result$med_cholesterol[1], 1L)
  expect_equal(result$med_blood_pressure[1], 1L)
  expect_equal(result$med_insulin[1], 0L)
  
  # Check row 4: drug1=[-1] -> NA
  expect_true(is.na(result$med_cholesterol[4]))
})


test_that("calculate_air_pollution works correctly", {
  skip_if_not_installed("data.table")
  
  # Create test data with pollution columns
  test_data <- data.table::data.table(
    eid = 1:5,
    p24016 = c(25, 30, 28, 22, 27),  # NO2 2005
    p24018 = c(24, 29, 27, 21, 26),  # NO2 2006
    p24017 = c(23, 28, 26, 20, 25),  # NO2 2007
    p24003 = c(22, 27, 25, 19, 24),  # NO2 2010
    p24006 = c(12, 15, 13, 11, 14)   # PM2.5
  )
  
  result <- calculate_air_pollution(test_data, c("NO2", "PM2.5"))
  
  expect_true("no2" %in% names(result))
  expect_true("pm25" %in% names(result))
  
  # NO2 should be average of 4 columns
  expected_no2 <- rowMeans(test_data[, .(p24016, p24018, p24017, p24003)])
  expect_equal(result$no2, expected_no2)
  
  # PM2.5 should be same as p24006
  expect_equal(result$pm25, test_data$p24006)
})


test_that("preprocess_baseline handles lifestyle variables", {
  skip_if_not_installed("data.table")
  
  test_data_path <- test_path("testdata", "raw_demo.csv")
  skip_if_not(file.exists(test_data_path), "Test data file raw_demo.csv not found")
  
  demo_data <- data.table::fread(test_data_path)
  
  # Check if smoking column exists before testing
  if ("p20116_i0" %in% names(demo_data)) {
    result <- preprocess_baseline(demo_data, variables = c("smoking"))
    
    expect_true("smoking" %in% names(result))
    expect_s3_class(result$smoking, "factor")
    expect_true(all(levels(result$smoking) %in% c("Never", "Previous", "Current")))
  }
})
