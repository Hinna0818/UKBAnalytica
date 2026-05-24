if (!exists("get_predefined_medications", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "medication_definitions.R"), local = FALSE)
}

test_that("predefined medication definitions expose field 20003 code groups", {
  meds <- get_predefined_medications()

  expect_true(all(c(
    "Blood_Pressure_Medication", "Diabetes_Medication", "ACE_Inhibitor",
    "ARB", "Beta_Blocker", "Calcium_Channel_Blocker", "Thiazide_Diuretic",
    "Insulin", "Metformin", "Sulfonylurea"
  ) %in% names(meds)))
  expect_true("1140883066" %in% meds$Insulin$codes)
  expect_true("1140874686" %in% meds$Metformin$codes)
  expect_true("1140851690" %in% meds$ACE_Inhibitor$codes)
  expect_equal(meds$Insulin$field_id, 20003L)
})

test_that("predefined medication codes are present in official coding 4", {
  coding4 <- load_ukb_medication_coding()
  meds <- get_predefined_medications()
  all_codes <- unique(unlist(lapply(meds, `[[`, "codes"), use.names = FALSE))

  expect_true(all(c("coding", "meaning") %in% names(coding4)))
  expect_true(nrow(coding4) > 6000)
  expect_setequal(setdiff(all_codes, coding4$coding), character(0))
})

test_that("extract_self_report_medications creates wide and long indicators", {
  dat <- data.frame(
    eid = 1:4,
    p20003_i0_a0 = c("1140883066", "1140874686", "1140851690", NA),
    p20003_i0_a1 = c(NA, "1140857494", NA, NA),
    stringsAsFactors = FALSE
  )

  wide <- extract_self_report_medications(
    dat,
    medications = c("Insulin", "Metformin", "ACE_Inhibitor", "Diabetes_Medication")
  )

  expect_equal(wide$med20003_insulin, c(1L, 0L, 0L, 0L))
  expect_equal(wide$med20003_metformin, c(0L, 1L, 0L, 0L))
  expect_equal(wide$med20003_ace_inhibitor, c(0L, 0L, 1L, 0L))
  expect_equal(wide$med20003_diabetes_medication, c(1L, 1L, 0L, 0L))

  long <- extract_self_report_medications(
    dat,
    medications = c("Insulin", "Metformin"),
    return_long = TRUE
  )
  expect_true(all(c("eid", "medication", "has_medication") %in% names(long)))
  expect_equal(nrow(long), 8L)
})
