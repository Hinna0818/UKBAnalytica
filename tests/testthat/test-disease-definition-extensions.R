if (!exists("get_predefined_diseases", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "first_occurrence.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "cancer_registry.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "disease_definitions.R"), local = FALSE)
}

test_that("predefined diseases include expanded chronic and cancer definitions", {
  defs <- get_predefined_diseases()

  expected <- c(
    "Dyspepsia", "Inflammatory_Bowel_Disease", "Diverticular_Disease",
    "Thyroid_Disorders", "Migraine", "Bronchiectasis", "Multiple_Sclerosis",
    "Glaucoma", "Cataract", "AMD", "Breast_Cancer", "Prostate_Cancer",
    "Colorectal_Cancer", "Stomach_Cancer"
  )
  expect_true(all(expected %in% names(defs)))
  expect_match(defs$Atrial_Fibrillation$icd10_pattern, "I48", fixed = TRUE)
  expect_false(grepl("J48", defs$Atrial_Fibrillation$icd10_pattern, fixed = TRUE))
  expect_equal(defs$Breast_Cancer$cancer_behaviour, 3L)
})
