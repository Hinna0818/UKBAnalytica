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
    "Colorectal_Cancer", "Colon_Cancer", "Rectal_Cancer", "Stomach_Cancer",
    "Uterus_Cancer", "MACE", "Hyperglycaemia", "PCOS",
    "Systemic_Lupus_Erythematosus"
  )
  expect_true(all(expected %in% names(defs)))
  expect_match(defs$Atrial_Fibrillation$icd10_pattern, "I48", fixed = TRUE)
  expect_false(grepl("J48", defs$Atrial_Fibrillation$icd10_pattern, fixed = TRUE))
  expect_equal(defs$Breast_Cancer$cancer_behaviour, 3L)
  expect_true(all(c(1276, 1468, 1607) %in% defs$Diabetes$sr_codes))
  expect_match(defs$Hyperglycaemia$icd10_pattern, "R73", fixed = TRUE)
  expect_equal(defs$PCOS$sr_codes, 1350)
  expect_match(defs$Systemic_Lupus_Erythematosus$icd10_pattern, "M32", fixed = TRUE)
})
