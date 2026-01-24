test_that("create_disease_definition works correctly", {
  def <- create_disease_definition(
    name = "Test Disease",
    icd10_pattern = "^I71",
    icd9_pattern = "^441",
    sr_codes = c(1065)
  )
  
  expect_type(def, "list")
  expect_equal(def$name, "Test Disease")
  expect_equal(def$icd10_pattern, "^I71")
  expect_equal(def$icd9_pattern, "^441")
  expect_equal(def$sr_codes, c(1065))
})

test_that("get_predefined_diseases returns valid list", {
  diseases <- get_predefined_diseases()
  
  expect_type(diseases, "list")
  expect_true("AA" %in% names(diseases))
  expect_true("Hypertension" %in% names(diseases))
  expect_true("Diabetes" %in% names(diseases))
  
  # Check structure of a disease definition
  aa <- diseases$AA
  expect_equal(aa$name, "Aortic Aneurysm")
  expect_equal(aa$icd10_pattern, "^I71")
})

test_that("combine_disease_definitions merges correctly", {
  diseases <- get_predefined_diseases()
  mace <- combine_disease_definitions(
    diseases$MI, diseases$Stroke,
    name = "MACE"
  )
  
  expect_equal(mace$name, "MACE")
  expect_true(grepl("I21", mace$icd10_pattern))
  expect_true(grepl("I60", mace$icd10_pattern))
})
