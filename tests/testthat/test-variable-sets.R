if (!exists("get_variable_sets", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "variable_sets.R"), local = FALSE)
}

test_that("curated variable sets expose common UKB categories", {
  vars <- get_variable_sets()

  expect_true(all(c(
    "set", "category", "variable", "field_id", "ukb_col", "label"
  ) %in% names(vars)))
  expect_true(nrow(vars) > 100)
  expect_true(all(c(
    "clinical_core", "air_pollution", "family_history",
    "algorithm_outcomes", "self_report_medication", "biomarkers_core"
  ) %in% unique(vars$set)))
  expect_true(all(c(93L, 94L, 4079L, 4080L) %in% get_variable_set("clinical_core", output = "field_id")))
  expect_true(all(c(30750L, 30760L, 30780L, 30870L) %in% get_variable_set("biomarkers_core", output = "field_id")))

  air_ids <- get_variable_set("air_pollution", output = "field_id")
  expect_true(all(c(24003L, 24004L, 24005L, 24006L, 24019L) %in% air_ids))

  family <- get_variable_set("family_history")
  expect_true(all(c(20107L, 20110L, 20111L) %in% family$field_id))
})

test_that("curated variable set filters validate inputs", {
  expect_error(get_variable_sets("not_a_set"), "Unknown variable set")
  expect_error(get_variable_sets(category = "not_a_category"), "Unknown category")
  expect_equal(get_variable_set("genetic_pcs_20", output = "ukb_col")[[1]], "p22009_a1")
})
