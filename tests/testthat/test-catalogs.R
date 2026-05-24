test_that("disease catalog exposes curated and pomegranate sources", {
  all_codes <- get_disease_catalog()
  curated <- get_disease_catalog(source = "curated")
  pomegranate <- get_disease_catalog(source = "pomegranate")

  expect_s3_class(all_codes, "data.frame")
  expect_gt(nrow(all_codes), 22000)
  expect_gt(nrow(curated), 100)
  expect_gt(nrow(pomegranate), 22000)
  expect_true(all(c("curated", "pomegranate") %in% unique(all_codes$source)))
  expect_true(all(c("definition_id", "disease_name", "code_system", "code") %in% names(all_codes)))
})

test_that("Pomegranate provenance and portal audit table are available", {
  manifest <- get_pomegranate_source_manifest()
  portal <- load_pomegranate_portal_coding()

  expect_s3_class(manifest, "data.frame")
  expect_true(all(c("pomegranate_github_yaml", "pomegranate_portal_csv") %in% manifest$source))
  expect_true(any(manifest$role == "canonical_disease_catalog"))
  expect_s3_class(portal, "data.frame")
  expect_gt(nrow(portal), 26000)
  expect_true(all(c("disease", "field_or_source", "codes", "phenotype_url") %in% names(portal)))
})

test_that("disease catalog can be searched and converted to pomegranate definitions", {
  copd <- get_disease_catalog(disease = "copd")
  expect_gt(nrow(copd), 0)
  expect_true(any(copd$source == "pomegranate"))

  asthma_defs <- get_pomegranate_diseases("asthma")
  expect_type(asthma_defs, "list")
  expect_true("asthma" %in% names(asthma_defs))
  expect_equal(asthma_defs$asthma$name, "Asthma")
  expect_match(asthma_defs$asthma$icd10_pattern, "J45")
})

test_that("medication catalog includes curated and official coding entries", {
  meds <- get_medication_catalog()
  metformin <- get_medication_catalog("metformin")

  expect_s3_class(meds, "data.frame")
  expect_gt(nrow(meds), 7000)
  expect_gt(nrow(metformin), 0)
  expect_true(all(c("medication_id", "medication_name", "code", "field_id") %in% names(meds)))
  expect_true(any(meds$is_default))
  expect_true(any(!meds$is_default))
})

test_that("existing predefined disease behavior remains unchanged", {
  defs <- get_predefined_diseases()
  expect_type(defs, "list")
  expect_equal(length(defs), 78)
  expect_true("COPD" %in% names(defs))
  expect_equal(defs$COPD$name, "Chronic Obstructive Pulmonary Disease")
})

test_that("get_predefined_diseases supports curated, pomegranate, and both sources", {
  curated <- get_predefined_diseases(source = "curated")
  pomegranate <- get_predefined_diseases(source = "pomegranate")
  both <- get_predefined_diseases(source = "both", merge_type = "intersection")
  both_union <- get_predefined_diseases(source = "both", merge_type = "union")

  expect_equal(length(curated), 78)
  expect_gt(length(pomegranate), 250)
  expect_gt(length(both), 40)
  expect_gte(length(both_union), length(pomegranate))
  expect_true(all(names(both) %in% names(both_union)))
  expect_true("aki" %in% names(both_union))
  expect_true("TAA" %in% names(both_union))
  expect_true(all(c("Asthma", "COPD") %in% names(both)))
  expect_equal(both$Asthma$name, "Asthma")
  expect_match(both$Asthma$icd10_pattern, "J45")
  expect_match(both$COPD$icd10_pattern, "J44")
  expect_match(both_union$MI$icd10_pattern, "I252")
  expect_match(both_union$CKD$icd10_pattern, "I120")
  expect_match(both_union$PAD$opcs4_pattern, "L50")
  expect_true(all(vapply(both, function(x) is.list(x) && !is.null(x$name), logical(1))))

  copd_union <- get_predefined_diseases(
    source = "both",
    merge_type = "union",
    disease = "COPD"
  )
  expect_equal(names(copd_union), "COPD")
  expect_match(copd_union$COPD$icd10_pattern, "J44")
  expect_true(1113 %in% copd_union$COPD$sr_codes)
})
