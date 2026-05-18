test_that("metabolite panel mapping and custom ORA workflow run", {
  panel <- load_ukb_metabolite_panel()
  expect_true(nrow(panel) > 0)
  expect_true(all(c("Description", "UKB_ID", "meta_ID") %in% names(panel)))

  cls <- classify_metabolites(panel$Description)
  expect_true(any(cls$category == "small_molecule"))
  expect_true(any(cls$category == "lipoprotein_lipid"))

  hits <- c(
    "Alanine", "Glutamine", "Glycine", "Histidine",
    "Isoleucine", "Leucine", "Valine", "Phenylalanine",
    "Tyrosine", "Glucose", "Lactate", "Pyruvate", "Citrate"
  )
  pathway_db <- data.frame(
    pathway = c(rep("Amino acid metabolism", 9), rep("Energy metabolism", 4)),
    metabolite = c(
      "L-Alanine", "L-Glutamine", "Glycine", "L-Histidine",
      "L-Isoleucine", "L-Leucine", "L-Valine", "L-Phenylalanine",
      "L-Tyrosine", "D-Glucose", "Lactic acid", "Pyruvic acid",
      "Citric acid"
    ),
    stringsAsFactors = FALSE
  )

  res <- run_metabolite_ora(
    metabolites = hits,
    pathway_db = pathway_db,
    universe = panel$Description,
    backend = "custom",
    min_metabolites = 3
  )

  expect_s3_class(res, "ukb_metabolite_ora")
  expect_true(nrow(res$ora_result) > 0)
  expect_true(all(c("pathway", "hits", "pvalue", "p_adjust", "fold_enrichment") %in% names(res$ora_result)))
  expect_s3_class(plot_metabolite_ora_dotplot(res), "ggplot")
  expect_s3_class(plot_metabolite_ora_barplot(res), "ggplot")
})
