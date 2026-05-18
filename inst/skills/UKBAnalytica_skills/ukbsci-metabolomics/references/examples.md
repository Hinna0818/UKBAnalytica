# Examples — ukbsci-metabolomics

## Example 1: End-to-end ORA with custom pathway library

Scenario: you ran `runmulti_cox()` on 249 Nightingale metabolites and obtained
a significant hit list at FDR < 5 %. You want to identify enriched metabolic
pathways.

```r
library(UKBAnalytica)

# --- Step 1: define the hit list ---
hits <- c(
  "Alanine", "Glutamine", "Glycine", "Isoleucine",
  "Lactate", "Pyruvate", "Citrate", "3-Hydroxybutyrate"
)

# --- Step 2: classify and check coverage ---
clf <- classify_metabolites(hits)
clf
# category "lipoprotein_lipid" and "protein" will be excluded from ORA

# --- Step 3: define pathway library ---
pathway_db <- data.frame(
  pathway = c(
    rep("Amino acid metabolism",  4),
    rep("Energy metabolism",      3),
    rep("Ketone body metabolism", 2),
    rep("TCA cycle",              2)
  ),
  metabolite = c(
    "L-Alanine", "L-Glutamine", "Glycine", "L-Isoleucine",
    "Lactic acid", "Pyruvic acid", "D-Glucose",
    "Acetoacetic acid", "3-Hydroxybutyric acid",
    "Citric acid", "Pyruvic acid"
  )
)

# --- Step 4: run ORA ---
ora <- run_metabolite_ora(
  metabolites     = hits,
  pathway_db      = pathway_db,
  backend         = "custom",
  min_metabolites = 2,
  p_adjust_method = "BH"
)

# Significant pathways
ora$ora_result[ora$ora_result$p_adjust < 0.05, ]

# --- Step 5: visualize ---
plot_metabolite_ora_dotplot(ora, top_n = 10)
plot_metabolite_ora_barplot(ora, top_n = 10)

# --- Step 6: save for manuscript ---
write.csv(ora$ora_result,
          "/mnt/project/<area>/04-results/metabolomics_ora.csv",
          row.names = FALSE)

ggsave("/mnt/project/<area>/04-results/metabolomics_dotplot.pdf",
       plot_metabolite_ora_dotplot(ora),
       width = 6, height = 4)
```

---

## Example 2: MetaboAnalystR backend with SMPDB pathways

Requires MetaboAnalystR installed and internet access.

```r
ora_ma <- run_metabolite_ora(
  metabolites    = hits,
  backend        = "metaboanalyst",
  library        = "smpdb_pathway",
  id_type        = "name",
  run_subprocess = TRUE
)

plot_metabolite_ora_dotplot(ora_ma, top_n = 15)
```

---

## Example 3: Custom name mapping for non-standard labels

```r
custom_map <- data.frame(
  metabolite         = c("3OHB", "AcAc"),
  metaboanalyst_name = c("3-Hydroxybutyric acid", "Acetoacetic acid")
)

hits_custom <- c("Alanine", "3OHB", "AcAc", "Lactate")

ora_custom <- run_metabolite_ora(
  metabolites   = hits_custom,
  pathway_db    = pathway_db,
  mapping_table = custom_map,
  backend       = "custom",
  min_metabolites = 2
)

ora_custom$mapping   # verify all names were mapped
ora_custom$unmatched # should be empty
```
