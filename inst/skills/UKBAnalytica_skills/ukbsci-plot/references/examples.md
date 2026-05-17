# `ukbsci-plot` — examples

```r
library(UKBAnalytica); library(ggplot2); library(patchwork); library(data.table)
# Load shared theme helpers (adjust path to your repo root)
source(file.path(getwd(), "inst/skills/UKBAnalytica_skills/ukbsci-plot/references/theme-and-palettes.R"))
```

## A. Forest plot from subgroup output

```r
sub <- fread("/mnt/project/<area>/04-results/03-cox_subgroup.csv")
p_forest <- plot_forest(
  results = sub,
  estimate_col = "estimate", lower_col = "lower95", upper_col = "upper95",
  label_col = "subgroup", pvalue_col = "pvalue",
  p_interaction_col = "p_interaction",
  null_value = 1, log_scale = TRUE,
  title = "Subgroup HRs for incident AA per 10-mmHg SBP",
  xlab  = "Hazard Ratio (95% CI)"
) + ukbsci_theme()
save_ukbsci_figure(p_forest, "Fig02-forest_subgroup",
                   width_mm = 180, height_mm = 140)
```

## B. Calibration plot from ML test predictions

```r
preds <- flow$final_test_predictions       # RAP-resident
p_cal <- plot_calibration(
  data      = preds,
  predicted = "prob",
  observed  = "truth",
  n_bins    = 10,
  smooth    = TRUE,
  conf_int  = TRUE
) + ukbsci_theme()
save_ukbsci_figure(p_cal, "Fig14-calibration",
                   width_mm = 89, height_mm = 89)

# Aggregate source CSV (the binned data)
src <- ggplot2::layer_data(p_cal)
fwrite(src, "/mnt/project/<area>/05-figs/data/Fig14-calibration.csv")
```

## C. Volcano of a protein screen

```r
prot <- fread("/mnt/project/<area>/04-results/protein_screen.csv")
p_volc <- plot_regression_volcano(
  prot,
  effect_col   = "HR",
  p_col        = "pvalue",
  adjusted_p_col = "padj",
  label_col    = "gene_symbol",
  significance_cutoff = 0.05,
  top_n_label_each = 10,
  null_effect  = 1,
  colors = c(neutral = ukbsci_clinical["control"],
              lower   = ukbsci_clinical["baseline"],
              higher  = ukbsci_clinical["exposure"]),
  show_cutoff = TRUE
)
save_ukbsci_figure(p_volc, "Fig05-volcano_proteins",
                   width_mm = 120, height_mm = 100)

# Aggregate source — use the plot_data attribute
fwrite(attr(p_volc, "plot_data"),
       "/mnt/project/<area>/05-figs/data/Fig05-volcano_proteins.csv")
```

## D. Composite Figure 1 (Forest + Volcano + Calibration + DCA)

```r
fig <- (p_forest | p_volc) / (p_cal | p_dca) +
        plot_annotation(tag_levels = "a") &
        ukbsci_theme()
save_ukbsci_figure(fig, "Fig01-composite",
                   width_mm = 180, height_mm = 200)
```
