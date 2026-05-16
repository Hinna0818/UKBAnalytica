# `ukbsci-plot` — function reference

## `plot_forest()`

```r
plot_forest(results,
            estimate_col      = "estimate",
            lower_col         = "lower95",
            upper_col         = "upper95",
            label_col         = "subgroup",
            pvalue_col        = "pvalue",
            p_interaction_col = "p_interaction",
            null_value        = 1,
            log_scale         = TRUE,
            colors            = NULL,
            title             = "Subgroup Analysis",
            xlab              = "Hazard Ratio (95% CI)",
            show_n            = TRUE,
            show_events       = TRUE)
```

Consumes data.frames with `estimate, lower95, upper95, label_col`. Reverses
row order so the top of the input is the bottom of the figure (typical for
"primary subgroup at top").

**Returns:** ggplot.

---

## `plot_calibration()`

```r
plot_calibration(data, predicted, observed,
                 n_bins   = 10,
                 smooth   = TRUE,
                 conf_int = TRUE)
```

Per-prediction `data` is binned; the rendered output is aggregate. CIs via
`binom.test()` per bin; LOESS smooth requires ≥ 3 bins.

---

## `plot_regression_volcano()`

```r
plot_regression_volcano(data,
                        effect_col          = NULL,    # auto: HR → OR → estimate
                        p_col               = "pvalue",
                        adjusted_p_col      = NULL,
                        label_col           = NULL,    # auto: gene_symbol → protein_clean → variable
                        significance_cutoff = 0.05,
                        top_n_label_each    = 5,
                        null_effect         = 1,
                        x_lab               = NULL,
                        y_lab               = NULL,
                        x_limits            = NULL,
                        y_limits            = NULL,
                        point_size          = 1.05,
                        label_size          = 2,
                        colors              = c(neutral = "#D8D8D8",
                                                 lower   = "#2F6FA3",
                                                 higher  = "#C74732"),
                        show_cutoff         = TRUE)
```

Attaches `attr(., "plot_data")` and `attr(., "label_data")` to the returned
ggplot for downstream CSV export.

---

## Shared helpers (defined by this skill, not exported by the package)

```r
ukbsci_theme(base_size = 7, base_family = "Arial")
ukbsci_clinical    # named char vector of hex codes (categorical)
ukbsci_diverging   # length-7 vector (heatmaps / signed effects)
ukbsci_sequential  # length-8 vector (density / continuous)
save_ukbsci_figure(plot, name, width_mm = 180, height_mm = 120, dpi = 300,
                   outdir = "/mnt/project/<area>/05-figs")
```

Copy from `references/theme-and-palettes.R` at the top of every figure
script.

---

## Multi-panel composition

`patchwork` is the recommended engine:

```r
library(patchwork)
fig <- (p_forest | p_volcano) / (p_calib | p_dca) +
        plot_annotation(tag_levels = "a") &
        ukbsci_theme()
```

`& ukbsci_theme()` propagates the theme to every patch.
