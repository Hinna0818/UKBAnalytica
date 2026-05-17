---
name: ukbsci-plot
description: >
  Publication-grade figure production for UK Biobank Research Analysis
  Platform (RAP) analyses built with UKBAnalytica. Wraps the package's
  built-in plotters — plot_forest (subgroup / regression forest),
  plot_calibration (clinical-model calibration), plot_regression_volcano
  (multi-exposure / multi-protein volcano), plus the family of survival,
  ML, SHAP, mediation, propensity, imputation, correlation, and
  enrichment plots from sibling skills — and adds a shared neutral theme
  (ukbsci_clinical / ukbsci_diverging / ukbsci_sequential palettes), a
  dual-format save helper, and a figure-contract checklist. Use this skill
  when the user asks for any manuscript figure derived from UKBAnalytica
  outputs: forest plots, volcano plots, KM curves, calibration plots,
  Love plots, SHAP summaries, multi-panel composites, or PDF / SVG
  high-resolution export. Triggers: UKB plotting, forest plot, volcano,
  calibration, manuscript figure, ggplot helper, ukbsci 画图, 论文级图,
  /ukbsci-plot. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-plot — Manuscript figure production for UKBAnalytica

## 0. RAP guardrails

**Core rule:** do not export UK Biobank participant-level raw data, direct
identifiers (`eid`, exact dates, raw RAP fields), or row-level source tables
that can be linked back to participants.

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

| Output | Sharing rule |
|--------|-------------|
| Rendered PDF / SVG / PNG figure | Shareable if it does not expose identifiers, exact dates, raw RAP fields, row labels, or data previews. |
| `05-figs/data/Fig*.csv` source CSV | Shareable **only if aggregate** (per-bin, per-strata, per-feature means/CIs). Never share a source CSV that contains raw per-participant values or identifiers. |
| Coefficient / metric table from `04-results/*.csv` | Shareable — aggregate only. |
| Per-participant raw values used to build a figure | Do not export. Keep intermediate data on RAP. |
| Individual-level interpretation figure (e.g. SHAP force plot) | Do not use a real participant row with the local agent. Generate synthetic or representative prototype-row figures inside RAP if needed. |

---

## 1. When to load

- Build the manuscript-ready forest plot of the main Cox / subgroup result.
- Produce calibration / DCA / ROC for a clinical model (paired with
  `ukbsci-ml`).
- Volcano of a protein / metabolite / batch regression screen.
- Polished multi-panel figures combining outputs from multiple skills.
- Standardized save (`.pdf + .png + .svg`) at 300 dpi.

## 2. When NOT to load

- The user wants a KM curve only → `ukbsci-survival::plot_km_curve` (this
  skill defers to it).
- The user wants Love / PS distribution → `ukbsci-propensity` plots.
- The user wants mediation / MI / SHAP plots → use the originating skill's
  plotter; this skill provides the **shared theme and save helper**.

---

## 3. Figure contract (always run first)

For any figure the user asks for, restate three things before plotting:

1. **Claim** — the one-sentence conclusion the figure must support.
2. **Evidence panel(s)** — which underlying CSV(s) from `04-results/` feed
   the figure.
3. **Export contract** — final size (mm), single-column / double-column,
   font family, palette family.

A figure that cannot be re-rendered from `05-figs/data/Fig*.csv` is **not
acceptable**.

---

## 4. Shared theme & palettes (define once, reuse everywhere)

`UKBAnalytica` does not ship custom themes. This skill defines the
following neutral helpers — copy verbatim into each plotting script:

```r
library(ggplot2)

ukbsci_theme <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line       = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.35, colour = "black"),
      legend.title    = element_text(size = base_size * 0.95),
      legend.text     = element_text(size = base_size * 0.85),
      strip.text      = element_text(size = base_size * 0.95, face = "bold"),
      plot.title      = element_text(size = base_size * 1.1, face = "bold"),
      panel.grid      = element_blank()
    )
}

# Categorical palette — clinical / two-arm + accent
ukbsci_clinical <- c(
  baseline   = "#2F6FA3",
  exposure   = "#C74732",
  control    = "#7F7F7F",
  accent     = "#F2A93B",
  highlight1 = "#3D8B66",
  highlight2 = "#7A4E9B"
)

# Diverging palette — for heatmaps / SMD / signed effects
ukbsci_diverging <- c("#2166AC", "#67A9CF", "#D1E5F0",
                      "#F7F7F7",
                      "#FDDBC7", "#EF8A62", "#B2182B")

# Sequential palette — for density / continuous gradients
ukbsci_sequential <- c("#FFF7EC", "#FEE8C8", "#FDD49E",
                       "#FDBB84", "#FC8D59", "#EF6548",
                       "#D7301F", "#990000")

# Apply the theme once, globally, at the top of a script
theme_set(ukbsci_theme())
```

These names are intentionally generic — they **do not** reference any
specific journal brand.

---

## 5. Standardized save helper

```r
save_ukbsci_figure <- function(plot, name,
                                width_mm = 180, height_mm = 120,
                                dpi = 300,
                                outdir = "/mnt/project/<area>/05-figs") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  w <- width_mm / 25.4; h <- height_mm / 25.4

  # SVG (editable)
  svglite::svglite(file.path(outdir, paste0(name, ".svg")), width = w, height = h)
  print(plot); dev.off()

  # PDF (vector, cairo for fonts)
  grDevices::cairo_pdf(file.path(outdir, paste0(name, ".pdf")),
                       width = w, height = h, family = "Arial")
  print(plot); dev.off()

  # PNG (raster, 300 dpi)
  ragg::agg_png(file.path(outdir, paste0(name, ".png")),
                width = w, height = h, units = "in", res = dpi)
  print(plot); dev.off()
}
```

Manuscript single-column ≈ `width_mm = 89`; double-column ≈ `180`. Use a
height that yields golden-ratio-ish aspect unless the figure dictates
otherwise.

---

## 6. UKBAnalytica plot families (where each one lives)

| Plot | Provided by | Notes |
|------|-------------|-------|
| Forest (subgroup / batch reg) | `UKBAnalytica::plot_forest` | This skill — see §7 |
| Calibration | `UKBAnalytica::plot_calibration` | This skill — see §8 |
| Volcano | `UKBAnalytica::plot_regression_volcano` | This skill — see §9 |
| KM survival | `plot_km_curve` | `ukbsci-survival` |
| Love / PS distribution | `plot_balance`, `plot_ps_distribution` | `ukbsci-propensity` |
| Mediation effects | `plot_mediation`, `plot_mediation_forest` | `ukbsci-mediation` |
| MI diagnostics | `plot_mi_pooled`, `plot_mi_diagnostics` | `ukbsci-imputation` |
| GO bar / lollipop | `plot_go_ora_bar`, `plot_enrichment_lollipop` | `ukbsci-proteomics` |
| Correlation heatmap | `plot_correlation` | this skill (re-export-friendly) |
| ML evaluation | `plot_ml_roc`, `plot_ml_calibration`, `plot_ml_dca`, … | `ukbsci-ml` |
| SHAP | `plot_shap_summary`, `plot_shap_beeswarm`, `plot_shap_dependence`, `plot_shap_force` | `ukbsci-ml` |

When mixing plots from sibling skills, **apply the same theme** via
`theme_set(ukbsci_theme())` at the top of the orchestrating script and
override colors using the `ukbsci_*` palettes.

---

## 7. Forest plot — `plot_forest()`

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

Input: any data.frame with at least `estimate, lower95, upper95, label`.
Naturally consumes outputs of `runmulti_cox / run_subgroup_analysis /
run_multi_subgroup / runmulti_competing / pool_mi_models$pooled / tidy()`.

---

## 8. Calibration — `plot_calibration()`

```r
plot_calibration(data, predicted, observed,
                 n_bins   = 10,
                 smooth   = TRUE,
                 conf_int = TRUE)
```

Input: per-prediction rows with predicted probabilities and binary
outcome. **Note:** this is participant-level so the `data` argument should
come from `flow$final_test_predictions` (which stays on RAP); the figure
itself summarises into `n_bins` bins, so the rendered output is aggregate.

---

## 9. Volcano — `plot_regression_volcano()`

```r
plot_regression_volcano(data,
                        effect_col          = NULL,
                        p_col               = "pvalue",
                        adjusted_p_col      = NULL,
                        label_col           = NULL,
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

Auto-detects `effect_col` (HR → OR → estimate) and `label_col`
(`gene_symbol → protein_clean → variable`). The output's
`attr(., "plot_data")` and `attr(., "label_data")` give the aggregate
underlying tables — write those to `05-figs/data/`.

---

## 10. Multi-panel composition

`patchwork` is the recommended layout engine:

```r
library(patchwork)
top  <- p_forest + p_volcano + plot_layout(widths = c(1, 1.2))
mid  <- p_calib  + p_dca     + plot_layout(widths = c(1, 1))
fig  <- top / mid +
        plot_annotation(tag_levels = "a") &
        ukbsci_theme()
save_ukbsci_figure(fig, "Fig01-composite",
                   width_mm = 180, height_mm = 200)
```

`plot_annotation(tag_levels = "a")` gives lowercase `a, b, c, …` panel
labels — common across observational epidemiology.

---

## 11. Common pitfalls

1. **Mixing `theme_classic()` and `theme_minimal()`** across panels of a
   composite. Always end the orchestrating script with `& ukbsci_theme()` so
   every patch inherits the same theme.
2. **Hard-coding hex colors per plot.** Use the `ukbsci_*` palettes and
   keep colour semantics consistent across figures (baseline = blue, exposure
   = red, accent = orange).
3. **Font availability.** `cairo_pdf` with `family = "Arial"` fails if Arial
   is not installed on the RAP node. Fallback: use `family = "sans"` or
   `family = "DejaVu Sans"`, or register a project-installed font with
   `grDevices::loadfonts()` before plotting.
3a. **High-DPI PNG for dense figures.** For complex multi-panel composites,
    heatmaps, or SHAP beeswarm plots with many points, add a 600 dpi PNG
    export alongside the standard 300 dpi one to improve readability in
    supplementary figures. The `save_ukbsci_figure()` helper uses 300 dpi by
    default; pass `dpi = 600` for these cases.
4. **`plot_forest()` reverses row order.** Top row in the data frame
   becomes the *bottom* of the figure (so the manuscript "primary" subgroup
   appears at the top). Adjust order with `dplyr::arrange()` accordingly.
5. **Volcano `effect_col` auto-detection** can pick the wrong column when
   multiple effect-like columns exist. Pass `effect_col` explicitly for the
   manuscript.
6. **Calibration n_bins = 10** can over-smooth in 100k+ cohorts. Increase
   to 20 for large samples.
7. **Survival figures** — generate via `ukbsci-survival`; this skill
   re-themes the result.
8. **Figure source CSVs**: every figure script should write its source
   data to `05-figs/data/Fig##-*.csv`. Reviewers expect it; reproducibility
   demands it.

---

## 12. Key functions

| Function | Returns |
|----------|---------|
| `plot_forest(...)` | ggplot |
| `plot_calibration(data, predicted, observed, n_bins, smooth, conf_int)` | ggplot |
| `plot_regression_volcano(...)` | ggplot (with `plot_data`, `label_data` attrs) |
| (helpers defined locally in this skill) `ukbsci_theme()`, `ukbsci_clinical`, `ukbsci_diverging`, `ukbsci_sequential`, `save_ukbsci_figure()` | utilities |

---

## 13. Related skills

| Skill | When |
|-------|------|
| `ukbsci-regression` | Source of forest / volcano inputs. |
| `ukbsci-survival` | KM plots (uses its own `plot_km_curve`). |
| `ukbsci-ml` | ROC / PR / calibration / DCA / SHAP plots. |
| `ukbsci-propensity` | Love + PS distribution. |
| `ukbsci-mediation` | Mediation forest. |
| `ukbsci-imputation` | Pooled forest + FMI. |
| `ukbsci-proteomics` | GO / KEGG figures. |
| `ukbsci-workflow` | Phase 8 — assemble the final manuscript figure pack. |

---

## 14. References

- [`references/functions.md`](references/functions.md)
- [`references/theme-and-palettes.R`](references/theme-and-palettes.R) —
  copy-paste helpers
- [`references/examples.md`](references/examples.md)
