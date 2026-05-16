## ukbsci-plot — shared theme & palettes
## Copy verbatim at the top of every figure script.

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

## Categorical palette — clinical / two-arm + accent (neutral-named on purpose)
ukbsci_clinical <- c(
  baseline   = "#2F6FA3",
  exposure   = "#C74732",
  control    = "#7F7F7F",
  accent     = "#F2A93B",
  highlight1 = "#3D8B66",
  highlight2 = "#7A4E9B"
)

## Diverging — heatmaps, SMD, signed effects
ukbsci_diverging <- c("#2166AC", "#67A9CF", "#D1E5F0",
                      "#F7F7F7",
                      "#FDDBC7", "#EF8A62", "#B2182B")

## Sequential — density, continuous gradients
ukbsci_sequential <- c("#FFF7EC", "#FEE8C8", "#FDD49E",
                       "#FDBB84", "#FC8D59", "#EF6548",
                       "#D7301F", "#990000")

## Apply globally
theme_set(ukbsci_theme())

## Dual-format save helper
save_ukbsci_figure <- function(plot, name,
                                width_mm = 180, height_mm = 120,
                                dpi = 300,
                                outdir = "/mnt/project/<area>/05-figs") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  ## SVG (editable)
  svglite::svglite(file.path(outdir, paste0(name, ".svg")),
                   width = w, height = h)
  print(plot); dev.off()

  ## PDF (vector)
  grDevices::cairo_pdf(file.path(outdir, paste0(name, ".pdf")),
                       width = w, height = h, family = "Arial")
  print(plot); dev.off()

  ## PNG (raster)
  ragg::agg_png(file.path(outdir, paste0(name, ".png")),
                width = w, height = h, units = "in", res = dpi)
  print(plot); dev.off()
}
