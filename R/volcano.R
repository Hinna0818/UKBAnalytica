#' Plot a volcano-style regression summary
#'
#' @description
#' Create a volcano-style plot from regression summary results such as
#' [runmulti_cox()] or [runmulti_logit()]. The x-axis shows the supplied effect
#' estimate column (for example `HR` or `OR`), and the y-axis shows
#' `-log10(P)`. Points can be highlighted by an adjusted p-value column, while
#' labels are selected from the largest and smallest highlighted effects.
#'
#' @param data A data.frame containing regression results.
#' @param effect_col Character. Column containing the effect estimate to plot on
#'   the x-axis. If `NULL`, the function uses `HR`, then `OR`, then `estimate`
#'   when available.
#' @param p_col Character. Column containing raw p-values. Default `"pvalue"`.
#' @param adjusted_p_col Optional character. Column used for highlighting
#'   significant points, such as `"p_bonferroni"` or `"p_bh"`. If `NULL`,
#'   `p_col` is used.
#' @param label_col Optional character. Column used for point labels. If `NULL`,
#'   the function uses `gene_symbol`, then `protein_clean`, then `variable` when
#'   available.
#' @param significance_cutoff Numeric cutoff applied to `adjusted_p_col`.
#'   Default `0.05`.
#' @param top_n_label_each Integer. Number of highlighted proteins to label from
#'   each direction. Direction is defined relative to `null_effect`.
#' @param null_effect Numeric null effect. Use `1` for ratio estimates such as
#'   HR/OR and `0` for beta estimates. Default `1`.
#' @param x_lab,y_lab Axis labels. If `NULL`, defaults are generated.
#' @param x_limits,y_limits Optional numeric vectors of length 2 for axis
#'   limits.
#' @param point_size Numeric point size. Default `1.05`.
#' @param label_size Numeric label size. Default `2`.
#' @param colors Named character vector for groups `neutral`, `lower`, and
#'   `higher`.
#' @param show_cutoff Logical. Whether to draw a horizontal significance cutoff
#'   line. Default `TRUE`.
#'
#' @return A `ggplot2` object with attributes `plot_data` and `label_data`.
#'
#'
#' @export
plot_regression_volcano <- function(data,
                                    effect_col = NULL,
                                    p_col = "pvalue",
                                    adjusted_p_col = NULL,
                                    label_col = NULL,
                                    significance_cutoff = 0.05,
                                    top_n_label_each = 5,
                                    null_effect = 1,
                                    x_lab = NULL,
                                    y_lab = NULL,
                                    x_limits = NULL,
                                    y_limits = NULL,
                                    point_size = 1.05,
                                    label_size = 2,
                                    colors = c(
                                      neutral = "#D8D8D8",
                                      lower = "#2F6FA3",
                                      higher = "#C74732"
                                    ),
                                    show_cutoff = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  if (is.null(effect_col)) {
    candidates <- c("HR", "OR", "estimate", "effect")
    effect_col <- candidates[candidates %in% names(data)][1]
    if (is.na(effect_col)) {
      stop("Could not infer `effect_col`. Please provide an effect estimate column.")
    }
  }
  if (!effect_col %in% names(data)) {
    stop("`effect_col` not found in `data`: ", effect_col)
  }
  if (!p_col %in% names(data)) {
    stop("`p_col` not found in `data`: ", p_col)
  }

  if (is.null(adjusted_p_col)) {
    adjusted_p_col <- p_col
  }
  if (!adjusted_p_col %in% names(data)) {
    stop("`adjusted_p_col` not found in `data`: ", adjusted_p_col)
  }

  if (is.null(label_col)) {
    candidates <- c("gene_symbol", "protein_clean", "variable", "term")
    label_col <- candidates[candidates %in% names(data)][1]
  }
  if (is.na(label_col) || is.null(label_col) || !label_col %in% names(data)) {
    label_col <- NULL
  }

  plot_data <- as.data.frame(data, stringsAsFactors = FALSE)
  plot_data$effect <- as.numeric(data[[effect_col]])
  plot_data$pvalue <- as.numeric(data[[p_col]])
  plot_data$adjusted_p <- as.numeric(data[[adjusted_p_col]])
  plot_data$label <- if (!is.null(label_col)) as.character(data[[label_col]]) else NA_character_
  if ("variable" %in% names(data)) {
    plot_data$variable <- as.character(data[["variable"]])
  } else {
    plot_data$variable <- seq_len(nrow(plot_data))
  }

  plot_data <- plot_data[is.finite(plot_data$effect) & is.finite(plot_data$pvalue), , drop = FALSE]
  if (!nrow(plot_data)) {
    stop("No finite effect and p-value rows available for plotting.")
  }

  plot_data$pvalue <- pmax(plot_data$pvalue, .Machine$double.xmin)
  plot_data$neg_log10_p <- -log10(plot_data$pvalue)
  plot_data$significant <- plot_data$adjusted_p < significance_cutoff
  plot_data$direction <- ifelse(plot_data$effect >= null_effect, "higher", "lower")
  plot_data$group <- ifelse(
    !plot_data$significant,
    "neutral",
    ifelse(plot_data$effect >= null_effect, "higher", "lower")
  )
  plot_data$group <- factor(plot_data$group, levels = c("neutral", "lower", "higher"))

  higher <- plot_data[plot_data$significant & plot_data$effect >= null_effect, , drop = FALSE]
  higher <- higher[order(-higher$effect, higher$pvalue), , drop = FALSE]
  lower <- plot_data[plot_data$significant & plot_data$effect < null_effect, , drop = FALSE]
  lower <- lower[order(lower$effect, lower$pvalue), , drop = FALSE]
  label_data <- unique(
    rbind(
      head(higher, top_n_label_each),
      head(lower, top_n_label_each)
    )
  )
  label_data <- label_data[!is.na(label_data$label) & nzchar(label_data$label), , drop = FALSE]

  if (is.null(x_lab)) {
    x_lab <- effect_col
  }
  if (is.null(y_lab)) {
    y_lab <- expression(-log[10](italic(P)))
  }

  if (is.null(x_limits)) {
    x_range <- range(plot_data$effect, na.rm = TRUE)
    x_pad <- diff(x_range) * 0.08
    if (!is.finite(x_pad) || x_pad == 0) {
      x_pad <- 0.1
    }
    x_limits <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  }
  if (is.null(y_limits)) {
    y_limits <- c(0, max(plot_data$neg_log10_p, na.rm = TRUE) * 1.04)
  }

  p <- ggplot(plot_data, aes(x = .data$effect, y = .data$neg_log10_p)) +
    geom_vline(xintercept = null_effect, linewidth = 0.25, colour = "#767676") +
    geom_point(
      aes(colour = .data$group),
      size = point_size,
      alpha = 0.82,
      stroke = 0
    ) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_x_continuous(
      limits = x_limits,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = y_limits,
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = 7) +
    theme(
      legend.position = "none",
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      panel.grid = element_blank()
    )

  if (show_cutoff) {
    p <- p + geom_hline(
      yintercept = -log10(significance_cutoff),
      linewidth = 0.28,
      linetype = "dashed",
      colour = "#767676"
    )
  }

  if (nrow(label_data) > 0) {
    p <- p + geom_text(
      data = label_data,
      aes(label = .data$label),
      size = label_size,
      check_overlap = TRUE,
      vjust = -0.5
    )
  }

  attr(p, "plot_data") <- plot_data
  attr(p, "label_data") <- label_data
  p
}
