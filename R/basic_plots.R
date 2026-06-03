#' Plot a publication-style heatmap
#'
#' @description
#' Draw a compact heatmap from long-format data. The function uses string
#' column names and `.data` pronouns internally, which makes it suitable for
#' scripted package workflows and CRAN checks.
#'
#' @param data A data.frame.
#' @param x Character column name for the x axis.
#' @param y Character column name for the y axis.
#' @param fill Character column name for the heatmap value.
#' @param label Optional character column name for tile labels.
#' @param show_values Logical. If TRUE, show values or `label` on tiles.
#' @param low Low-end color for the diverging scale.
#' @param mid Midpoint color for the diverging scale.
#' @param high High-end color for the diverging scale.
#' @param midpoint Midpoint for the diverging scale.
#' @param title Optional title. If NULL, no title is shown.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param fill_lab Optional fill legend label.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 labs
#'   theme element_blank element_text theme_classic element_line margin
#' @importFrom rlang .data
#' @export
plot_heatmap <- function(data,
                         x,
                         y,
                         fill,
                         label = NULL,
                         show_values = FALSE,
                         low = "#2F6FA3",
                         mid = "#F7F7F7",
                         high = "#C74732",
                         midpoint = 0,
                         title = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         fill_lab = NULL,
                         base_size = 7) {
  .ukb_check_plot_data(data, c(x, y, fill))
  d <- as.data.table(copy(data))
  if (!is.null(label)) {
    .ukb_check_plot_data(d, label)
    d[[".label"]] <- as.character(d[[label]])
  } else {
    d[[".label"]] <- sprintf("%.2f", d[[fill]])
  }

  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_tile(color = "white", linewidth = 0.28) +
    scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = midpoint,
      name = fill_lab
    ) +
    labs(title = title, x = xlab, y = ylab) +
    .ukb_theme_publication(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.border = element_blank(),
      legend.position = "right"
    )

  if (isTRUE(show_values)) {
    p <- p + geom_text(aes(label = .data[[".label"]]), size = base_size * 0.32, colour = "#202020")
  }
  p
}

#' Plot a publication-style stacked bar chart
#'
#' @description
#' Summarize observations by `x` and `fill`, then draw either proportional or
#' count-based stacked bars.
#'
#' @param data A data.frame.
#' @param x Character column name for bar groups.
#' @param fill Character column name for stack groups.
#' @param weight Optional numeric column name for weighted summaries.
#' @param position Either `"fill"` for proportions or `"stack"` for counts.
#' @param palette Optional vector of fill colors.
#' @param title Optional title. If NULL, no title is shown.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param legend_title Optional legend title.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col scale_y_continuous scale_fill_manual
#'   labs theme expansion theme_classic element_line element_text element_blank
#'   margin
#' @importFrom rlang .data
#' @importFrom stats ave
#' @export
plot_stacked_bar <- function(data,
                             x,
                             fill,
                             weight = NULL,
                             position = c("fill", "stack"),
                             palette = NULL,
                             title = NULL,
                             xlab = NULL,
                             ylab = NULL,
                             legend_title = NULL,
                             base_size = 7) {
  position <- match.arg(position)
  .ukb_check_plot_data(data, c(x, fill))
  if (!is.null(weight)) {
    .ukb_check_plot_data(data, weight)
  }
  d <- as.data.table(copy(data))
  keep <- !is.na(d[[x]]) & !is.na(d[[fill]])
  d <- d[keep]
  if (!is.null(weight)) {
    keep_weight <- !is.na(d[[weight]])
    d <- d[keep_weight]
    d[[".ukb_weight"]] <- d[[weight]]
    summary_dt <- d[, .(n = sum(get(".ukb_weight"), na.rm = TRUE)), by = c(x, fill)]
  } else {
    summary_dt <- d[, .(n = .N), by = c(x, fill)]
  }
  summary_dt[["total"]] <- ave(summary_dt[["n"]], summary_dt[[x]], FUN = sum)
  summary_dt[["proportion"]] <- summary_dt[["n"]] / summary_dt[["total"]]
  summary_dt[[".y"]] <- if (position == "fill") summary_dt[["proportion"]] else summary_dt[["n"]]
  if (is.null(ylab)) {
    ylab <- if (position == "fill") "Proportion" else "Count"
  }

  p <- ggplot(summary_dt, aes(x = .data[[x]], y = .data[[".y"]], fill = .data[[fill]])) +
    geom_col(width = 0.72, color = "white", linewidth = 0.18) +
    labs(title = title, x = xlab, y = ylab, fill = legend_title) +
    .ukb_theme_publication(base_size = base_size) +
    theme(legend.position = "right")
  if (position == "fill") {
    p <- p + scale_y_continuous(labels = function(z) paste0(round(z * 100), "%"), expand = expansion(mult = c(0, 0.03)))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.04)))
  }
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  } else {
    p <- p + scale_fill_manual(values = .ukb_palette_discrete(length(unique(summary_dt[[fill]]))))
  }
  p
}

#' Plot a publication-style violin plot
#'
#' @description
#' Draw grouped distributions using violin layers with optional boxplot overlay.
#'
#' @param data A data.frame.
#' @param x Character column name for groups.
#' @param y Character numeric column name.
#' @param fill Optional fill grouping column. Defaults to `x`.
#' @param palette Optional vector of fill colors.
#' @param add_boxplot Logical. Overlay a narrow boxplot.
#' @param add_points Logical. Overlay jittered observations.
#' @param title Optional title. If NULL, no title is shown.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_jitter
#'   scale_fill_manual labs theme theme_classic element_line element_text
#'   element_blank margin
#' @importFrom rlang .data
#' @export
plot_violin <- function(data,
                        x,
                        y,
                        fill = NULL,
                        palette = NULL,
                        add_boxplot = TRUE,
                        add_points = FALSE,
                        title = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        base_size = 7) {
  fill <- if (is.null(fill)) x else fill
  .ukb_check_plot_data(data, c(x, y, fill))
  d <- as.data.table(copy(data))
  keep <- !is.na(d[[x]]) & !is.na(d[[y]]) & !is.na(d[[fill]])
  d <- d[keep]

  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_violin(width = 0.82, alpha = 0.82, linewidth = 0.25, color = "#333333", trim = FALSE) +
    labs(title = title, x = xlab, y = ylab, fill = NULL) +
    .ukb_theme_publication(base_size = base_size) +
    theme(legend.position = if (identical(fill, x)) "none" else "right")
  if (isTRUE(add_boxplot)) {
    p <- p + geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25, fill = "white", alpha = 0.75)
  }
  if (isTRUE(add_points)) {
    p <- p + geom_jitter(width = 0.08, height = 0, size = 0.55, alpha = 0.35, color = "#222222")
  }
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  } else {
    p <- p + scale_fill_manual(values = .ukb_palette_discrete(length(unique(d[[fill]]))))
  }
  p
}

#' Plot a publication-style scatter plot
#'
#' @description
#' Draw a scatter plot with optional color grouping, linear smooth, and reference
#' line. This is intended for compact association or validation panels.
#'
#' @param data A data.frame.
#' @param x Character numeric column name for the x axis.
#' @param y Character numeric column name for the y axis.
#' @param color Optional grouping column for point colors.
#' @param palette Optional vector of colors.
#' @param add_smooth Logical. Add a linear smooth line.
#' @param add_identity Logical. Add a dashed y = x reference line.
#' @param alpha Point alpha.
#' @param point_size Point size.
#' @param title Optional title. If NULL, no title is shown.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_smooth
#'   scale_color_manual labs theme_classic theme element_line element_text
#'   element_blank margin
#' @importFrom rlang .data
#' @export
plot_scatter <- function(data,
                         x,
                         y,
                         color = NULL,
                         palette = NULL,
                         add_smooth = TRUE,
                         add_identity = FALSE,
                         alpha = 0.72,
                         point_size = 1.2,
                         title = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         base_size = 7) {
  required <- c(x, y, color)
  required <- required[!is.null(required)]
  .ukb_check_plot_data(data, required)
  d <- as.data.table(copy(data))
  keep <- !is.na(d[[x]]) & !is.na(d[[y]])
  d <- d[keep]
  if (!is.null(color)) {
    keep_color <- !is.na(d[[color]])
    d <- d[keep_color]
  }

  p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]]))
  if (is.null(color)) {
    p <- p + geom_point(size = point_size, alpha = alpha, color = "#2F6FA3", stroke = 0)
  } else {
    p <- p + geom_point(aes(color = .data[[color]]), size = point_size, alpha = alpha, stroke = 0)
  }
  if (isTRUE(add_identity)) {
    p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.32, color = "#6F6F6F")
  }
  if (isTRUE(add_smooth)) {
    p <- p + geom_smooth(method = "lm", se = TRUE, linewidth = 0.42, color = "#C74732", fill = "#C74732", alpha = 0.14)
  }
  p <- p +
    labs(title = title, x = xlab, y = ylab, color = NULL) +
    .ukb_theme_publication(base_size = base_size)
  if (!is.null(color)) {
    if (!is.null(palette)) {
      p <- p + scale_color_manual(values = palette)
    } else {
      p <- p + scale_color_manual(values = .ukb_palette_discrete(length(unique(d[[color]]))))
    }
  }
  p
}

.ukb_check_plot_data <- function(data, vars) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  vars <- vars[!is.na(vars) & nzchar(vars)]
  missing <- setdiff(vars, names(data))
  if (length(missing) > 0L) {
    stop("`data` is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.ukb_theme_publication <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.28, colour = "black"),
      axis.text = element_text(size = base_size * 0.88, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      plot.title = element_text(size = base_size * 1.02, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size * 0.88),
      legend.text = element_text(size = base_size * 0.82),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

.ukb_palette_discrete <- function(n) {
  base <- c("#2F6FA3", "#C74732", "#4D8B57", "#7A5DA8", "#D69C2F", "#5B8C8F", "#8C6D62", "#6F6F6F")
  rep(base, length.out = n)
}
