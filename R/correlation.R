#' Calculate correlation between variables
#' @description Before the formal regression analysis, it can be useful to check the correlation between variables. 
#' This function calculates the correlation matrix for a set of specified variables, which can help identify potential multicollinearity issues or inform variable selection.
#' @param df A data.frame or data.table containing the variables of interest.
#' @param vars A character vector of column names for which to calculate the correlation matrix.
#' @param method The method to use for calculating correlation. Options are "pearson", "spearman", or "kendall". Default is "pearson".
#' @param threshold Numeric value between 0 and 1. If specified, the variables with absolute correlation above this threshold will be highlighted in the output. Default is 0.7.
#' @importFrom stats cor
#' @return A correlation matrix of the specified variables.
#' @export
#' @name run_correlation

# Suppress R CMD check NOTEs for ggplot2 NSE
if(getRversion() >= "2.15.1") utils::globalVariables(c("Var1", "Var2", "Correlation", "label"))

run_correlation <- function(df,
                            vars,
                            method = "pearson",
                            threshold = 0.7) {
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or data.table.")
  }
  
  if (missing(vars) || length(vars) == 0) {
    stop("`vars` must be a non-empty character vector.")
  }
  
  if (!all(vars %in% names(df))) {
    missing_vars <- setdiff(vars, names(df))
    stop(sprintf("The following variables are not in the data frame: %s", paste(missing_vars, collapse = ", ")))
  }
  
  # Subset the data frame to the specified variables
  # Use as.data.frame to ensure compatibility with both data.frame and data.table
  df_subset <- as.data.frame(df)[, vars, drop = FALSE]
  
  # Calculate the correlation matrix
  corr_matrix <- cor(df_subset, use = "pairwise.complete.obs", method = method)

  # Highlight correlations above the threshold
  if (!is.null(threshold) && is.numeric(threshold) && threshold > 0){
    outliers <- which(abs(corr_matrix) > threshold, arr.ind = TRUE)
    if (length(outliers) > 0) {
      cat("Correlations above the threshold of", threshold, ":\n")
      for (i in seq_len(nrow(outliers))) {
        var1 <- rownames(corr_matrix)[outliers[i, 1]]
        var2 <- colnames(corr_matrix)[outliers[i, 2]]
        corr_value <- corr_matrix[outliers[i, 1], outliers[i, 2]]
        cat(sprintf("%s and %s: %.2f\n", var1, var2, corr_value))
      }
    } else {
      cat("No correlations above the threshold of", threshold, "found.\n")
    }
  }
  
  return(corr_matrix)
}


#' Visualize correlation matrix as a heatmap
#'
#' @description
#' Create an annotated heatmap of a correlation matrix with customizable appearance.
#' This helps identify patterns, multicollinearity, and variable relationships visually.
#'
#' @param corr_matrix A numeric correlation matrix (from \code{run_correlation()} or \code{cor()}).
#' @param title Character string. Plot title. Default: \code{"Correlation Matrix"}.
#' @param show_values Logical. If \code{TRUE}, display correlation values on tiles. Default: \code{TRUE}.
#' @param digits Integer. Number of decimal places for correlation values. Default: 2.
#' @param text_size Numeric. Size of text labels on tiles. Default: 3.
#' @param color_low Character. Color for negative correlations. Default: \code{"#3B4CC0"} (blue).
#' @param color_mid Character. Color for zero correlation. Default: \code{"white"}.
#' @param color_high Character. Color for positive correlations. Default: \code{"#B40426"} (red).
#' @param upper_triangle Logical. If \code{TRUE}, show only upper triangle. Default: \code{FALSE}.
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#'   theme_minimal theme element_text element_blank coord_fixed labs
#' @return A \code{ggplot2} object. Can be further customized with \code{ggplot2} functions.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' dt <- data.table(x = rnorm(100), y = rnorm(100), z = rnorm(100))
#' corr_mat <- run_correlation(dt, vars = c("x", "y", "z"))
#'
#' # Basic heatmap
#' plot_correlation(corr_mat)
#'
#' # Upper triangle only without values
#' plot_correlation(
#'   corr_mat,
#'   show_values = FALSE,
#'   upper_triangle = TRUE
#' )
#' }
#'
#' @export

plot_correlation <- function(corr_matrix,
                              title = "Correlation Matrix",
                              show_values = TRUE,
                              digits = 2,
                              text_size = 3,
                              color_low = "#3B4CC0",
                              color_mid = "white",
                              color_high = "#B40426",
                              upper_triangle = FALSE) {
  
  # Input validation
  if (!is.matrix(corr_matrix)) {
    stop("Input must be a correlation matrix.")
  }
  
  if (!is.numeric(corr_matrix)) {
    stop("Correlation matrix must contain numeric values.")
  }
  
  if (nrow(corr_matrix) != ncol(corr_matrix)) {
    stop("Correlation matrix must be square (same number of rows and columns).")
  }
  
  if (is.null(rownames(corr_matrix)) || is.null(colnames(corr_matrix))) {
    rownames(corr_matrix) <- colnames(corr_matrix) <- paste0("Var", seq_len(nrow(corr_matrix)))
  }
  
  # Optional: mask lower triangle
  if (upper_triangle) {
    corr_matrix[lower.tri(corr_matrix)] <- NA
  }
  
  # Convert to long format for ggplot2
  corr_df <- as.data.frame(as.table(corr_matrix))
  names(corr_df) <- c("Var1", "Var2", "Correlation")
  
  # Remove NA values (from masking lower triangle)
  if (upper_triangle) {
    corr_df <- corr_df[!is.na(corr_df$Correlation), ]
  }
  
  # Create base plot
  p <- ggplot(
    corr_df,
    aes(x = Var1, y = Var2, fill = Correlation)
  ) +
    geom_tile(color = "gray90", linewidth = 0.5) +
    scale_fill_gradient2(
      low = color_low,
      high = color_high,
      mid = color_mid,
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      name = "Correlation",
      na.value = "gray95"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    coord_fixed() +
    labs(title = title)
  
  # Add correlation values as text
  if (show_values) {
    corr_df$label <- sprintf(paste0("%.", digits, "f"), corr_df$Correlation)
    
    # Choose text color based on correlation strength
    corr_df$text_color <- ifelse(
      abs(corr_df$Correlation) > 0.5,
      "white",
      "black"
    )
    
    p <- p + ggplot2::geom_text(
      data = corr_df,
      ggplot2::aes(label = label),
      color = corr_df$text_color,
      size = text_size,
      show.legend = FALSE
    )
  }
  
  return(p)
}