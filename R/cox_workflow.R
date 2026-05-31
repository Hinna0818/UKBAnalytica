#' Standardize variables using training-set parameters
#'
#' @description
#' Standardize a set of variables in the training data and optionally apply the
#' same centering and scaling parameters to a validation data set. This is useful
#' for omics analyses where all downstream association estimates should be
#' expressed per one training-set standard deviation.
#'
#' @param train_data Training data.
#' @param validation_data Optional validation data.
#' @param variables Character vector of variables to standardize.
#' @param center Logical. If TRUE, subtract the training-set mean.
#' @param scale Logical. If TRUE, divide by the training-set standard deviation.
#'
#' @return A list with `train`, `validation`, and `parameters`.
#' @export
#'
#' @examples
#' dat <- data.frame(x = 1:5, y = c(2, 3, 5, 7, 11))
#' ukb_standardize_by_train(dat, variables = c("x", "y"))$parameters
ukb_standardize_by_train <- function(train_data,
                                     validation_data = NULL,
                                     variables,
                                     center = TRUE,
                                     scale = TRUE) {
  if (!is.data.frame(train_data)) {
    stop("`train_data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.null(validation_data) && !is.data.frame(validation_data)) {
    stop("`validation_data` must be NULL or a data.frame/data.table.", call. = FALSE)
  }
  if (!is.character(variables) || length(variables) == 0L || anyNA(variables)) {
    stop("`variables` must be a non-empty character vector.", call. = FALSE)
  }
  missing_train <- setdiff(variables, names(train_data))
  if (length(missing_train) > 0L) {
    stop("Training data is missing variable(s): ", paste(missing_train, collapse = ", "), call. = FALSE)
  }
  if (!is.null(validation_data)) {
    missing_valid <- setdiff(variables, names(validation_data))
    if (length(missing_valid) > 0L) {
      stop("Validation data is missing variable(s): ", paste(missing_valid, collapse = ", "), call. = FALSE)
    }
  }

  train <- as.data.table(copy(train_data))
  validation <- if (!is.null(validation_data)) as.data.table(copy(validation_data)) else NULL

  params <- data.frame(
    variable = variables,
    center = NA_real_,
    scale = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(variables)) {
    var <- variables[[i]]
    x <- train[[var]]
    if (!is.numeric(x)) {
      stop("Variable `", var, "` must be numeric for standardization.", call. = FALSE)
    }
    mu <- if (isTRUE(center)) mean(x, na.rm = TRUE) else 0
    sig <- if (isTRUE(scale)) sd(x, na.rm = TRUE) else 1
    if (!is.finite(mu)) mu <- 0
    if (!is.finite(sig) || sig == 0) sig <- 1

    train[[var]] <- (train[[var]] - mu) / sig
    if (!is.null(validation)) {
      validation[[var]] <- (validation[[var]] - mu) / sig
    }
    params$center[[i]] <- mu
    params$scale[[i]] <- sig
  }

  list(train = train, validation = validation, parameters = params)
}

#' Standardize variables using existing scaling parameters
#'
#' @description
#' Apply previously estimated centering and scaling parameters to a data set.
#' The parameter table can use either the native output from
#' [ukb_standardize_by_train()] (`variable`, `center`, `scale`) or the legacy
#' long format used by early case-study scripts (`protein`, `statistic`,
#' `value`).
#'
#' @param data Data frame to transform.
#' @param parameters Scaling parameter table.
#' @param variables Optional variables to transform. Defaults to all variables
#'   available in the parameter table.
#'
#' @return A data.table with standardized variables.
#' @export
ukb_scale_with_parameters <- function(data, parameters, variables = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.data.frame(parameters)) {
    stop("`parameters` must be a data.frame or data.table.", call. = FALSE)
  }
  params <- as.data.table(copy(parameters))
  if (all(c("protein", "statistic", "value") %in% names(params))) {
    params <- dcast(params, protein ~ statistic, value.var = "value")
    names(params)[names(params) == "protein"] <- "variable"
    names(params)[names(params) == "mean"] <- "center"
    names(params)[names(params) == "sd"] <- "scale"
  }
  required <- c("variable", "center", "scale")
  missing_cols <- setdiff(required, names(params))
  if (length(missing_cols) > 0L) {
    stop("`parameters` is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (is.null(variables)) {
    variables <- params$variable
  }
  variables <- intersect(variables, params$variable)
  missing_data <- setdiff(variables, names(data))
  if (length(missing_data) > 0L) {
    stop("`data` is missing variable(s): ", paste(missing_data, collapse = ", "), call. = FALSE)
  }
  out <- as.data.table(copy(data))
  params <- params[match(variables, params$variable), , drop = FALSE]
  params$scale[!is.finite(params$scale) | params$scale == 0] <- 1
  params$center[!is.finite(params$center)] <- 0
  for (i in seq_along(variables)) {
    var <- variables[[i]]
    out[[var]] <- (out[[var]] - params$center[[i]]) / params$scale[[i]]
  }
  out
}

#' Compare Cox results between training and validation sets
#'
#' @description
#' Merge two Cox result tables by variable, summarize replication of
#' training-set significant variables in validation, and compute log(HR)
#' correlations.
#'
#' @param train_results Cox result table for the training set.
#' @param validation_results Cox result table for the validation set.
#' @param variable_col Variable column name.
#' @param hr_col Hazard-ratio column name.
#' @param p_col Raw p-value column name.
#' @param train_prefix Prefix for training-set columns in the comparison table.
#' @param validation_prefix Prefix for validation-set columns.
#' @param p_adjust_methods Multiple-testing correction methods to add when
#'   adjusted p-value columns are absent. Defaults to BH and Bonferroni.
#' @param alpha Significance threshold.
#'
#' @return A list with `train_results`, `validation_results`, `comparison`,
#'   `replication_summary`, and `correlation_summary`.
#' @export
ukb_compare_cox_results <- function(train_results,
                                    validation_results,
                                    variable_col = "variable",
                                    hr_col = "HR",
                                    p_col = "pvalue",
                                    train_prefix = "train",
                                    validation_prefix = "validation",
                                    p_adjust_methods = c("BH", "bonferroni"),
                                    alpha = 0.05) {
  train <- .ukb_prepare_cox_table(train_results, variable_col, hr_col, p_col, p_adjust_methods, alpha)
  validation <- .ukb_prepare_cox_table(validation_results, variable_col, hr_col, p_col, p_adjust_methods, alpha)

  train_keep <- .ukb_prefix_cox_columns(train, variable_col, train_prefix)
  valid_keep <- .ukb_prefix_cox_columns(validation, variable_col, validation_prefix)
  comparison <- merge(train_keep, valid_keep, by = variable_col, all.x = TRUE)

  train_log <- paste0(train_prefix, "_logHR")
  valid_log <- paste0(validation_prefix, "_logHR")
  comparison[["same_direction"]] <- sign(comparison[[train_log]]) == sign(comparison[[valid_log]])

  for (method in p_adjust_methods) {
    method_id <- .ukb_p_adjust_id(method)
    valid_sig_col <- paste0(validation_prefix, "_significant_", method_id)
    comparison[[paste0(validation_prefix, "_", method_id)]] <- comparison[[valid_sig_col]]
  }
  comparison[[paste0(validation_prefix, "_nominal")]] <- comparison[[paste0(validation_prefix, "_", p_col)]] < alpha
  comparison <- comparison[order(comparison[[paste0(train_prefix, "_", p_col)]]), , drop = FALSE]

  replication_summary <- .ukb_cox_replication_summary(
    comparison = comparison,
    train_prefix = train_prefix,
    validation_prefix = validation_prefix,
    p_adjust_methods = p_adjust_methods,
    alpha = alpha
  )
  correlation_summary <- .ukb_cox_correlation_summary(
    comparison = comparison,
    train_prefix = train_prefix,
    validation_prefix = validation_prefix,
    p_adjust_methods = p_adjust_methods
  )

  list(
    train_results = train,
    validation_results = validation,
    comparison = comparison,
    replication_summary = replication_summary,
    correlation_summary = correlation_summary
  )
}

#' Compare sensitivity Cox results against a main analysis
#'
#' @description
#' Merge one or more sensitivity-analysis Cox result tables with a main Cox
#' result table, then summarize concordance by sensitivity analysis.
#'
#' @param main_results Main Cox result table.
#' @param sensitivity_results Sensitivity Cox result table containing one row
#'   per variable and sensitivity analysis.
#' @param sensitivity_col Column identifying the sensitivity analysis.
#' @param variable_col Variable column.
#' @param hr_col Hazard-ratio column.
#' @param p_col Raw p-value column.
#' @param main_prefix Prefix for main-analysis columns.
#' @param sensitivity_prefix Prefix for sensitivity-analysis columns.
#' @param p_adjust_methods Multiple-testing correction methods to add if absent.
#' @param alpha Significance threshold.
#'
#' @return A list with standardized result tables, comparison table, and
#'   correlation summary.
#' @export
ukb_compare_sensitivity_cox <- function(main_results,
                                        sensitivity_results,
                                        sensitivity_col = "sensitivity",
                                        variable_col = "variable",
                                        hr_col = "HR",
                                        p_col = "pvalue",
                                        main_prefix = "main",
                                        sensitivity_prefix = "sensitivity",
                                        p_adjust_methods = c("BH", "bonferroni"),
                                        alpha = 0.05) {
  main <- .ukb_prepare_cox_table(main_results, variable_col, hr_col, p_col, p_adjust_methods, alpha)
  sens <- .ukb_prepare_cox_table(sensitivity_results, variable_col, hr_col, p_col, p_adjust_methods, alpha)
  if (!sensitivity_col %in% names(sens)) {
    stop("`sensitivity_results` is missing `", sensitivity_col, "`.", call. = FALSE)
  }

  main_keep <- .ukb_prefix_cox_columns(main, variable_col, main_prefix)
  sens_keep <- as.data.frame(sens)
  keep_cols <- setdiff(names(sens_keep), "dataset")
  sens_keep <- sens_keep[, keep_cols, drop = FALSE]
  prefix_cols <- setdiff(names(sens_keep), c(variable_col, sensitivity_col))
  names(sens_keep)[match(prefix_cols, names(sens_keep))] <- paste0(sensitivity_prefix, "_", prefix_cols)

  comparison <- merge(main_keep, sens_keep, by = variable_col, allow.cartesian = TRUE)
  main_log <- paste0(main_prefix, "_logHR")
  sens_log <- paste0(sensitivity_prefix, "_logHR")
  comparison[["same_direction"]] <- sign(comparison[[main_log]]) == sign(comparison[[sens_log]])

  for (method in p_adjust_methods) {
    method_id <- .ukb_p_adjust_id(method)
    comparison[[paste0(main_prefix, "_", method_id, "_retained")]] <-
      comparison[[paste0(main_prefix, "_significant_", method_id)]] %in% TRUE &
      comparison[[paste0(sensitivity_prefix, "_significant_", method_id)]] %in% TRUE
  }

  correlation_summary <- .ukb_cox_sensitivity_summary(
    comparison = comparison,
    sensitivity_col = sensitivity_col,
    main_prefix = main_prefix,
    sensitivity_prefix = sensitivity_prefix,
    p_adjust_methods = p_adjust_methods
  )

  list(
    main_results = main,
    sensitivity_results = sens,
    comparison = comparison,
    correlation_summary = correlation_summary
  )
}

#' Run Cox models in training and validation sets
#'
#' @description
#' Fit the same multivariable Cox model series in a training set and validation
#' set, optionally standardizing the main variables using training-set
#' parameters, then summarize replication and log(HR) concordance.
#'
#' @param train_data Training data.
#' @param validation_data Validation data.
#' @param main_vars Main variables to evaluate.
#' @param covariates Adjustment covariates.
#' @param endpoint Two-column endpoint passed to [runmulti_cox()].
#' @param standardize_main_vars Logical. If TRUE, standardize `main_vars` using
#'   training-set means and SDs.
#' @param add_protein_annotation Logical. If TRUE, add parsed protein names and
#'   gene symbols for Olink-style protein columns.
#' @param protein_prefix Regular expression prefix removed from protein columns.
#' @param train_label Training-set label.
#' @param validation_label Validation-set label.
#' @param comparison_train_prefix Prefix for training columns in the comparison
#'   table.
#' @param comparison_validation_prefix Prefix for validation columns in the
#'   comparison table.
#' @param p_adjust_methods P-value adjustment methods.
#' @param alpha Significance threshold.
#' @param ... Additional arguments passed to [runmulti_cox()].
#'
#' @return A list containing scaled data, scaling parameters, Cox results, and
#'   comparison summaries.
#' @export
ukb_train_validation_cox <- function(train_data,
                                     validation_data,
                                     main_vars,
                                     covariates,
                                     endpoint,
                                     standardize_main_vars = TRUE,
                                     add_protein_annotation = FALSE,
                                     protein_prefix = "^olink_instance_0[.]",
                                     train_label = "train",
                                     validation_label = "validation",
                                     comparison_train_prefix = "train",
                                     comparison_validation_prefix = "valid",
                                     p_adjust_methods = c("BH", "bonferroni"),
                                     alpha = 0.05,
                                     ...) {
  if (isTRUE(standardize_main_vars)) {
    scaled <- ukb_standardize_by_train(train_data, validation_data, main_vars)
    train_fit_data <- scaled$train
    validation_fit_data <- scaled$validation
    scaling_parameters <- scaled$parameters
  } else {
    train_fit_data <- as.data.table(copy(train_data))
    validation_fit_data <- as.data.table(copy(validation_data))
    scaling_parameters <- NULL
  }

  train_res <- runmulti_cox(
    data = train_fit_data,
    main_var = main_vars,
    covariates = covariates,
    endpoint = endpoint,
    ...
  )
  validation_res <- runmulti_cox(
    data = validation_fit_data,
    main_var = main_vars,
    covariates = covariates,
    endpoint = endpoint,
    ...
  )

  train_res <- as.data.table(train_res)
  validation_res <- as.data.table(validation_res)
  train_res[["dataset"]] <- train_label
  validation_res[["dataset"]] <- validation_label

  cmp <- ukb_compare_cox_results(
    train_results = train_res,
    validation_results = validation_res,
    train_prefix = comparison_train_prefix,
    validation_prefix = comparison_validation_prefix,
    p_adjust_methods = p_adjust_methods,
    alpha = alpha
  )

  train_res <- as.data.table(cmp$train_results)
  validation_res <- as.data.table(cmp$validation_results)
  train_res[["dataset"]] <- train_label
  validation_res[["dataset"]] <- validation_label

  if (isTRUE(add_protein_annotation)) {
    annotation <- ukb_protein_annotation(main_vars, protein_prefix = protein_prefix)
    train_res <- .ukb_add_annotation(train_res, annotation)
    validation_res <- .ukb_add_annotation(validation_res, annotation)
    cmp <- ukb_compare_cox_results(
      train_results = train_res,
      validation_results = validation_res,
      train_prefix = comparison_train_prefix,
      validation_prefix = comparison_validation_prefix,
      p_adjust_methods = p_adjust_methods,
      alpha = alpha
    )
  } else {
    annotation <- NULL
  }

  list(
    train_data = train_fit_data,
    validation_data = validation_fit_data,
    scaling_parameters = scaling_parameters,
    annotation = annotation,
    train_results = train_res,
    validation_results = validation_res,
    comparison = cmp$comparison,
    replication_summary = cmp$replication_summary,
    correlation_summary = cmp$correlation_summary
  )
}

#' Annotate Olink-style protein variables
#'
#' @param variables Protein variable names.
#' @param protein_prefix Regular expression prefix removed from variables.
#' @param drop_unmapped Passed to [protein_to_gene_symbol()].
#'
#' @return A data.frame with variable, protein_clean, gene_symbol, and
#'   mapping_source.
#' @export
ukb_protein_annotation <- function(variables,
                                   protein_prefix = "^olink_instance_0[.]",
                                   drop_unmapped = FALSE) {
  if (!is.character(variables) || length(variables) == 0L || anyNA(variables)) {
    stop("`variables` must be a non-empty character vector.", call. = FALSE)
  }
  manifest <- data.frame(
    variable = variables,
    protein_clean = toupper(sub(protein_prefix, "", variables)),
    stringsAsFactors = FALSE
  )
  map <- protein_to_gene_symbol(manifest$protein_clean, drop_unmapped = drop_unmapped)
  names(map)[names(map) == "protein"] <- "protein_clean"
  merged <- merge(manifest, map, by = "protein_clean", all.x = TRUE)
  split_rows <- split(merged, paste(merged$variable, merged$protein_clean, sep = "\r"))
  out <- do.call(rbind, lapply(split_rows, function(x) {
    collapse_unique <- function(v) {
      v <- unique(v[!is.na(v) & nzchar(v)])
      if (length(v) == 0L) NA_character_ else paste(v, collapse = ";")
    }
    data.frame(
      variable = x$variable[[1]],
      protein_clean = x$protein_clean[[1]],
      gene_symbol = collapse_unique(x$gene_symbol),
      mapping_source = collapse_unique(x$mapping_source),
      stringsAsFactors = FALSE
    )
  }))
  out[order(match(out$variable, variables)), , drop = FALSE]
}

#' Select top Cox associations by hazard ratio
#'
#' @param results Cox result table.
#' @param n_each_direction Number of HR > 1 and HR < 1 rows to keep.
#' @param p_col Adjusted p-value column used for filtering.
#' @param alpha Significance threshold.
#' @param hr_col Hazard-ratio column.
#' @param label_cols Candidate label columns.
#' @param dataset Optional dataset label added to output.
#'
#' @return A data.frame.
#' @export
ukb_top_hr_results <- function(results,
                               n_each_direction = 10,
                               p_col = "p_bonferroni",
                               alpha = 0.05,
                               hr_col = "HR",
                               label_cols = c("gene_symbol", "protein_clean", "variable"),
                               dataset = NULL) {
  d <- as.data.table(copy(results))
  required <- c(hr_col, p_col)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0L) {
    stop("`results` is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- d[is.finite(d[[hr_col]]) & !is.na(d[[p_col]]) & d[[p_col]] < alpha]
  if (nrow(d) == 0L) {
    return(as.data.frame(d))
  }
  label_col <- label_cols[label_cols %in% names(d)][1]
  d[["label"]] <- if (!is.na(label_col)) as.character(d[[label_col]]) else as.character(seq_len(nrow(d)))
  d[["label"]][is.na(d[["label"]]) | !nzchar(d[["label"]])] <- as.character(d[["variable"]][is.na(d[["label"]]) | !nzchar(d[["label"]])])

  order_p <- if ("pvalue" %in% names(d)) "pvalue" else p_col
  high <- d[d[[hr_col]] > 1, , drop = FALSE]
  high <- high[order(-high[[hr_col]], high[[order_p]]), , drop = FALSE]
  low <- d[d[[hr_col]] < 1, , drop = FALSE]
  low <- low[order(low[[hr_col]], low[[order_p]]), , drop = FALSE]
  out <- rbind(
    head(high, n_each_direction),
    head(low, n_each_direction)
  )
  out[["direction"]] <- ifelse(out[[hr_col]] > 1, "HR > 1", "HR < 1")
  if (!is.null(dataset)) {
    out[["dataset"]] <- dataset
  }
  row.names(out) <- NULL
  as.data.frame(out)
}

#' Plot top positive and inverse Cox associations
#'
#' @param top_results A data.frame from [ukb_top_hr_results()] or equivalent.
#' @param facet_col Optional column used for faceting, commonly `"dataset"`.
#' @param hr_col HR column.
#' @param lower_col Lower confidence-limit column.
#' @param upper_col Upper confidence-limit column.
#' @param label_col Label column.
#'
#' @return A ggplot object.
#' @export
plot_top_hr_bars <- function(top_results,
                             facet_col = "dataset",
                             hr_col = "HR",
                             lower_col = "lower95",
                             upper_col = "upper95",
                             label_col = "label") {
  d <- as.data.table(copy(top_results))
  required <- c(hr_col, lower_col, upper_col, label_col)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0L) {
    stop("`top_results` is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d[["signed_log_hr"]] <- log(d[[hr_col]])
  d[["signed_log_lower"]] <- log(d[[lower_col]])
  d[["signed_log_upper"]] <- log(d[[upper_col]])
  d[["gradient_value"]] <- d[["signed_log_hr"]]
  if (!"direction" %in% names(d)) {
    d[["direction"]] <- ifelse(d[[hr_col]] > 1, "HR > 1", "HR < 1")
  }
  if (!is.null(facet_col) && facet_col %in% names(d)) {
    d[[".facet"]] <- d[[facet_col]]
  } else {
    d[[".facet"]] <- "Cox results"
  }
  d[[".y"]] <- paste(d[[".facet"]], d[["direction"]], d[[label_col]], sep = "__")
  d <- d[order(d[[".facet"]], d[["signed_log_hr"]]), , drop = FALSE]
  d[[".y"]] <- factor(d[[".y"]], levels = rev(unique(d[[".y"]])))
  max_abs <- max(abs(d[["gradient_value"]]), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  ggplot(d, aes(x = .data[["signed_log_hr"]], y = .data[[".y"]], fill = .data[["gradient_value"]])) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "#767676") +
    geom_col(width = 0.72, alpha = 0.95) +
    geom_errorbar(
      aes(xmin = .data[["signed_log_lower"]], xmax = .data[["signed_log_upper"]]),
      orientation = "y",
      width = 0.18,
      linewidth = 0.25,
      colour = "#272727"
    ) +
    facet_wrap(vars(.data[[".facet"]]), nrow = 1, scales = "free_y") +
    scale_fill_gradient2(low = "#2F6FA3", mid = "#F7F7F7", high = "#C74732", midpoint = 0, limits = c(-max_abs, max_abs), guide = "none") +
    scale_x_continuous(breaks = log(c(0.75, 0.85, 1, 1.2, 1.5)), labels = c("0.75", "0.85", "1", "1.2", "1.5"), expand = expansion(mult = c(0.04, 0.08))) +
    scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
    labs(x = "HR", y = NULL) +
    theme_classic(base_size = 7) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 6.5, face = "bold"),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )
}

#' Plot training-validation Cox log(HR) concordance
#'
#' @param comparison Comparison table from [ukb_compare_cox_results()].
#' @param train_loghr_col Training log(HR) column.
#' @param validation_loghr_col Validation log(HR) column.
#' @param highlight_col Optional logical column used to highlight proteins.
#' @param highlight_label Highlight legend label.
#'
#' @return A ggplot object.
#' @export
plot_cox_loghr_correlation <- function(comparison,
                                       train_loghr_col = "train_logHR",
                                       validation_loghr_col = "validation_logHR",
                                       highlight_col = "train_significant_bonferroni",
                                       highlight_label = "Train Bonferroni significant") {
  d <- as.data.table(copy(comparison))
  if (!validation_loghr_col %in% names(d) && "valid_logHR" %in% names(d)) {
    validation_loghr_col <- "valid_logHR"
  }
  required <- c(train_loghr_col, validation_loghr_col)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0L) {
    stop("`comparison` is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- d[is.finite(d[[train_loghr_col]]) & is.finite(d[[validation_loghr_col]]), , drop = FALSE]
  if (nrow(d) < 3L) {
    stop("At least three finite paired log(HR) values are required.", call. = FALSE)
  }
  d[[".group"]] <- "Other variables"
  if (!is.null(highlight_col) && highlight_col %in% names(d)) {
    d[[".group"]] <- ifelse(isTRUE(d[[highlight_col]]) | d[[highlight_col]] %in% TRUE, highlight_label, "Other variables")
  }
  d[[".group"]] <- factor(d[[".group"]], levels = c("Other variables", highlight_label))

  pearson <- cor.test(d[[train_loghr_col]], d[[validation_loghr_col]], method = "pearson")
  spearman <- cor.test(d[[train_loghr_col]], d[[validation_loghr_col]], method = "spearman")
  label_text <- sprintf("r = %.3f\nrho = %.3f\nn = %d", unname(pearson$estimate), unname(spearman$estimate), nrow(d))
  axis_lim <- range(c(d[[train_loghr_col]], d[[validation_loghr_col]]), na.rm = TRUE)
  pad <- diff(axis_lim) * 0.06
  if (is.finite(pad) && pad > 0) axis_lim <- axis_lim + c(-pad, pad)

  ggplot(d, aes(x = .data[[train_loghr_col]], y = .data[[validation_loghr_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#6F6F6F") +
    geom_point(aes(colour = .data[[".group"]]), size = 0.8, alpha = 0.7, stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.48, colour = "#D55E00") +
    annotate("text", x = axis_lim[[1]], y = axis_lim[[2]], label = label_text, hjust = 0, vjust = 1, size = 2.1, lineheight = 0.95) +
    scale_colour_manual(values = c("Other variables" = "#B8B8B8", highlight_label = "#0072B2"), drop = FALSE) +
    coord_equal(xlim = axis_lim, ylim = axis_lim) +
    labs(x = "Training-set log(HR)", y = "Validation-set log(HR)", colour = NULL) +
    theme_classic(base_size = 7) +
    theme(legend.position = "bottom")
}

#' Plot sensitivity-analysis Cox log(HR) concordance
#'
#' @param comparison Comparison table from [ukb_compare_sensitivity_cox()].
#' @param sensitivity_col Sensitivity-analysis label column.
#' @param main_loghr_col Main-analysis log(HR) column.
#' @param sensitivity_loghr_col Sensitivity-analysis log(HR) column.
#' @param highlight_col Optional logical column used to highlight variables.
#' @param highlight_label Highlight legend label.
#' @param nrow,ncol Facet layout.
#'
#' @return A ggplot object.
#' @export
plot_cox_sensitivity_correlation <- function(comparison,
                                             sensitivity_col = "sensitivity",
                                             main_loghr_col = "main_logHR",
                                             sensitivity_loghr_col = "sensitivity_logHR",
                                             highlight_col = "main_significant_bonferroni",
                                             highlight_label = "Main Bonferroni significant",
                                             nrow = NULL,
                                             ncol = NULL) {
  d <- as.data.table(copy(comparison))
  required <- c(sensitivity_col, main_loghr_col, sensitivity_loghr_col)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0L) {
    stop("`comparison` is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  d <- d[is.finite(d[[main_loghr_col]]) & is.finite(d[[sensitivity_loghr_col]]), , drop = FALSE]
  if (nrow(d) < 3L) {
    stop("At least three finite paired log(HR) values are required.", call. = FALSE)
  }
  d[[".group"]] <- "Other variables"
  if (!is.null(highlight_col) && highlight_col %in% names(d)) {
    d[[".group"]] <- ifelse(d[[highlight_col]] %in% TRUE, highlight_label, "Other variables")
  }
  d[[".group"]] <- factor(d[[".group"]], levels = c("Other variables", highlight_label))

  label_dt <- .ukb_sensitivity_plot_labels(d, sensitivity_col, main_loghr_col, sensitivity_loghr_col)
  d <- merge(d, label_dt, by = sensitivity_col, all.x = TRUE)
  axis_lim <- range(c(d[[main_loghr_col]], d[[sensitivity_loghr_col]]), na.rm = TRUE)
  pad <- diff(axis_lim) * 0.05
  if (is.finite(pad) && pad > 0) axis_lim <- axis_lim + c(-pad, pad)

  ggplot(d, aes(x = .data[[main_loghr_col]], y = .data[[sensitivity_loghr_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#6F6F6F") +
    geom_point(aes(colour = .data[[".group"]]), size = 0.75, alpha = 0.68, stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.45, colour = "#D55E00") +
    geom_text(
      data = label_dt,
      aes(x = axis_lim[[1]], y = axis_lim[[2]], label = .data[["label_text"]]),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.0,
      colour = "black",
      lineheight = 0.95
    ) +
    facet_wrap(vars(.data[[sensitivity_col]]), nrow = nrow, ncol = ncol) +
    scale_colour_manual(values = c("Other variables" = "#B8B8B8", highlight_label = "#0072B2"), drop = FALSE) +
    coord_equal(xlim = axis_lim, ylim = axis_lim) +
    labs(x = "Main analysis log(HR)", y = "Sensitivity analysis log(HR)", colour = NULL) +
    theme_classic(base_size = 7) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 6.5, face = "bold"),
      legend.position = "bottom",
      panel.spacing.x = unit(8, "pt"),
      panel.spacing.y = unit(8, "pt")
    )
}

.ukb_prepare_cox_table <- function(results, variable_col, hr_col, p_col, p_adjust_methods, alpha) {
  d <- as.data.table(copy(results))
  required <- c(variable_col, hr_col, p_col)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0L) {
    stop("Cox result table is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!"logHR" %in% names(d)) {
    d[["logHR"]] <- log(d[[hr_col]])
  }
  for (method in p_adjust_methods) {
    method_id <- .ukb_p_adjust_id(method)
    p_adj_col <- paste0("p_", method_id)
    sig_col <- paste0("significant_", method_id)
    if (!p_adj_col %in% names(d)) {
      d[[p_adj_col]] <- p.adjust(d[[p_col]], method = method)
    }
    if (!sig_col %in% names(d)) {
      d[[sig_col]] <- d[[p_adj_col]] < alpha
    }
  }
  d[order(d[[p_col]]), , drop = FALSE]
}

.ukb_prefix_cox_columns <- function(d, variable_col, prefix) {
  keep <- setdiff(names(d), "dataset")
  out <- as.data.frame(d)[, keep, drop = FALSE]
  names(out)[names(out) != variable_col] <- paste0(prefix, "_", names(out)[names(out) != variable_col])
  out
}

.ukb_p_adjust_id <- function(method) {
  method <- tolower(method)
  out <- method
  out[method %in% c("bh", "fdr")] <- "bh"
  out
}

.ukb_cox_replication_summary <- function(comparison, train_prefix, validation_prefix, p_adjust_methods, alpha) {
  rows <- list()
  for (method in p_adjust_methods) {
    method_id <- .ukb_p_adjust_id(method)
    train_col <- paste0(train_prefix, "_significant_", method_id)
    d <- comparison[comparison[[train_col]] %in% TRUE & !is.na(comparison[[train_col]]), , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      train_threshold = paste0(method_id, "_", alpha),
      n_train_significant = nrow(d),
      n_validation_nominal = sum(d[[paste0(validation_prefix, "_nominal")]], na.rm = TRUE),
      pct_validation_nominal = .ukb_pct(sum(d[[paste0(validation_prefix, "_nominal")]], na.rm = TRUE), nrow(d)),
      n_validation_bh = if ("bh" %in% .ukb_p_adjust_id(p_adjust_methods)) sum(d[[paste0(validation_prefix, "_bh")]], na.rm = TRUE) else NA_integer_,
      pct_validation_bh = if ("bh" %in% .ukb_p_adjust_id(p_adjust_methods)) .ukb_pct(sum(d[[paste0(validation_prefix, "_bh")]], na.rm = TRUE), nrow(d)) else NA_real_,
      n_validation_bonferroni = if ("bonferroni" %in% .ukb_p_adjust_id(p_adjust_methods)) sum(d[[paste0(validation_prefix, "_bonferroni")]], na.rm = TRUE) else NA_integer_,
      pct_validation_bonferroni = if ("bonferroni" %in% .ukb_p_adjust_id(p_adjust_methods)) .ukb_pct(sum(d[[paste0(validation_prefix, "_bonferroni")]], na.rm = TRUE), nrow(d)) else NA_real_,
      n_same_direction = sum(d[["same_direction"]], na.rm = TRUE),
      pct_same_direction = .ukb_pct(sum(d[["same_direction"]], na.rm = TRUE), nrow(d)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.ukb_pct <- function(num, den) {
  if (den > 0) round(100 * num / den, 2) else NA_real_
}

.ukb_cox_correlation_summary <- function(comparison, train_prefix, validation_prefix, p_adjust_methods) {
  rows <- list(.ukb_one_cox_cor_row(comparison, "all_variables", train_prefix, validation_prefix))
  for (method in p_adjust_methods) {
    method_id <- .ukb_p_adjust_id(method)
    train_col <- paste0(train_prefix, "_significant_", method_id)
    d <- comparison[comparison[[train_col]] %in% TRUE & !is.na(comparison[[train_col]]), , drop = FALSE]
    rows[[length(rows) + 1L]] <- .ukb_one_cox_cor_row(d, paste0(train_prefix, "_", method_id, "_", "significant"), train_prefix, validation_prefix)
  }
  do.call(rbind, rows)
}

.ukb_one_cox_cor_row <- function(d, set_name, train_prefix, validation_prefix) {
  x <- d[[paste0(train_prefix, "_logHR")]]
  y <- d[[paste0(validation_prefix, "_logHR")]]
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L) {
    return(data.frame(protein_set = set_name, n = length(x), method = c("pearson", "spearman"), estimate = NA_real_, pvalue = NA_real_))
  }
  pearson <- cor.test(x, y, method = "pearson")
  spearman <- cor.test(x, y, method = "spearman")
  data.frame(
    protein_set = set_name,
    n = length(x),
    method = c("pearson", "spearman"),
    estimate = c(unname(pearson$estimate), unname(spearman$estimate)),
    pvalue = c(pearson$p.value, spearman$p.value),
    stringsAsFactors = FALSE
  )
}

.ukb_cox_sensitivity_summary <- function(comparison, sensitivity_col, main_prefix, sensitivity_prefix, p_adjust_methods) {
  pieces <- split(comparison, comparison[[sensitivity_col]])
  rows <- lapply(names(pieces), function(sensitivity) {
    d <- pieces[[sensitivity]]
    x <- d[[paste0(main_prefix, "_logHR")]]
    y <- d[[paste0(sensitivity_prefix, "_logHR")]]
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    same <- d[["same_direction"]][keep]
    out <- data.frame(
      sensitivity = sensitivity,
      n_variables = length(x),
      n_same_direction = sum(same, na.rm = TRUE),
      pct_same_direction = .ukb_pct(sum(same, na.rm = TRUE), length(x)),
      pearson_r = NA_real_,
      pearson_pvalue = NA_real_,
      spearman_rho = NA_real_,
      spearman_pvalue = NA_real_,
      stringsAsFactors = FALSE
    )
    if (length(x) >= 3L) {
      pearson <- cor.test(x, y, method = "pearson")
      spearman <- cor.test(x, y, method = "spearman")
      out$pearson_r <- unname(pearson$estimate)
      out$pearson_pvalue <- pearson$p.value
      out$spearman_rho <- unname(spearman$estimate)
      out$spearman_pvalue <- spearman$p.value
    }
    for (method in p_adjust_methods) {
      method_id <- .ukb_p_adjust_id(method)
      main_col <- paste0(main_prefix, "_significant_", method_id)
      retain_col <- paste0(main_prefix, "_", method_id, "_retained")
      den <- sum(d[[main_col]] %in% TRUE, na.rm = TRUE)
      num <- sum(d[[retain_col]] %in% TRUE, na.rm = TRUE)
      out[[paste0(main_prefix, "_", method_id, "_n")]] <- den
      out[[paste0(main_prefix, "_", method_id, "_retained_n")]] <- num
      out[[paste0(main_prefix, "_", method_id, "_retained_pct")]] <- .ukb_pct(num, den)
    }
    out
  })
  do.call(rbind, rows)
}

.ukb_sensitivity_plot_labels <- function(d, sensitivity_col, main_loghr_col, sensitivity_loghr_col) {
  pieces <- split(d, d[[sensitivity_col]])
  rows <- lapply(names(pieces), function(sensitivity) {
    x <- pieces[[sensitivity]][[main_loghr_col]]
    y <- pieces[[sensitivity]][[sensitivity_loghr_col]]
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(x) >= 3L) {
      pearson <- cor.test(x, y, method = "pearson")
      spearman <- cor.test(x, y, method = "spearman")
      label <- sprintf("Pearson r = %.3f\nSpearman rho = %.3f\nn = %d", unname(pearson$estimate), unname(spearman$estimate), length(x))
    } else {
      label <- sprintf("n = %d", length(x))
    }
    data.frame(sensitivity = sensitivity, label_text = label, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  names(out)[names(out) == "sensitivity"] <- sensitivity_col
  out
}

.ukb_add_annotation <- function(results, annotation) {
  merge(as.data.frame(results), annotation, by = "variable", all.x = TRUE)
}
