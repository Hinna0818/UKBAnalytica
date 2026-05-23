#' Build a participant flow table
#'
#' @description
#' Apply sequential inclusion or exclusion rules and record the number of
#' participants retained and removed at each step. Rules can be supplied as
#' one-sided formulas, functions, logical vectors, or character vectors of
#' variables requiring complete-case data.
#'
#' @param data A data.frame or data.table.
#' @param steps A named list of rules. Each rule can be:
#'   a one-sided formula such as `~ !is.na(age)`, a function returning a logical
#'   vector, a logical vector, or a character vector of variable names to retain
#'   complete cases.
#' @param id_col Optional participant identifier column. If supplied, duplicate
#'   non-missing IDs are reported as an error.
#' @param outcome_col Optional 0/1 outcome column used to count events after
#'   each step.
#' @param event_value Value in `outcome_col` indicating an event. Defaults to 1.
#' @param start_label Label for the first row.
#'
#' @return A data.frame with class `ukb_participant_flow`. The kept row index is
#'   stored in `attr(result, "kept_index")`.
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   eid = 1:5,
#'   age = c(50, 60, NA, 55, 70),
#'   status = c(0, 1, 0, 1, 0)
#' )
#' flow <- ukb_participant_flow(
#'   dat,
#'   steps = list("Complete age" = "age"),
#'   id_col = "eid",
#'   outcome_col = "status"
#' )
ukb_participant_flow <- function(data,
                                 steps,
                                 id_col = NULL,
                                 outcome_col = NULL,
                                 event_value = 1,
                                 start_label = "Initial population") {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.list(steps) || length(steps) == 0) {
    stop("`steps` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(steps)) || any(!nzchar(names(steps)))) {
    stop("Every item in `steps` must have a step label.", call. = FALSE)
  }
  if (!is.null(id_col)) {
    .ukb_flow_check_column(data, id_col, "id_col")
    ids <- data[[id_col]]
    duplicated_ids <- ids[!is.na(ids) & duplicated(ids)]
    if (length(duplicated_ids) > 0) {
      stop("`id_col` contains duplicated non-missing participant IDs.", call. = FALSE)
    }
  }
  if (!is.null(outcome_col)) {
    .ukb_flow_check_column(data, outcome_col, "outcome_col")
  }

  n_total <- nrow(data)
  keep <- rep(TRUE, n_total)
  rows <- list(.ukb_flow_row(
    step_order = 0L,
    step = start_label,
    rule_type = "start",
    n_before = n_total,
    n_removed = 0L,
    n_after = n_total,
    n_total = n_total,
    events_after = .ukb_flow_events(data, keep, outcome_col, event_value),
    events_removed = 0L
  ))

  for (i in seq_along(steps)) {
    label <- names(steps)[[i]]
    rule <- steps[[i]]
    before <- keep
    n_before <- sum(before)
    current_data <- data[before, , drop = FALSE]
    rule_keep <- .ukb_eval_flow_rule(rule, current_data, data, before)
    rule_keep[is.na(rule_keep)] <- FALSE

    after <- rep(FALSE, n_total)
    after[which(before)] <- rule_keep
    removed <- before & !after
    keep <- after

    rows[[length(rows) + 1L]] <- .ukb_flow_row(
      step_order = i,
      step = label,
      rule_type = .ukb_flow_rule_type(rule),
      n_before = n_before,
      n_removed = sum(removed),
      n_after = sum(keep),
      n_total = n_total,
      events_after = .ukb_flow_events(data, keep, outcome_col, event_value),
      events_removed = .ukb_flow_events(data, removed, outcome_col, event_value)
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$retained_from_previous <- ifelse(out$n_before > 0, out$n_after / out$n_before, NA_real_)
  out$retained_from_initial <- if (n_total > 0) out$n_after / n_total else NA_real_
  out$event_rate_after <- ifelse(out$n_after > 0, out$events_after / out$n_after, NA_real_)
  attr(out, "kept_index") <- keep
  class(out) <- c("ukb_participant_flow", class(out))
  out
}

#' Plot a participant flow table
#'
#' @param flow A `ukb_participant_flow` object.
#' @param show_removed Logical. If `TRUE`, annotate removals at each step.
#' @param show_events Logical. If `TRUE`, include event counts in labels when
#'   `outcome_col` was supplied to `ukb_participant_flow()`.
#' @param fill Fill color for the retained-participant bars.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' dat <- data.frame(eid = 1:5, age = c(50, 60, NA, 55, 70), status = c(0, 1, 0, 1, 0))
#' flow <- ukb_participant_flow(dat, list("Complete age" = "age"), outcome_col = "status")
#' plot_participant_flow(flow)
plot_participant_flow <- function(flow,
                                  show_removed = TRUE,
                                  show_events = TRUE,
                                  fill = "#2C7FB8") {
  if (!inherits(flow, "ukb_participant_flow")) {
    stop("`flow` must be created by ukb_participant_flow().", call. = FALSE)
  }
  if (!is.logical(show_removed) || length(show_removed) != 1 || is.na(show_removed)) {
    stop("`show_removed` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_events) || length(show_events) != 1 || is.na(show_events)) {
    stop("`show_events` must be TRUE or FALSE.", call. = FALSE)
  }

  plot_data <- flow
  plot_data$step <- factor(plot_data$step, levels = rev(plot_data$step))
  plot_data$label <- paste0("n=", format(plot_data$n_after, big.mark = ","))
  if (isTRUE(show_events) && any(plot_data$events_after > 0)) {
    plot_data$label <- paste0(plot_data$label, "; events=", format(plot_data$events_after, big.mark = ","))
  }
  if (isTRUE(show_removed)) {
    plot_data$removed_label <- ifelse(
      plot_data$step_order == 0,
      "",
      paste0("-", format(plot_data$n_removed, big.mark = ","))
    )
  } else {
    plot_data$removed_label <- ""
  }

  ggplot(plot_data, aes(x = .data$step, y = .data$n_after)) +
    geom_col(width = 0.68, fill = fill, alpha = 0.9) +
    geom_text(aes(label = .data$label), hjust = -0.05, size = 3.4, color = "#243447") +
    geom_text(aes(label = .data$removed_label), hjust = 1.08, size = 3.2, color = "#8A1F1F") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = "Participants retained") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(color = "#243447"),
      axis.text.x = element_text(color = "#243447"),
      axis.title.x = element_text(color = "#243447", margin = margin(t = 8)),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

.ukb_flow_row <- function(step_order, step, rule_type, n_before, n_removed,
                          n_after, n_total, events_after, events_removed) {
  data.frame(
    step_order = as.integer(step_order),
    step = step,
    rule_type = rule_type,
    n_before = as.integer(n_before),
    n_removed = as.integer(n_removed),
    n_after = as.integer(n_after),
    n_total = as.integer(n_total),
    events_removed = as.integer(events_removed),
    events_after = as.integer(events_after),
    stringsAsFactors = FALSE
  )
}

.ukb_flow_check_column <- function(data, col, arg) {
  if (!is.character(col) || length(col) != 1 || is.na(col) || !nzchar(col)) {
    stop(sprintf("`%s` must be a single non-empty column name.", arg), call. = FALSE)
  }
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' was not found in `data`.", col), call. = FALSE)
  }
}

.ukb_flow_events <- function(data, keep, outcome_col, event_value) {
  if (is.null(outcome_col)) {
    return(0L)
  }
  sum(data[[outcome_col]][keep] == event_value, na.rm = TRUE)
}

.ukb_eval_flow_rule <- function(rule, current_data, full_data, before) {
  if (inherits(rule, "formula")) {
    if (length(rule) != 2) {
      stop("Flow formulas must be one-sided, e.g. ~ !is.na(age).", call. = FALSE)
    }
    out <- eval(rule[[2]], envir = current_data, enclos = parent.frame())
  } else if (is.function(rule)) {
    out <- rule(current_data)
  } else if (is.logical(rule)) {
    if (length(rule) == nrow(full_data)) {
      out <- rule[before]
    } else {
      out <- rule
    }
  } else if (is.character(rule)) {
    missing_vars <- setdiff(rule, names(current_data))
    if (length(missing_vars) > 0) {
      stop(
        "Complete-case flow variable(s) not found: ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }
    out <- stats::complete.cases(current_data[, rule, drop = FALSE])
  } else {
    stop("Unsupported flow rule type.", call. = FALSE)
  }

  if (!is.logical(out) || length(out) != nrow(current_data)) {
    stop("Each flow rule must return a logical vector matching the current data.", call. = FALSE)
  }
  out
}

.ukb_flow_rule_type <- function(rule) {
  if (inherits(rule, "formula")) {
    return("formula")
  }
  if (is.function(rule)) {
    return("function")
  }
  if (is.logical(rule)) {
    return("logical")
  }
  if (is.character(rule)) {
    return("complete_case")
  }
  "unknown"
}
