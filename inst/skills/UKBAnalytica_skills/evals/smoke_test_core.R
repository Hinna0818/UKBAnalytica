#!/usr/bin/env Rscript

# Minimal executable checks for the core UKBAnalytica skills.
# Run from the repository root:
#   Rscript inst/skills/UKBAnalytica_skills/evals/smoke_test_core.R

args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(args_file)) {
  script_path <- normalizePath(sub("^--file=", "", args_file[[1]]), mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
  setwd(repo_root)
}

quiet_load <- function() {
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(UKBAnalytica)
  }
}

quiet_load()
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

pass <- function(name) message("[PASS] ", name)
set.seed(20260516)

n <- 120L
toy <- data.table(
  eid = seq_len(n),
  p53_i0 = as.Date("2010-01-01") + sample(0:30, n, replace = TRUE),
  age = round(rnorm(n, 58, 7), 1),
  sex = factor(sample(c("Female", "Male"), n, replace = TRUE)),
  bmi = round(rnorm(n, 27, 4), 1),
  smoking = factor(sample(c("Never", "Former", "Current"), n, replace = TRUE)),
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n)
)
toy[, outcome := rbinom(.N, 1, plogis(-0.8 + 0.8 * x1 - 0.5 * x2))]
toy[, `:=`(p41270 = NA_character_, p41280_a0 = as.Date(NA))]
toy[1:3, `:=`(p41270 = "['J44']", p41280_a0 = p53_i0 - 365)]
toy[4:15, `:=`(p41270 = "['J44']", p41280_a0 = p53_i0 + 365)]

diseases <- get_predefined_diseases()
stopifnot(length(diseases) >= 70, "COPD" %in% names(diseases))
pass("predefined disease catalog")

cohort <- build_survival_dataset(
  dt = copy(toy),
  disease_definitions = diseases["COPD"],
  prevalent_sources = "ICD10",
  outcome_sources = "ICD10",
  censor_date = as.Date("2020-01-01"),
  baseline_col = "p53_i0",
  primary_disease = "COPD",
  show_flow = FALSE
)
stopifnot(all(c("COPD_history", "outcome_status", "outcome_surv_time") %in% names(cohort)))
stopifnot(sum(cohort$outcome_status, na.rm = TRUE) >= 10)
pass("cohort survival dataset")

tab <- create_baseline_table(
  data = toy,
  case_col = "outcome",
  factor_cols = c("sex", "smoking"),
  continuous_cols = c("age", "bmi"),
  test = TRUE
)
stopifnot(inherits(tab, "TableOne"))
pass("baseline table")

reg_dt <- copy(toy)
reg_dt[, `:=`(time = rexp(.N, rate = 0.1) + 0.1, status = outcome)]
cox_res <- runmulti_cox(
  data = reg_dt,
  main_var = c("age", "bmi"),
  covariates = "sex",
  endpoint = c("time", "status")
)
stopifnot(nrow(cox_res) == 2, all(c("HR", "pvalue") %in% names(cox_res)))
pass("multi-Cox regression")

split <- ukb_ml_split_data(
  df = toy[, .(outcome, x1, x2, x3)],
  outcome = "outcome",
  outcome_type = "binary",
  split = "train_test",
  train_ratio = 0.7,
  stratify_by = "outcome",
  seed = 20260516,
  verbose = FALSE
)
flow <- ukb_ml_flow(
  formula = outcome ~ x1 + x2 + x3,
  split = split,
  model = "logistic",
  outcome_type = "binary",
  tune = FALSE,
  threshold_method = "youden",
  seed = 20260516,
  verbose = FALSE
)
stopifnot(inherits(flow, "ukb_ml_flow"), nrow(flow$metrics) > 0, nrow(flow$roc) > 0)
stopifnot(!"eid" %in% names(flow$metrics), !"eid" %in% names(flow$roc))
pass("classification ML flow")

source("inst/skills/UKBAnalytica_skills/ukbsci-plot/references/theme-and-palettes.R")
stopifnot(all(c("ukbsci_clinical", "ukbsci_diverging", "ukbsci_sequential") %in% ls()))
p <- ggplot(toy, aes(x = age, y = bmi, colour = sex)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = unname(ukbsci_clinical[c("baseline", "exposure")])) +
  ukbsci_theme()
stopifnot(inherits(p, "ggplot"))
pass("plot theme and palettes")

message("All core UKBAnalytica skill smoke tests passed.")
