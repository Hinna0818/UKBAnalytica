# `ukbsci-baseline` — examples

```r
library(UKBAnalytica); library(tableone); library(data.table)
```

## A. Vanilla — stratified by disease history

```r
factor_cols <- c("sex", "ethnic_background", "smoking_status",
                 "Hypertension_history", "Diabetes_history")
cont_cols   <- c("age", "bmi", "sbp", "dbp", "ldl_cholesterol")

tab1 <- create_baseline_table(
  data            = cohort,
  case_col        = "AA_history",
  factor_cols     = factor_cols,
  continuous_cols = cont_cols,
  test            = TRUE
)

mat <- print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
write.csv(mat, "/mnt/project/<area>/04-results/01-baseline_table.csv")
```

## B. Skewed biomarkers → Wilcoxon + median[IQR]

```r
tab1 <- create_baseline_table(
  data = cohort, case_col = "AA_history",
  factor_cols = factor_cols,
  continuous_cols = c(cont_cols, "crp", "alt", "trigly"),
  test = TRUE
)

mat <- print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE,
             nonnormal = c("crp", "alt", "trigly"))
write.csv(mat, "/mnt/project/<area>/04-results/01-baseline_table_nonnormal.csv")
```

## C. After PSM — SMD instead of p-values

```r
# `matched` is the data.table from match_propensity()
tab1 <- create_baseline_table(
  data = matched, case_col = "treated",
  factor_cols = factor_cols, continuous_cols = cont_cols,
  test = FALSE
)

mat <- print(tab1, smd = TRUE, quote = FALSE, noSpaces = TRUE)
write.csv(mat, "/mnt/project/<area>/04-results/01-baseline_table_matched.csv")
```
