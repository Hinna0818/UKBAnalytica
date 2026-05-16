# `ukbsci-baseline` — function reference

## `create_baseline_table()`

```r
create_baseline_table(data, case_col,
                      factor_cols     = NULL,
                      continuous_cols = NULL,
                      test            = FALSE)
```

| Arg | Type | Meaning |
|-----|------|---------|
| `data` | data.frame / data.table | cohort |
| `case_col` | char | binary stratifier (0/1 or 2-level factor) |
| `factor_cols` | char vector | categorical variables |
| `continuous_cols` | char vector | continuous variables |
| `test` | logical | run group comparison tests |

**Returns:** `tableone::TableOne` S3 object.

**Wraps:** `tableone::CreateTableOne(strata = case_col, ...)`.

**Test logic (automatic via tableone):**

| Variable kind | Default test |
|---------------|--------------|
| Continuous, two groups | t-test |
| Continuous, ≥ 3 groups | one-way ANOVA |
| Continuous, declared `nonnormal` in `print()` | Wilcoxon / Kruskal-Wallis |
| Categorical | Chi-square (Fisher's exact if cell counts ≤ 5 by default in tableone) |

**Caveats:**

- Errors if `tableone` is not installed.
- Errors if `case_col` is missing from `data`.
- Coerces `case_col` to factor internally.
- Rows with `NA` in `case_col` are dropped by `tableone`.

---

## Recommended `print()` calling convention

```r
mat <- print(tab1,
             showAllLevels = TRUE,   # show every factor level
             quote         = FALSE,
             noSpaces      = TRUE,   # CSV-friendly
             smd           = FALSE,  # set TRUE for matched cohorts
             missing       = TRUE,   # missingness column
             nonnormal     = c(),    # Wilcoxon vars
             exact         = c(),    # Fisher's exact vars
             test          = TRUE)
write.csv(mat, "/mnt/project/<area>/04-results/01-baseline_table.csv")
```

`mat` is a character matrix; preserve row/column names for the manuscript
Table 1 layout.
