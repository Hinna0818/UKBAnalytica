#' Create a baseline table comparing cases and controls under different conditions.
#' @references https://github.com/kaz-yos/tableone
#' @param data a data.table containing the baseline characteristics and case/control status.
#' @param case_col the name of the column indicating case/control status (1 for cases, 0 for controls).
#' @param factor_cols a vector of column names that are factors (categorical variables)
#' @param continuous_cols a vector of column names that are continuous variables.
#' @param test whether to perform statistical tests comparing cases and controls for each variable (default: FALSE).
#' @return a list containing table one information.
#' @import data.table
#' @importFrom tableone CreateTableOne
#' @export
create_baseline_table <- function(data, case_col, factor_cols = NULL, continuous_cols = NULL, test = FALSE) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }
  
  # Ensure case_col is a factor
  data[[case_col]] <- as.factor(data[[case_col]])
  
  # Ensure factor columns are factors, and continuous columns are numeric
  if (!is.null(factor_cols)) {
    for (col in factor_cols) {
      data[[col]] <- as.factor(data[[col]])
    }
  }
  if (!is.null(continuous_cols)) {
    for (col in continuous_cols) {
      data[[col]] <- as.numeric(data[[col]])  
    }
  }
  
  # Create the tableone object
  table_one <- CreateTableOne(
    vars = c(factor_cols, continuous_cols),
    strata = case_col,
    data = data,
    factorVars = factor_cols,
    test = test
  )
  
  return(table_one)
}
