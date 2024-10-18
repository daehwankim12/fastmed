#' Perform Mediation Analysis
#'
#' This function performs mediation analysis using bootstrap resampling to estimate direct, indirect, and total effects.  It processes data in chunks to handle large datasets efficiently and writes results directly to a CSV file.
#'
#' @param data A `data.table` or `data.frame` containing the data.  If a `data.frame` is provided, it will be converted to a `data.table` internally.
#' @param columns A list with three named elements: `exposure`, `mediator`, and `outcome`. Each element should be a character vector containing the prefixes of the column names for the corresponding variables. For example, if your exposure variables are named "Exposure_1", "Exposure_2", etc., the `exposure` element should be `"Exposure"`.
#' @param nrep An integer specifying the number of bootstrap replicates to perform.  Higher values generally lead to more stable estimates but increase computation time. Default is 1000.
#' @param output_file A character string specifying the path to the output CSV file.  Results are written to this file in real-time.  The file will be overwritten if it already exists.
#' @param num_threads An integer specifying the number of threads to use for parallel processing.  Default is the number of available cores detected by `parallel::detectCores()`.
#' @param pert A character string specifying the method to use for calculating standard errors and confidence intervals.  Options are `"asymptotic"` and `"bootstrap"`.  The `"bootstrap"` method is more computationally intensive but may provide more accurate results for small sample sizes. If your data is non-parametric, the `"bootstrap"` method is recommended. Default is `"asymptotic"`.
#' @return None.  The results are written to the specified `output_file` in CSV format.
#'
#' @details This function estimates the following effects:
#' * **Indirect Effect (t1):**  The effect of the exposure on the outcome through the mediator when the exposure is set to its value at t=1 (typically meaning its presence or a higher level).
#' * **Indirect Effect (t0):** The effect of the exposure on the outcome through the mediator when the exposure is set to its value at t=0 (typically meaning its absence or a lower level).
#' * **Direct Effect (t1):** The effect of the exposure on the outcome not through the mediator when the exposure is set to t=1.
#' * **Direct Effect (t0):** The effect of the exposure on the outcome not through the mediator when the exposure is set to t=0.
#' * **Total Effect:**  The total effect of the exposure on the outcome, both direct and indirect.
#'
#' The output CSV file includes estimates, standard errors, 95% confidence intervals (LCB, UCB), and p-values for each effect and combination of variables.
#'
#' @examples
#' \dontrun{
#' # Example data (replace with your own data)
#' my_data <- data.table(
#'   Exposure_A = rnorm(1000),
#'   Exposure_B = rnorm(1000),
#'   Mediator_X = rnorm(1000),
#'   Mediator_Y = rnorm(1000),
#'   Outcome_1 = rnorm(1000),
#'   Outcome_2 = rnorm(1000)
#' )
#'
#' # Perform mediation analysis
#' mediation_analysis(
#'   data = my_data,
#'   columns = list(exposure = "Exposure", mediator = "Mediator", outcome = "Outcome"),
#'   nrep = 500, # Reduced for example speed
#'   output_file = "mediation_results.csv"
#' )
#' }
#' 
#' @import data.table Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel setThreadOptions
#' @useDynLib fastmed, .registration = TRUE
#' @export
mediation_analysis <- function(data,
                               columns,
                               nrep = 1000,
                               output_file,
                               num_threads = parallel::detectCores(),
                               pert = "asymptotic") {
  stopifnot(nrep > 0, 
            is.character(output_file), 
            is.list(columns) && length(columns) == 3,
            pert %in% c("asymptotic", "bootstrap"))
  
  # Ensure data is a data.table
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }
  
  data <- data[, lapply(.SD, function(col) {
    if (is.integer(col)) {
      as.numeric(col)
    } else {
      col
    }
  })]
  
  # Set number of threads for RcppParallel
  RcppParallel::setThreadOptions(numThreads = num_threads)
  
  # Function to find columns that start with given prefixes
  find_columns <- function(prefixes, all_columns) {
    unique(unlist(lapply(prefixes, function(prefix) {
      grep(paste0("^", prefix), all_columns, value = TRUE)
    })))
  }
  
  # Find actual column names based on prefixes
  all_columns <- names(data)
  exposure_cols <- find_columns(columns$exposure, all_columns)
  mediator_cols <- find_columns(columns$mediator, all_columns)
  outcome_cols <- find_columns(columns$outcome, all_columns)
  
  # Check if columns were found
  if (length(exposure_cols) == 0 ||
      length(mediator_cols) == 0 || length(outcome_cols) == 0) {
    stop("No columns found for one or more of exposure, mediator, or outcome.")
  }
  
  # Generate all combinations
  combinations <- expand.grid(
    exposure = exposure_cols,
    mediator = mediator_cols,
    outcome = outcome_cols,
    stringsAsFactors = FALSE
  )
  
  # Check uniqueness of column names
  unique_cols <- unique(c(exposure_cols, mediator_cols, outcome_cols))
  if (length(unique_cols) < length(c(exposure_cols, mediator_cols, outcome_cols))) {
    stop("Column names must be unique across exposure, mediator, and outcome variables.")
  }
  
  # Call C++ function
  mediation_analysis_cpp(
    as.matrix(data),
    colnames(data),
    combinations,
    nrep,
    output_file,
    pert
  )
  
  cat(
    "Mediation analysis completed. Results saved to",
    output_file,
    "\n"
  )
}
