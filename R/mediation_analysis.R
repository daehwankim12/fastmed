#' Perform Mediation Analysis
#'
#' This function performs mediation analysis using bootstrap resampling to estimate direct, indirect, and total effects.  It processes data in chunks to handle large datasets efficiently and writes results directly to a CSV file.
#'
#' @param data A `data.table` or `data.frame` containing the data.  If a `data.frame` is provided, it will be converted to a `data.table` internally.
#' @param columns A list with three named elements: `exposure`, `mediator`, and `outcome`. Each element should be a character vector containing the prefixes of the column names for the corresponding variables. For example, if your exposure variables are named "Exposure_1", "Exposure_2", etc., the `exposure` element should be `"Exposure"`.
#' @param nrep An integer specifying the number of bootstrap replicates to perform.  Higher values generally lead to more stable estimates but increase computation time. Default is 1000.
#' @param output_file A character string specifying the path to the output CSV file.  Results are written to this file in real-time.  The file will be overwritten if it already exists.
#' @param num_threads An integer specifying the number of threads to use for parallel processing.  Default is the number of available cores detected by `parallel::detectCores()`.
#' @return None.  The results are written to the specified `output_file` in CSV format.
#'
#' @details
#' ... (rest of the documentation remains the same)
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
                               num_threads = parallel::detectCores()) {
  # Ensure data is a data.table
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }
  
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
    output_file
  )
  
  cat(
    "Mediation analysis completed. Results saved to",
    output_file,
    "\n"
  )
}
