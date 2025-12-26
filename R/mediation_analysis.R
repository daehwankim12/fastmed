#' Perform Mediation Analysis
#'
#' This function performs mediation analysis using multiple linear regressions
#' and either an asymptotic (normal) perturbation or bootstrap resampling to
#' estimate direct, indirect, and total effects. Results are written to a CSV
#' file.
#'
#' @param data A `data.table` or `data.frame` containing the data.  If a `data.frame` is provided, it will be converted to a `data.table` internally.
#' @param columns A list with three named elements: `exposure`, `mediator`, and `outcome`. Each element should be a character vector containing the prefixes of the column names for the corresponding variables. For example, if your exposure variables are named "Exposure_1", "Exposure_2", etc., the `exposure` element should be `"Exposure"`.
#' @param nrep An integer specifying the number of bootstrap replicates to perform.  Higher values generally lead to more stable estimates but increase computation time. Default is 1000.
#' @param output_file A character string specifying the path to the output CSV file.  Results are written to this file in real-time.  The file will be overwritten if it already exists.
#' @param num_threads An integer specifying the number of threads to use for parallel processing.  Default is the number of available cores detected by `parallel::detectCores()`.
#' @param pert A character string specifying the method to use for uncertainty
#'   estimation. Options are `"asymptotic"` and `"bootstrap"`. Default is
#'   `"asymptotic"`.
#' @param seed Optional non-negative integer. If provided, results are
#'   reproducible across thread counts within the same build/runtime
#'   environment.
#' @param chunk_size Maximum number of (exposure, mediator, outcome) combinations
#'   to process per call into the C++ backend. For large analyses, smaller
#'   values reduce peak memory usage. Default is 10000.
#' @return None.  The results are written to the specified `output_file` in CSV format.
#'
#' @details This function estimates the following effects:
#' * **ACME(0):** Average causal mediation (indirect) effect at baseline treatment.
#' * **ADE(0):** Average direct effect at baseline treatment.
#' * **Total Effect:**  The total effect of the exposure on the outcome, both direct and indirect.
#'
#' The output CSV file includes mean estimates, 95% percentile confidence
#' intervals, and p-values for each effect and combination of variables.
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
                               pert = "asymptotic",
                               seed = NULL,
                               chunk_size = 10000) {
  if (!is.numeric(nrep) || length(nrep) != 1 || is.na(nrep) || nrep <= 0) {
    stop("nrep must be a positive integer.")
  }
  nrep <- as.integer(nrep)

  if (!is.character(output_file) || length(output_file) != 1 || is.na(output_file) || !nzchar(output_file)) {
    stop("output_file must be a non-empty character scalar.")
  }

  if (!is.numeric(num_threads) || length(num_threads) != 1 || is.na(num_threads) || num_threads <= 0) {
    stop("num_threads must be a positive integer.")
  }
  num_threads <- as.integer(num_threads)

  if (!is.character(pert) || length(pert) != 1 || is.na(pert) || !pert %in% c("asymptotic", "bootstrap")) {
    stop("pert must be one of: 'asymptotic', 'bootstrap'.")
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || is.na(chunk_size) || chunk_size <= 0) {
    stop("chunk_size must be a positive integer.")
  }
  chunk_size <- as.integer(chunk_size)

  validate_data(data)
  validate_columns(columns)

  if (data.table::is.data.table(data)) {
    data <- data.table::copy(data)
  } else {
    data <- data.table::as.data.table(data)
  }

  data[] <- lapply(data, function(col) {
    if (is.integer(col)) {
      as.numeric(col)
    } else {
      col
    }
  })

  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "double"

  # Set number of threads for RcppParallel
  RcppParallel::setThreadOptions(numThreads = num_threads)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) || !is.finite(seed)) {
      stop("seed must be a single, non-negative integer.")
    }
    if (seed < 0 || seed > .Machine$integer.max) {
      stop("seed must be between 0 and .Machine$integer.max.")
    }
    if (seed != floor(seed)) {
      stop("seed must be an integer.")
    }
    base_seed <- as.integer(seed)
  } else {
    base_seed <- sample.int(.Machine$integer.max, 1)
  }

  all_columns <- names(data)
  exposure_cols <- find_matching_columns(columns$exposure, all_columns, "exposure")
  mediator_cols <- find_matching_columns(columns$mediator, all_columns, "mediator")
  outcome_cols <- find_matching_columns(columns$outcome, all_columns, "outcome")

  # Check uniqueness of column names
  unique_cols <- unique(c(exposure_cols, mediator_cols, outcome_cols))
  if (length(unique_cols) < length(c(exposure_cols, mediator_cols, outcome_cols))) {
    stop("Column names must be unique across exposure, mediator, and outcome variables.")
  }

  n_exposure <- length(exposure_cols)
  n_mediator <- length(mediator_cols)
  n_outcome <- length(outcome_cols)
  total_combinations <- as.double(n_exposure) * as.double(n_mediator) * as.double(n_outcome)

  if (!is.finite(total_combinations) || total_combinations < 1) {
    stop("No (exposure, mediator, outcome) combinations to process.")
  }
  if (total_combinations > .Machine$integer.max) {
    stop("Too many combinations for chunking in R (exceeds .Machine$integer.max). Reduce the number of columns.")
  }
  total_combinations <- as.integer(total_combinations)

  if (total_combinations <= chunk_size) {
    combinations <- expand.grid(
      exposure = exposure_cols,
      mediator = mediator_cols,
      outcome = outcome_cols,
      stringsAsFactors = FALSE
    )
    combinations$global_idx <- seq_len(nrow(combinations)) - 1L

    mediation_analysis_cpp(
      data_mat,
      colnames(data_mat),
      combinations,
      nrep,
      output_file,
      pert,
      base_seed,
      append = FALSE
    )
  } else {
    processed <- 0L
    first_chunk <- TRUE

    message(sprintf("Processing %d combinations in chunks of %d", total_combinations, chunk_size))

    e <- n_exposure
    m <- n_mediator
    em <- e * m

    while (processed < total_combinations) {
      n_chunk <- min(chunk_size, total_combinations - processed)
      idx <- processed + seq_len(n_chunk) - 1L

      exposure_idx <- (idx %% e) + 1L
      mediator_idx <- ((idx %/% e) %% m) + 1L
      outcome_idx <- (idx %/% em) + 1L

      combinations <- data.frame(
        exposure = exposure_cols[exposure_idx],
        mediator = mediator_cols[mediator_idx],
        outcome = outcome_cols[outcome_idx],
        global_idx = idx,
        stringsAsFactors = FALSE
      )

      mediation_analysis_cpp(
        data_mat,
        colnames(data_mat),
        combinations,
        nrep,
        output_file,
        pert,
        base_seed,
        append = !first_chunk
      )

      first_chunk <- FALSE
      processed <- processed + n_chunk

      if (interactive() && processed %% (10L * chunk_size) == 0L) {
        message(sprintf("Progress: %d/%d (%.1f%%)", processed, total_combinations, 100 * processed / total_combinations))
      }
    }
  }

  cat("Mediation analysis completed. Results saved to", output_file, "\n")
}

validate_data <- function(data) {
  if (!data.table::is.data.table(data) && !is.data.frame(data)) {
    stop("Data must be a data.table or data.frame.")
  }

  if (nrow(data) == 0) {
    stop("Data is empty.")
  }

  if (anyNA(data)) {
    stop("Data contains missing values. Remove or impute NAs before analysis.")
  }

  non_numeric_cols <- names(data)[!vapply(data, is.numeric, logical(1))]
  if (length(non_numeric_cols) > 0) {
    stop(paste0(
      "Non-numeric columns found: ",
      paste(non_numeric_cols, collapse = ", "),
      ". All columns must be numeric."
    ))
  }

  invisible(TRUE)
}

validate_columns <- function(columns) {
  if (!is.list(columns) || length(columns) != 3) {
    stop("columns must be a list with exactly 3 elements: exposure, mediator, outcome")
  }

  required <- c("exposure", "mediator", "outcome")
  missing <- setdiff(required, names(columns))
  if (length(missing) > 0) {
    stop(paste0("Missing required column types: ", paste(missing, collapse = ", ")))
  }

  for (type_name in required) {
    prefixes <- columns[[type_name]]
    if (!is.character(prefixes) || length(prefixes) < 1) {
      stop(paste0("columns$", type_name, " must be a non-empty character vector of prefixes."))
    }
  }

  invisible(TRUE)
}

find_matching_columns <- function(prefixes, all_columns, type_name) {
  matches <- unique(unlist(lapply(prefixes, function(prefix) {
    all_columns[startsWith(all_columns, prefix)]
  })))

  if (length(matches) == 0) {
    stop(paste0(
      "No columns found matching ", type_name, " prefixes: ",
      paste(prefixes, collapse = ", ")
    ))
  }

  matches
}
