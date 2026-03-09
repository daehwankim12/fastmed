#' Perform Mediation Analysis
#'
#' This function performs mediation analysis using generalized linear models
#' (Gaussian/Binomial/Poisson) and either an asymptotic (quasi-Bayesian)
#' perturbation or bootstrap resampling to estimate direct, indirect, and total
#' effects. Results are written to a CSV file.
#'
#' @param data A `data.table` or `data.frame` containing the data.  If a `data.frame` is provided, it will be converted to a `data.table` internally.
#' @param columns A list with three named elements: `exposure`, `mediator`, and `outcome`. Each element should be a character vector containing the prefixes of the column names for the corresponding variables. For example, if your exposure variables are named "Exposure_1", "Exposure_2", etc., the `exposure` element should be `"Exposure"`.
#' @param weights Optional numeric vector of prior weights (length `nrow(data)`). Non-negative, finite, and must sum to a positive value. Default is `NULL` (unweighted).
#' @param treat.value Numeric scalar giving the treatment value for `T=1` potential outcomes. Default is 1.
#' @param control.value Numeric scalar giving the control value for `T=0` potential outcomes. Default is 0.
#' @param nrep An integer specifying the number of simulation draws (asymptotic)
#'   or bootstrap replicates (bootstrap). Higher values generally lead to more
#'   stable estimates but increase computation time. Default is 1000.
#' @param output_file A character string specifying the path to the output CSV
#'   file. Results are written to this file in real-time.
#' @param overwrite Logical; if `FALSE` and `output_file` already exists, throw
#'   an error instead of overwriting. Default is `TRUE` (preserves existing
#'   behavior).
#' @param excel_safe_csv Logical; if `TRUE`, prefix potentially dangerous
#'   spreadsheet formula strings (e.g. values starting with `=`, `+`, `-`, `@`)
#'   with a leading `'` in the output CSV. Default is `FALSE` (preserves
#'   existing output).
#' @param num_threads An integer specifying the number of threads to use for parallel processing.  Default is the number of available cores detected by `parallel::detectCores()`.
#' @param pert A character string specifying the method to use for uncertainty
#'   estimation. Options are `"asymptotic"` and `"bootstrap"`. Default is
#'   `"asymptotic"`.
#' @param seed Optional non-negative integer. If provided, results are
#'   reproducible across thread counts within the same build/runtime
#'   environment. When `rng_mode = "mediate"` and `seed = NULL`, fastmed uses
#'   the current R RNG stream (it does not call `set.seed()` internally).
#' @param match_mediation Deprecated alias for `rng_mode`. If supplied, it
#'   overrides `rng_mode`.
#' @param rng_mode Random number mode. `"fast"` uses fast deterministic
#'   per-combination seeding; `"mediate"` aligns random number usage more
#'   closely with `mediation::mediate()`. Default is `"fast"`.
#' @param chunk_size Maximum number of (exposure, mediator, outcome) combinations
#'   to buffer per chunk inside the C++ backend. Smaller values reduce peak
#'   memory usage for very large analyses. Default is 10000.
#' @param grain_size Number of (exposure, mediator, outcome) combinations to
#'   process per parallel task inside the C++ backend. Larger values reduce
#'   scheduling overhead when individual combinations are cheap. Default is 100.
#' @param loop_order Controls the order in which (exposure, mediator, outcome)
#'   combinations are processed/written. Provide either a length-3 character
#'   vector (outer-to-inner) containing a permutation of `c("exposure",
#'   "mediator", "outcome")`, or a shorthand string like `"EMO"` (the default).
#' @param mediator.family Model family for the mediator regression. One of
#'   `"auto"`, `"gaussian"`, `"binomial"`, `"poisson"`. Default is `"auto"`.
#' @param outcome.family Model family for the outcome regression. One of
#'   `"auto"`, `"gaussian"`, `"binomial"`, `"poisson"`. Default is `"auto"`.
#' @param replace.outcome Logical; if TRUE, replace the simulated `Y(0,M(0))`
#'   for observed controls and `Y(1,M(1))` for observed treated units with the
#'   observed outcome. Default is FALSE.
#' @param output.format Output CSV schema. `"mediate"` writes mediate-style
#'   effect columns (`d0`, `d1`, `z0`, `z1`, `tau`). `"legacy"` writes the
#'   previous schema (`ACME_*`, `ADE_*`, `Total_Effect_*`). Default is
#'   `"mediate"`.
#' @param failure_mode How to handle failing combinations. `"na_row"` writes an
#'   `NA` row (default). `"error"` stops immediately with an error.
#' @param include_failure_reason Logical; if `TRUE`, add a `failure_reason`
#'   column to the output CSV (`NA` for successful rows, reason text for
#'   failed rows). Default is `FALSE`.
#' @return None.  The results are written to the specified `output_file` in CSV format.
#'
#' @details This function estimates the following effects:
#' * **d0:** ACME(control) average causal mediation effect.
#' * **d1:** ACME(treated) average causal mediation effect.
#' * **z0:** ADE(control) average direct effect.
#' * **z1:** ADE(treated) average direct effect.
#' * **tau:** Total effect.
#'
#' The output CSV file includes mean estimates, 95% percentile confidence
#' intervals, and p-values for each effect and combination of variables.
#'
#' Family auto-detection is based on the response values:
#' * All values are 0/1 -> Binomial
#' * All values are non-negative integers -> Poisson
#' * Otherwise -> Gaussian
#'
#' Missing values are handled per combination: for each (exposure, mediator,
#' outcome) triplet, rows with `NA`/`NaN` in any of those three variables are
#' dropped prior to fitting and simulation/bootstrapping (complete-case analysis
#' for that combination). If a combination has too few complete cases to fit the
#' models, the corresponding output row is `NA`.
#'
#' Infinite values (`Inf`, `-Inf`) are rejected at input validation time.
#'
#' Setting `replace.outcome = TRUE` replaces some simulated outcomes with
#' observed outcomes and may reduce agreement with `mediation::mediate()`.
#'
#' If `failure_mode = "na_row"`, failing combinations are kept in output as
#' `NA` rows. If `failure_mode = "error"`, analysis stops on first failure.
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
                               match_mediation = NULL,
                               rng_mode = c("fast", "mediate"),
                               chunk_size = 10000,
                               grain_size = 100,
                               loop_order = c("exposure", "mediator", "outcome"),
                               mediator.family = "auto",
                               outcome.family = "auto",
                               replace.outcome = FALSE,
                               output.format = c("mediate", "legacy"),
                               weights = NULL,
                               treat.value = 1,
                               control.value = 0,
                               overwrite = TRUE,
                               excel_safe_csv = FALSE,
                               failure_mode = c("na_row", "error"),
                               include_failure_reason = FALSE) {
  if (!is.numeric(nrep) || length(nrep) != 1 || is.na(nrep) || nrep <= 0) {
    stop("nrep must be a positive integer.")
  }
  nrep <- as.integer(nrep)

  if (!is.character(output_file) || length(output_file) != 1 || is.na(output_file) || !nzchar(output_file)) {
    stop("output_file must be a non-empty character scalar.")
  }

  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  if (!is.logical(excel_safe_csv) || length(excel_safe_csv) != 1L || is.na(excel_safe_csv)) {
    stop("excel_safe_csv must be TRUE or FALSE.")
  }

  if (!is.numeric(num_threads) || length(num_threads) != 1 || is.na(num_threads) || num_threads <= 0) {
    stop("num_threads must be a positive integer.")
  }
  num_threads <- as.integer(num_threads)

  if (!is.character(pert) || length(pert) != 1 || is.na(pert) || !pert %in% c("asymptotic", "bootstrap")) {
    stop("pert must be one of: 'asymptotic', 'bootstrap'.")
  }

  rng_mode <- match.arg(rng_mode)
  if (!is.null(match_mediation)) {
    if (!is.logical(match_mediation) || length(match_mediation) != 1L || is.na(match_mediation)) {
      stop("match_mediation must be TRUE, FALSE, or NULL.")
    }
    warning(
      "match_mediation is deprecated; use rng_mode = 'mediate' or 'fast' instead.",
      call. = FALSE
    )
    rng_mode <- if (isTRUE(match_mediation)) "mediate" else "fast"
  }
  match_mediation_cpp <- identical(rng_mode, "mediate")

  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || is.na(chunk_size) || chunk_size <= 0) {
    stop("chunk_size must be a positive integer.")
  }
  chunk_size <- as.integer(chunk_size)

  if (!is.numeric(grain_size) || length(grain_size) != 1 || is.na(grain_size) || grain_size <= 0) {
    stop("grain_size must be a positive integer.")
  }
  grain_size <- as.integer(grain_size)

  normalize_loop_order <- function(loop_order) {
    if (!is.character(loop_order) || anyNA(loop_order)) {
      stop("loop_order must be a character vector with no NA values.")
    }

    if (length(loop_order) == 1L) {
      if (!nzchar(loop_order)) stop("loop_order must be non-empty.")
      key <- tolower(loop_order)
      key <- gsub("[[:space:],./-]+", "_", key)

      if (nchar(key) == 3L && grepl("^[emo]{3}$", key)) {
        tokens <- strsplit(key, "", fixed = TRUE)[[1]]
      } else {
        tokens <- strsplit(key, "_", fixed = TRUE)[[1]]
        tokens <- tokens[nzchar(tokens)]
      }
    } else if (length(loop_order) == 3L) {
      tokens <- tolower(loop_order)
    } else {
      stop("loop_order must be a length-1 string (e.g. 'EMO') or length-3 character vector.")
    }

    map_token <- function(tok) {
      if (tok %in% c("e", "exp", "exposure")) return("exposure")
      if (tok %in% c("m", "med", "mediator")) return("mediator")
      if (tok %in% c("o", "out", "outcome")) return("outcome")
      stop("Invalid loop_order token: ", tok)
    }

    tokens <- vapply(tokens, map_token, FUN.VALUE = character(1))

    if (length(unique(tokens)) != 3L) {
      stop("loop_order must contain each of: exposure, mediator, outcome exactly once.")
    }

    paste(tokens, collapse = "_")
  }

  loop_order_cpp <- normalize_loop_order(loop_order)

  valid_families <- c("gaussian", "binomial", "poisson", "auto")
  if (!is.character(mediator.family) || length(mediator.family) != 1L || is.na(mediator.family)) {
    stop("mediator.family must be a single character string.")
  }
  if (!is.character(outcome.family) || length(outcome.family) != 1L || is.na(outcome.family)) {
    stop("outcome.family must be a single character string.")
  }
  mediator.family <- tolower(mediator.family)
  outcome.family <- tolower(outcome.family)
  if (!mediator.family %in% valid_families) {
    stop("mediator.family must be one of: ", paste(valid_families, collapse = ", "))
  }
  if (!outcome.family %in% valid_families) {
    stop("outcome.family must be one of: ", paste(valid_families, collapse = ", "))
  }
  if (!is.logical(replace.outcome) || length(replace.outcome) != 1L || is.na(replace.outcome)) {
    stop("replace.outcome must be TRUE or FALSE.")
  }

  output.format <- match.arg(output.format)
  failure_mode <- match.arg(failure_mode)
  if (!is.logical(include_failure_reason) || length(include_failure_reason) != 1L || is.na(include_failure_reason)) {
    stop("include_failure_reason must be TRUE or FALSE.")
  }

  validate_data(data)
  validate_columns(columns)

  if (!is.null(weights)) {
    if (!is.numeric(weights) || is.matrix(weights) || length(weights) != nrow(data)) {
      stop("weights must be a numeric vector of length nrow(data), or NULL.")
    }
    weights <- as.numeric(weights)
    if (any(!is.finite(weights))) {
      stop("weights must be finite (no NA/Inf).")
    }
    if (any(weights < 0)) {
      stop("weights must be non-negative (zeros allowed).")
    }
    wsum <- sum(weights)
    if (!is.finite(wsum) || wsum <= 0) {
      stop("weights must sum to a positive finite value.")
    }
  }

  if (!is.numeric(treat.value) || length(treat.value) != 1L || is.na(treat.value) || !is.finite(treat.value)) {
    stop("treat.value must be a single, finite numeric value.")
  }
  if (!is.numeric(control.value) || length(control.value) != 1L || is.na(control.value) || !is.finite(control.value)) {
    stop("control.value must be a single, finite numeric value.")
  }
  treat.value <- as.numeric(treat.value)
  control.value <- as.numeric(control.value)
  if (treat.value == control.value) {
    stop("treat.value and control.value must be different.")
  }

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
    if (match_mediation_cpp) {
      set.seed(base_seed)
    }
  } else if (!match_mediation_cpp) {
    base_seed <- sample.int(.Machine$integer.max, 1)
  } else {
    # In mediate-parity mode with no explicit seed, consume the current R RNG
    # stream exactly as-is (no internal reseeding).  base_seed is a sentinel
    # (ignored by the C++ match_mediation path, which uses precomputed R-side
    # RNG draws instead of derive_seed()).
    base_seed <- 0L
  }

  overwrite <- isTRUE(overwrite)
  excel_safe_csv <- isTRUE(excel_safe_csv)

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

  validate_response_family <- function(y, family_name) {
    eps <- 1e-8
    y <- y[!is.na(y)]

    if (family_name == "binomial") {
      if (any(!is.finite(y))) stop("Binomial y contains non-finite")
      ok <- all(abs(y) <= eps | abs(y - 1) <= eps)
      if (!ok) stop("Binomial requires y in {0,1} (found non-binary).")
      return(invisible(TRUE))
    }

    if (family_name == "poisson") {
      if (any(!is.finite(y))) stop("Poisson y contains non-finite")
      if (any(y < -eps)) stop("Poisson requires y >= 0.")
      ok <- all(abs(y - round(y)) <= eps)
      if (!ok) stop("Poisson requires integer y.")
      return(invisible(TRUE))
    }

    invisible(TRUE)
  }

  if (mediator.family %in% c("binomial", "poisson")) {
    for (col in mediator_cols) {
      validate_response_family(data[[col]], mediator.family)
    }
  }
  if (outcome.family %in% c("binomial", "poisson")) {
    for (col in outcome_cols) {
      validate_response_family(data[[col]], outcome.family)
    }
  }

  exp_idx <- match(exposure_cols, colnames(data_mat)) - 1L
  med_idx <- match(mediator_cols, colnames(data_mat)) - 1L
  out_idx <- match(outcome_cols, colnames(data_mat)) - 1L

  mediation_analysis_cpp(
    data = data_mat,
    column_names = colnames(data_mat),
    exposure_col_idx = exp_idx,
    mediator_col_idx = med_idx,
    outcome_col_idx = out_idx,
    weights = weights,
    nrep = nrep,
    output_file = output_file,
    pert = pert,
    base_seed = base_seed,
    mediator_family = mediator.family,
    outcome_family = outcome.family,
    replace_outcome = isTRUE(replace.outcome),
    output_format = output.format,
    treat_value = treat.value,
    control_value = control.value,
    chunk_size = chunk_size,
    grain_size = grain_size,
    overwrite = overwrite,
    excel_safe_csv = excel_safe_csv,
    match_mediation = match_mediation_cpp,
    loop_order = loop_order_cpp,
    fail_fast = identical(failure_mode, "error"),
    include_failure_reason = isTRUE(include_failure_reason)
  )

  cat("Mediation analysis completed. Results saved to", output_file, "\n")
}

validate_data <- function(data) {
  if (!data.table::is.data.table(data) && !is.data.frame(data)) {
    stop("Data must be a data.table or data.frame.")
  }

  if (nrow(data) == 0) {
    stop("Data is empty.")
  }

  non_numeric_cols <- names(data)[!vapply(data, is.numeric, logical(1))]
  if (length(non_numeric_cols) > 0) {
    stop(paste0(
      "Non-numeric columns found: ",
      paste(non_numeric_cols, collapse = ", "),
      ". All columns must be numeric."
    ))
  }

  non_finite_cols <- names(data)[vapply(data, function(col) {
    any(is.infinite(col))
  }, logical(1))]
  if (length(non_finite_cols) > 0) {
    stop(paste0(
      "Infinite values found in columns: ",
      paste(non_finite_cols, collapse = ", "),
      ". Replace Inf/-Inf with NA or finite values."
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
    if (anyNA(prefixes) || any(!nzchar(prefixes))) {
      stop(paste0("columns$", type_name, " must not contain NA or empty strings."))
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
