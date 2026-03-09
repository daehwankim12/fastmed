#' Serial mediation path analysis (lavaan-like)
#'
#' Fits three nested observed-variable path models for each \eqn{(X, M_1, \ldots, M_m, Y)}
#' combination, where \eqn{m} (the number of ordered mediator stages) can be any
#' positive integer.
#'
#' - **Partial:** each \eqn{M_j} is regressed on \eqn{X} and all earlier mediators,
#'   and \eqn{Y} is regressed on \eqn{X} and all mediators.
#' - **Semifull:** same, but without the direct path \eqn{Y ~ X}
#' - **Full:** adjacent serial chain only:
#'   \eqn{M_1 ~ X}, \eqn{M_j ~ M_{j-1}} for \eqn{j=2..m}, \eqn{Y ~ M_m}
#'
#' The function computes standard SEM-style fit measures (chi-square, df, p-value,
#' CFI, TLI, RMSEA, AIC) from the model-implied covariance matrix, performs
#' chi-square difference tests between nested models, selects a model using the
#' same decision logic as the common `lavaan::anova()` workflow, and then reports
#' path coefficients and defined indirect effects.
#'
#' Standard errors are estimated via a fast parametric bootstrap on the covariance
#' matrix (Wishart resampling) when `nrep > 0`.
#'
#' @param data A `data.table` or `data.frame` containing the data. All referenced
#'   columns must be numeric.
#' @param columns Column specification. Supported formats:
#'   - **New (stage list):** `list(x = ..., mediators = list(m1 = ..., m2 = ..., ...), y = ...)`.
#'     If `m` is omitted, `m` is inferred from `length(columns$mediators)`.
#'   - **Convenience (shared pool):** `list(x = ..., mediators = <character>, y = ...)`.
#'     The mediator pool is replicated across `m` stages (`m` defaults to 2).
#'   - **Legacy (deprecated):** `list(x = ..., m1 = ..., m2 = ..., y = ...)` (m=2),
#'     or `list(x = ..., m1 = ..., m2 = ..., m3 = ..., y = ...)` (m=3), etc.
#' @param output_file Deprecated alias for `output_dir` + `output_prefix`. If
#'   provided, results are written to `dirname(output_file)` with
#'   `output_prefix = tools::file_path_sans_ext(basename(output_file))`.
#' @param m Number of mediator stages when mediator stages are not explicitly
#'   listed in `columns` (default 2).
#' @param output_dir Output directory (default `"."`).
#' @param output_prefix Base name for output files (default `"serial_path"`).
#' @param nrep Integer. Number of parametric bootstrap replicates used to
#'   estimate standard errors. Use `0` to skip SE/z/p computation (faster).
#' @param alpha Significance level for nested model selection (default 0.05).
#' @param num_threads Number of threads for parallel processing.
#' @param seed Optional non-negative integer seed for reproducibility across
#'   thread counts (within the same build/runtime environment).
#' @param match How to interpret the strings in `columns`:
#'   - `"exact"`: treat as exact column names.
#'   - `"prefix"`: treat as prefixes (like `mediation_analysis()`).
#'   - `"auto"`: exact match when possible, otherwise prefix match.
#' @param combination_start 0-based global combination index start (inclusive).
#' @param combination_end 0-based global combination index end (exclusive).
#'   Use `NULL` for the full range.
#' @param shard_id Shard index in `[0, shard_count)` (default 0).
#' @param shard_count Number of shards (default 1).
#' @param chunk_size Maximum number of combinations per chunk inside the C++
#'   backend (controls peak memory).
#' @param grain_size Number of combinations per parallel task inside the C++
#'   backend.
#' @param warn_combinations Warn when combinations in this call exceed this.
#' @param max_combinations Stop when combinations in this call exceed this
#'   (default `Inf` to disable).
#' @param overwrite If `FALSE` and `output_file` already exists, error instead of
#'   overwriting.
#' @param excel_safe_csv If `TRUE`, prefix potentially dangerous spreadsheet
#'   formula strings (values starting with `=`, `+`, `-`, `@`) with a leading `'`
#'   in text fields.
#' @param include_failure_reason Logical; if `TRUE`, add `failure_reason` to fit
#'   output (`NA` for successful combinations, reason text for failed rows).
#' @param output Which output files to write. One or more of
#'   `c("fit", "params", "effects")`.
#'
#' @return None. Results are written to CSV files in `output_dir`.
#' @export
serial_path_analysis <- function(data,
                                 columns,
                                 output_file = NULL,
                                 m = 2,
                                 output_dir = ".",
                                 output_prefix = "serial_path",
                                 nrep = 200,
                                 alpha = 0.05,
                                 num_threads = parallel::detectCores(),
                                 seed = NULL,
                                 match = c("auto", "exact", "prefix"),
                                 combination_start = 0,
                                 combination_end = NULL,
                                 shard_id = 0,
                                 shard_count = 1,
                                 chunk_size = 10000,
                                 grain_size = 100,
                                 warn_combinations = 1e6,
                                 max_combinations = Inf,
                                 overwrite = TRUE,
                                 excel_safe_csv = FALSE,
                                 include_failure_reason = FALSE,
                                 output = c("fit", "params", "effects")) {
  validate_data(data)

  match <- match.arg(match)

  if (!is.numeric(nrep) || length(nrep) != 1 || is.na(nrep) || nrep < 0) {
    stop("nrep must be a single non-negative integer.")
  }
  nrep <- as.integer(nrep)

  if (!is.numeric(num_threads) || length(num_threads) != 1 || is.na(num_threads) || num_threads <= 0) {
    stop("num_threads must be a positive integer.")
  }
  num_threads <- as.integer(num_threads)
  RcppParallel::setThreadOptions(numThreads = num_threads)

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0,1).")
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || is.na(chunk_size) || chunk_size <= 0) {
    stop("chunk_size must be a positive integer.")
  }
  chunk_size <- as.integer(chunk_size)

  if (!is.numeric(grain_size) || length(grain_size) != 1 || is.na(grain_size) || grain_size <= 0) {
    stop("grain_size must be a positive integer.")
  }
  grain_size <- as.integer(grain_size)

  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  overwrite <- isTRUE(overwrite)

  if (!is.logical(excel_safe_csv) || length(excel_safe_csv) != 1L || is.na(excel_safe_csv)) {
    stop("excel_safe_csv must be TRUE or FALSE.")
  }
  excel_safe_csv <- isTRUE(excel_safe_csv)
  if (!is.logical(include_failure_reason) || length(include_failure_reason) != 1L || is.na(include_failure_reason)) {
    stop("include_failure_reason must be TRUE or FALSE.")
  }
  include_failure_reason <- isTRUE(include_failure_reason)

  resolve_columns <- function(spec, all_columns, type_name) {
    if (match == "exact") {
      bad <- setdiff(spec, all_columns)
      if (length(bad) > 0) {
        stop(paste0("Unknown ", type_name, " columns: ", paste(bad, collapse = ", ")))
      }
      return(unique(spec))
    }

    if (match == "prefix") {
      matches <- unique(unlist(lapply(spec, function(prefix) {
        all_columns[startsWith(all_columns, prefix)]
      })))
      if (length(matches) == 0) {
        stop(paste0("No columns found matching ", type_name, " prefixes: ", paste(spec, collapse = ", ")))
      }
      return(matches)
    }

    # auto
    out <- character()
    for (token in spec) {
      if (token %in% all_columns) {
        out <- c(out, token)
      } else {
        out <- c(out, all_columns[startsWith(all_columns, token)])
      }
    }
    out <- unique(out)
    if (length(out) == 0) {
      stop(paste0("No columns found for ", type_name, " (match='auto'): ", paste(spec, collapse = ", ")))
    }
    out
  }

  normalize_columns <- function(columns, all_columns, match, m = 2, m_is_missing = FALSE) {
    if (!is.list(columns) || is.null(names(columns))) {
      stop("columns must be a named list.")
    }

    if (!is.numeric(m) || length(m) != 1 || is.na(m) || m < 1 || m != floor(m)) {
      stop("m must be a positive integer.")
    }
    m <- as.integer(m)

    # Legacy format: columns$m1..columns$mK without columns$mediators
    if (!is.null(columns$m1) && is.null(columns$mediators)) {
      mediators <- list()
      i <- 1L
      while (!is.null(columns[[paste0("m", i)]])) {
        mediators[[paste0("m", i)]] <- columns[[paste0("m", i)]]
        i <- i + 1L
      }
      m_inferred <- length(mediators)
      if (m_inferred < 1L) stop("At least one mediator stage required.")
      if (!m_is_missing && m != m_inferred) {
        stop(sprintf(
          "m=%d does not match legacy columns (inferred m=%d from m1..mK).",
          m, m_inferred
        ))
      }
      m <- m_inferred
      columns <- list(x = columns$x, mediators = mediators, y = columns$y)
      warning("Legacy columns format (x/m1..mK/y) is deprecated; use columns$mediators instead.", call. = FALSE)
    }

    # Convenience: shared mediator pool -> expand to m stages (default m=2)
    if (!is.null(columns$mediators) && is.character(columns$mediators)) {
      mediators <- replicate(m, columns$mediators, simplify = FALSE)
      names(mediators) <- paste0("m", seq_len(m))
      columns$mediators <- mediators
    }

    # Stage list: infer m from list length if m omitted; error if conflicts when m explicit
    if (!is.null(columns$mediators) && is.list(columns$mediators)) {
      m_inferred <- length(columns$mediators)
      if (m_inferred < 1L) stop("At least one mediator stage required.")
      if (!m_is_missing && m != m_inferred) {
        stop(sprintf("m=%d does not match columns$mediators length=%d.", m, m_inferred))
      }
      m <- m_inferred
    }

    if (is.null(columns$x) || is.null(columns$mediators) || is.null(columns$y)) {
      stop("columns must contain x, mediators, and y.")
    }

    if (!is.character(columns$x) || length(columns$x) < 1 || anyNA(columns$x)) {
      stop("columns$x must be a non-empty character vector with no NA values.")
    }
    if (!is.character(columns$y) || length(columns$y) < 1 || anyNA(columns$y)) {
      stop("columns$y must be a non-empty character vector with no NA values.")
    }
    if (!is.list(columns$mediators) || length(columns$mediators) < 1) {
      stop("columns$mediators must be a non-empty list of mediator stages.")
    }
    for (i in seq_along(columns$mediators)) {
      v <- columns$mediators[[i]]
      if (!is.character(v) || length(v) < 1 || anyNA(v)) {
        stop(sprintf("columns$mediators[[%d]] must be a non-empty character vector with no NA values.", i))
      }
    }

    list(
      x = resolve_columns(columns$x, all_columns, "x"),
      mediators = lapply(seq_along(columns$mediators), function(i) {
        resolve_columns(columns$mediators[[i]], all_columns, paste0("m", i))
      }),
      y = resolve_columns(columns$y, all_columns, "y"),
      m = m
    )
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

  all_columns <- colnames(data_mat)
  m_is_missing <- missing(m)
  norm <- normalize_columns(columns, all_columns, match, m = m, m_is_missing = m_is_missing)

  x_cols <- norm$x
  mediator_cols <- norm$mediators
  y_cols <- norm$y
  m <- norm$m

  if (!is.character(output_dir) || length(output_dir) != 1 || is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be a non-empty character scalar.")
  }
  if (!is.character(output_prefix) || length(output_prefix) != 1 || is.na(output_prefix) || !nzchar(output_prefix)) {
    stop("output_prefix must be a non-empty character scalar.")
  }

  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1 || is.na(output_file) || !nzchar(output_file)) {
      stop("output_file must be NULL or a non-empty character scalar.")
    }
    output_dir <- dirname(output_file)
    output_prefix <- tools::file_path_sans_ext(basename(output_file))
    warning("output_file is deprecated; use output_dir + output_prefix.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("output_dir does not exist: ", output_dir)
  }

  if (!is.character(output) || anyNA(output)) {
    stop("output must be a character vector.")
  }
  allowed_outputs <- c("fit", "params", "effects")
  bad_outputs <- setdiff(unique(output), allowed_outputs)
  if (length(bad_outputs) > 0) {
    stop("Unknown output value(s): ", paste(bad_outputs, collapse = ", "))
  }
  write_fit <- "fit" %in% output
  write_params <- "params" %in% output
  write_effects <- "effects" %in% output

  derive_output_files <- function(output_dir, output_prefix) {
    list(
      fit = file.path(output_dir, paste0(output_prefix, "_fit.csv")),
      params = file.path(output_dir, paste0(output_prefix, "_params.csv")),
      effects = file.path(output_dir, paste0(output_prefix, "_effects.csv"))
    )
  }
  out_files <- derive_output_files(output_dir, output_prefix)

  selected_files <- c(
    if (write_fit) out_files$fit,
    if (write_params) out_files$params,
    if (write_effects) out_files$effects
  )

  if (!overwrite) {
    existing <- selected_files[file.exists(selected_files)]
    if (length(existing) > 0) {
      stop("Output file(s) already exist and overwrite=FALSE: ", paste(existing, collapse = ", "))
    }
  }

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

  if (!is.numeric(combination_start) || length(combination_start) != 1 || is.na(combination_start) ||
      !is.finite(combination_start) || combination_start < 0 || combination_start != floor(combination_start)) {
    stop("combination_start must be a single, non-negative whole number.")
  }
  max_exact <- 2^53 - 1
  if (combination_start > max_exact) {
    stop("combination_start must be <= 2^53-1 for exact integer behavior in R.")
  }
  if (!is.null(combination_end)) {
    if (!is.numeric(combination_end) || length(combination_end) != 1 || is.na(combination_end) ||
        !is.finite(combination_end) || combination_end < 0 || combination_end != floor(combination_end)) {
      stop("combination_end must be NULL or a single, non-negative whole number.")
    }
    if (combination_end > max_exact) {
      stop("combination_end must be <= 2^53-1 for exact integer behavior in R.")
    }
  }
  if (!is.null(combination_end) && combination_end < combination_start) {
    stop("combination_end must be >= combination_start.")
  }
  if (!is.numeric(shard_id) || length(shard_id) != 1 || is.na(shard_id) || shard_id < 0 || shard_id != floor(shard_id)) {
    stop("shard_id must be a single non-negative integer.")
  }
  if (!is.numeric(shard_count) || length(shard_count) != 1 || is.na(shard_count) || shard_count < 1 || shard_count != floor(shard_count)) {
    stop("shard_count must be a single positive integer.")
  }
  shard_id <- as.numeric(shard_id)
  shard_count <- as.numeric(shard_count)
  if (shard_id >= shard_count) {
    stop("shard_id must be in [0, shard_count).")
  }

  if (!is.numeric(warn_combinations) || length(warn_combinations) != 1 || is.na(warn_combinations) || warn_combinations < 0) {
    stop("warn_combinations must be a single non-negative number.")
  }
  if (!is.numeric(max_combinations) || length(max_combinations) != 1 || is.na(max_combinations) || max_combinations < 0) {
    stop("max_combinations must be a single non-negative number or Inf.")
  }

  dims <- c(length(x_cols), vapply(mediator_cols, length, integer(1)), length(y_cols))
  global_total <- prod(as.double(dims))

  # Preflight: estimate number of combinations processed by this call (range + sharding).
  end <- if (is.null(combination_end)) global_total else as.double(combination_end)
  start <- as.double(combination_start)

  n_to_process <- 0
  if (is.finite(global_total) && start <= end && end > start) {
    start_mod <- start %% shard_count
    delta <- (shard_id - start_mod) %% shard_count
    first <- start + delta
    if (first < end) {
      n_to_process <- 1 + floor((end - 1 - first) / shard_count)
    }
  } else if (!is.finite(global_total) && is.null(combination_end)) {
    stop("Unable to compute combination counts: global_total overflowed; specify combination_end explicitly.")
  }

  if (is.finite(n_to_process) && n_to_process > warn_combinations) {
    warning(sprintf(
      "Large analysis: global_total=%.0f; n_to_process=%.0f; dims=[%s].",
      global_total, n_to_process, paste(dims, collapse = " x ")
    ), call. = FALSE)
  }
  if (is.finite(n_to_process) && is.finite(max_combinations) && n_to_process > max_combinations) {
    stop(sprintf(
      "Refusing to run: n_to_process=%.0f exceeds max_combinations=%.0f (set max_combinations=Inf to override).",
      n_to_process, max_combinations
    ))
  }

  serial_path_analysis_cpp(
    data = data_mat,
    column_names = colnames(data_mat),
    x_col_idx = match(x_cols, colnames(data_mat)) - 1L,
    mediator_col_idx_list = lapply(mediator_cols, function(cols) match(cols, colnames(data_mat)) - 1L),
    y_col_idx = match(y_cols, colnames(data_mat)) - 1L,
    m = m,
    nrep = nrep,
    output_fit_file = out_files$fit,
    output_params_file = out_files$params,
    output_effects_file = out_files$effects,
    base_seed = base_seed,
    alpha = alpha,
    combination_start = as.double(combination_start),
    combination_end = if (is.null(combination_end)) NA_real_ else as.double(combination_end),
    shard_id = as.integer(shard_id),
    shard_count = as.integer(shard_count),
    chunk_size = chunk_size,
    grain_size = grain_size,
    overwrite = overwrite,
    excel_safe_csv = excel_safe_csv,
    write_fit = write_fit,
    write_params = write_params,
    write_effects = write_effects,
    include_failure_reason = include_failure_reason
  )

  created <- selected_files[file.exists(selected_files)]
  if (length(created) > 0) {
    cat("Serial path analysis completed. Results saved to:\n", paste(created, collapse = "\n "), "\n")
  } else {
    cat("Serial path analysis completed.\n")
  }
}
