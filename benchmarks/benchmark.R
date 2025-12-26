# Benchmark script for fastmed performance tuning
#
# Run from repo root (after installation):
#   Rscript --vanilla benchmarks/benchmark.R
#
# Or profile a single configuration:
#   Rscript --vanilla benchmarks/benchmark.R --profile

args <- commandArgs(trailingOnly = TRUE)
do_profile <- any(args %in% "--profile")

if (!requireNamespace("fastmed", quietly = TRUE)) {
  stop("fastmed must be installed. Run: R CMD INSTALL .")
}

get_rss_kb <- function() {
  pid <- Sys.getpid()
  cmd <- sprintf("ps -o rss= -p %d", pid)
  out <- suppressWarnings(system(cmd, intern = TRUE))
  rss <- suppressWarnings(as.numeric(trimws(out[1])))
  if (!is.finite(rss)) {
    return(NA_real_)
  }
  rss
}

generate_benchmark_data <- function(n, n_exposure, n_mediator, n_outcome, seed = 1) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  dat <- data.frame(matrix(0, nrow = n, ncol = 0))

  for (j in seq_len(n_exposure)) {
    dat[[sprintf("EXP%d", j)]] <- rbinom(n, 1, 0.5)
  }
  for (j in seq_len(n_mediator)) {
    exp_j <- dat[[sprintf("EXP%d", ((j - 1) %% n_exposure) + 1)]]
    dat[[sprintf("MED%d", j)]] <- 0.2 + 0.6 * exp_j + rnorm(n)
  }
  for (j in seq_len(n_outcome)) {
    exp_j <- dat[[sprintf("EXP%d", ((j - 1) %% n_exposure) + 1)]]
    med_j <- dat[[sprintf("MED%d", ((j - 1) %% n_mediator) + 1)]]
    dat[[sprintf("OUT%d", j)]] <- 0.1 + 0.5 * med_j + 0.3 * exp_j + rnorm(n)
  }

  dat
}

run_once <- function(label,
                     n,
                     n_exposure,
                     n_mediator,
                     n_outcome,
                     nrep,
                     num_threads,
                     chunk_size,
                     grain_size) {
  data <- generate_benchmark_data(
    n = n,
    n_exposure = n_exposure,
    n_mediator = n_mediator,
    n_outcome = n_outcome,
    seed = 1
  )

  out_file <- tempfile(fileext = ".csv")

  rss_before <- get_rss_kb()

  prof_file <- tempfile(fileext = ".out")
  if (isTRUE(do_profile)) {
    Rprof(prof_file, interval = 0.01)
    on.exit(Rprof(NULL), add = TRUE)
  }

  t <- system.time({
    fastmed::mediation_analysis(
      data = data,
      columns = list(exposure = "EXP", mediator = "MED", outcome = "OUT"),
      nrep = nrep,
      output_file = out_file,
      num_threads = num_threads,
      pert = "asymptotic",
      seed = 1,
      chunk_size = chunk_size,
      grain_size = grain_size,
      mediator.family = "gaussian",
      outcome.family = "gaussian"
    )
  })

  if (isTRUE(do_profile)) {
    Rprof(NULL)
  }

  rss_after <- get_rss_kb()

  res <- data.frame(
    label = label,
    n = n,
    n_exposure = n_exposure,
    n_mediator = n_mediator,
    n_outcome = n_outcome,
    combinations = n_exposure * n_mediator * n_outcome,
    nrep = nrep,
    num_threads = num_threads,
    chunk_size = chunk_size,
    grain_size = grain_size,
    elapsed = unname(t["elapsed"]),
    rss_kb_before = rss_before,
    rss_kb_after = rss_after
  )

  if (isTRUE(do_profile)) {
    cat("\nTop functions by total time (Rprof):\n")
    print(utils::head(summaryRprof(prof_file)$by.total, 20))
  }

  res
}

run_benchmarks <- function() {
  results <- list()
  i <- 0L

  configs <- list(
    list(label = "small_many", n = 200, n_exposure = 10, n_mediator = 10, n_outcome = 10, nrep = 200),
    list(label = "med", n = 1000, n_exposure = 5, n_mediator = 5, n_outcome = 5, nrep = 500),
    list(label = "large_n", n = 5000, n_exposure = 2, n_mediator = 2, n_outcome = 2, nrep = 500)
  )

  for (cfg in configs) {
    i <- i + 1L
    results[[i]] <- run_once(
      label = cfg$label,
      n = cfg$n,
      n_exposure = cfg$n_exposure,
      n_mediator = cfg$n_mediator,
      n_outcome = cfg$n_outcome,
      nrep = cfg$nrep,
      num_threads = max(1L, parallel::detectCores(logical = FALSE)),
      chunk_size = 10000L,
      grain_size = 100L
    )
  }

  do.call(rbind, results)
}

cat("Running fastmed benchmarks...\n")
bench <- run_benchmarks()
print(bench)