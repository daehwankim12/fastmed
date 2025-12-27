test_that("fuzz-style randomized runs do not crash or hang", {
  trials <- suppressWarnings(as.integer(Sys.getenv("FASTMED_FUZZ_TRIALS", "200")))
  if (is.na(trials) || trials <= 0) {
    trials <- 200L
  }

  gen_response <- function(n, family_name) {
    if (family_name == "gaussian") {
      return(stats::rnorm(n))
    }
    if (family_name == "binomial") {
      return(stats::rbinom(n, size = 1L, prob = 0.5))
    }
    if (family_name == "poisson") {
      return(stats::rpois(n, lambda = 1.0))
    }
    stop("Unknown family: ", family_name)
  }

  for (i in seq_len(trials)) {
    set.seed(1000 + i)

    n <- sample.int(40, 1) + 10
    fam_m <- sample(c("gaussian", "binomial", "poisson"), 1, prob = c(0.4, 0.3, 0.3))
    fam_y <- sample(c("gaussian", "binomial", "poisson"), 1, prob = c(0.5, 0.25, 0.25))

    dt <- data.table::data.table(
      EXP1 = stats::rbinom(n, 1, 0.5),
      MED1 = gen_response(n, fam_m),
      OUT1 = gen_response(n, fam_y)
    )

    weights <- NULL
    if (stats::runif(1) < 0.25) {
      weights <- stats::runif(n)
      weights[sample.int(n, size = floor(n / 3))] <- 0
      if (sum(weights) <= 0) {
        weights[1] <- 1
      }
    }

    output_csv <- withr::local_tempfile(fileext = ".csv")
    expect_error(
      mediation_analysis(
        data = dt,
        columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
        nrep = 3,
        output_file = output_csv,
        num_threads = 1,
        pert = "asymptotic",
        seed = i,
        weights = weights,
        mediator.family = "auto",
        outcome.family = "auto",
        chunk_size = 1,
        grain_size = 1
      ),
      NA
    )

    results <- data.table::fread(output_csv)
    expect_equal(nrow(results), 1)
    expect_equal(results$Combination, "EXP1_MED1_OUT1")

    numeric_cols <- names(results)[vapply(results, is.numeric, logical(1))]
    for (col in numeric_cols) {
      x <- results[[col]]
      expect_true(all(is.finite(x) | is.na(x)))
    }
  }
})

