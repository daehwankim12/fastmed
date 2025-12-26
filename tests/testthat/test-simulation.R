test_that("mediate-style simulation matches an R reference implementation (gaussian mediator, binomial outcome)", {
  detect_family <- function(y) {
    eps <- 1e-8
    if (all(abs(y) <= eps | abs(y - 1) <= eps)) return("binomial")
    if (all(y >= -eps & abs(y - round(y)) <= eps)) return("poisson")
    "gaussian"
  }

  safe_exp <- function(x) exp(pmin(pmax(x, -700), 700))

  linkinv <- function(eta, fam) {
    switch(fam,
      gaussian = eta,
      binomial = stats::plogis(eta),
      poisson = safe_exp(eta),
      stop("Unknown family: ", fam)
    )
  }

  simulate_reference <- function(data, sims, seed) {
    fam_m <- detect_family(data$MED1)
    fam_y <- detect_family(data$OUT1)

    fit_m <- stats::glm(MED1 ~ EXP1,
      data = data,
      family = switch(fam_m,
        gaussian = stats::gaussian(),
        binomial = stats::binomial(),
        poisson = stats::poisson()
      )
    )
    fit_y <- stats::glm(OUT1 ~ MED1 + EXP1,
      data = data,
      family = switch(fam_y,
        gaussian = stats::gaussian(),
        binomial = stats::binomial(),
        poisson = stats::poisson()
      )
    )

    beta_m <- stats::coef(fit_m)
    beta_y <- stats::coef(fit_y)

    vcov_m <- stats::vcov(fit_m)
    vcov_y <- stats::vcov(fit_y)

    Lm <- t(chol(vcov_m))
    Ly <- t(chol(vcov_y))

    sigma2_m <- summary(fit_m)$dispersion
    n <- nrow(data)

    d0 <- numeric(sims)
    z0 <- numeric(sims)
    tau <- numeric(sims)

    set.seed(seed)
    for (s in seq_len(sims)) {
      bm <- as.numeric(beta_m + Lm %*% stats::rnorm(length(beta_m)))
      by <- as.numeric(beta_y + Ly %*% stats::rnorm(length(beta_y)))

      eta_m0 <- bm[1] + bm[2] * 0
      eta_m1 <- bm[1] + bm[2] * 1

      if (fam_m == "gaussian") {
        e <- stats::rnorm(n, 0, sqrt(max(0, sigma2_m)))
        M0 <- eta_m0 + e
        M1 <- eta_m1 + e
      } else if (fam_m == "binomial") {
        p0 <- stats::plogis(eta_m0)
        p1 <- stats::plogis(eta_m1)
        M0 <- stats::rbinom(n, 1, p0)
        M1 <- stats::rbinom(n, 1, p1)
      } else {
        lam0 <- safe_exp(eta_m0)
        lam1 <- safe_exp(eta_m1)
        M0 <- stats::rpois(n, lam0)
        M1 <- stats::rpois(n, lam1)
      }

      eta00 <- by[1] + by[2] * M0 + by[3] * 0
      eta01 <- by[1] + by[2] * M1 + by[3] * 0
      eta10 <- by[1] + by[2] * M0 + by[3] * 1
      eta11 <- by[1] + by[2] * M1 + by[3] * 1

      y00 <- linkinv(eta00, fam_y)
      y01 <- linkinv(eta01, fam_y)
      y10 <- linkinv(eta10, fam_y)
      y11 <- linkinv(eta11, fam_y)

      mean_y00 <- mean(y00)
      mean_y01 <- mean(y01)
      mean_y10 <- mean(y10)
      mean_y11 <- mean(y11)

      d0[s] <- mean_y01 - mean_y00
      z0[s] <- mean_y10 - mean_y00
      tau[s] <- mean_y11 - mean_y00
    }

    list(
      acme = mean(d0),
      ade = mean(z0),
      total = mean(tau)
    )
  }

  set.seed(42)
  n <- 300
  T <- stats::rbinom(n, 1, 0.5)
  M <- 0.2 + 0.6 * T + stats::rnorm(n, sd = 1.0)
  pY <- stats::plogis(-0.3 + 0.7 * M + 0.4 * T)
  Y <- stats::rbinom(n, 1, pY)

  test_data <- data.frame(EXP1 = T, MED1 = M, OUT1 = Y)

  sims <- 400
  seed <- 123

  ref <- simulate_reference(test_data, sims = sims, seed = seed)

  output_csv <- withr::local_tempfile(fileext = ".csv")
  mediation_analysis(
    data = test_data,
    columns = list(exposure = c("EXP"), mediator = c("MED"), outcome = c("OUT")),
    nrep = sims,
    output_file = output_csv,
    num_threads = 1,
    pert = "asymptotic",
    seed = seed
  )

  results <- data.table::fread(output_csv)
  expect_equal(nrow(results), 1)

  expect_equal(results$ACME_Mean, ref$acme, tolerance = 0.10)
  expect_equal(results$ADE_Mean, ref$ade, tolerance = 0.10)
  expect_equal(results$Total_Effect_Mean, ref$total, tolerance = 0.10)
})

