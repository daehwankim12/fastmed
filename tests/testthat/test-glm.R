test_that("glm_fit_cpp matches R glm for Gaussian", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- 2 + 3 * x + rnorm(n)

  r_fit <- glm(y ~ x, family = gaussian())
  cpp_fit <- fastmed:::glm_fit_cpp(cbind(1, x), y, "gaussian")

  expect_equal(unname(cpp_fit$coef), unname(coef(r_fit)), tolerance = 1e-6)
  expect_equal(unname(cpp_fit$vcov), unname(vcov(r_fit)), tolerance = 1e-6)
})

test_that("glm_fit_cpp matches R glm for Binomial", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  p <- plogis(1 + 2 * x)
  y <- rbinom(n, 1, p)

  r_fit <- glm(y ~ x, family = binomial())
  cpp_fit <- fastmed:::glm_fit_cpp(cbind(1, x), y, "binomial")

  expect_equal(unname(cpp_fit$coef), unname(coef(r_fit)), tolerance = 1e-4)
  expect_equal(unname(cpp_fit$vcov), unname(vcov(r_fit)), tolerance = 1e-4)
})

test_that("glm_fit_cpp matches R glm for Poisson", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  mu <- exp(0.2 + 0.5 * x)
  y <- rpois(n, mu)

  r_fit <- glm(y ~ x, family = poisson())
  cpp_fit <- fastmed:::glm_fit_cpp(cbind(1, x), y, "poisson")

  expect_equal(unname(cpp_fit$coef), unname(coef(r_fit)), tolerance = 1e-5)
  expect_equal(unname(cpp_fit$vcov), unname(vcov(r_fit)), tolerance = 1e-4)
})

test_that("glm_fit_cpp auto-detects Binomial from y in {0,1}", {
  set.seed(1)
  n <- 150
  x <- rnorm(n)
  p <- plogis(-0.3 + 1.1 * x)
  y <- rbinom(n, 1, p)

  r_fit <- glm(y ~ x, family = binomial())
  cpp_fit <- fastmed:::glm_fit_cpp(cbind(1, x), y, "auto")

  expect_equal(unname(cpp_fit$coef), unname(coef(r_fit)), tolerance = 1e-4)
})

test_that("glm_fit_cpp errors on rank deficiency", {
  set.seed(1)
  n <- 50
  x <- rnorm(n)
  X <- cbind(1, x, x) # duplicate column
  y <- 1 + 2 * x + rnorm(n)

  expect_error(fastmed:::glm_fit_cpp(X, y, "gaussian"), "rank deficient")
})