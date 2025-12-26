generate_mediation_data <- function(n = 100,
                                    treat_family = c("binary", "continuous"),
                                    mediator_family = c("gaussian", "binomial", "poisson"),
                                    outcome_family = c("gaussian", "binomial", "poisson"),
                                    true_acme = 0.3,
                                    true_ade = 0.4,
                                    seed = NULL) {
  treat_family <- match.arg(treat_family)
  mediator_family <- match.arg(mediator_family)
  outcome_family <- match.arg(outcome_family)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (treat_family == "binary") {
    T <- stats::rbinom(n, 1, 0.5)
  } else {
    T <- stats::rnorm(n)
  }

  if (mediator_family == "gaussian") {
    M <- 1 + 0.5 * T + stats::rnorm(n)
  } else if (mediator_family == "binomial") {
    p <- stats::plogis(-0.2 + 0.6 * T)
    M <- stats::rbinom(n, 1, p)
  } else if (mediator_family == "poisson") {
    lambda <- exp(0.2 + 0.3 * T)
    M <- stats::rpois(n, lambda)
  }

  standardize <- function(x) {
    sx <- stats::sd(x)
    if (!is.finite(sx) || sx <= 0) {
      return(x - mean(x))
    }
    (x - mean(x)) / sx
  }

  if (outcome_family == "gaussian") {
    Y <- 2 + true_acme * M + true_ade * T + stats::rnorm(n)
  } else if (outcome_family == "binomial") {
    lp <- -0.2 + true_acme * standardize(M) + true_ade * T
    p <- stats::plogis(lp)
    Y <- stats::rbinom(n, 1, p)
  } else if (outcome_family == "poisson") {
    lp <- 0.1 + true_acme * standardize(M) + true_ade * T
    lambda <- exp(lp)
    Y <- stats::rpois(n, lambda)
  }

  data.frame(T = T, M = M, Y = Y)
}

