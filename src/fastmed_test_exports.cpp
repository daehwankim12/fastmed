// fastmed_test_exports.cpp
#include <Rcpp.h>
#include <RcppEigen.h>

#include "fastmed_rng.h"
#include "fastmed_stats.h"

#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::IntegerMatrix;
using Rcpp::List;
using Rcpp::NumericVector;

static VectorXd linear_regression(const MatrixXd& X, const VectorXd& y) {
    if (X.rows() != y.rows()) {
        throw std::runtime_error("Mismatch in number of rows between X and y");
    }

    Eigen::JacobiSVD<MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return svd.solve(y);
}

// [[Rcpp::export]]
double fastmed_test_p_value_cpp(NumericVector samples) {
    std::vector<double> vec = Rcpp::as<std::vector<double>>(samples);
    if (vec.empty()) {
        return 1.0;
    }
    const double est = mean_cpp(vec);
    return pval_mediate(vec, est);
}

// [[Rcpp::export]]
List fastmed_test_calculate_statistics_cpp(NumericVector samples) {
    std::vector<double> vec = Rcpp::as<std::vector<double>>(samples);
    StatisticsSummary stats = calculate_statistics_inplace(vec);
    return List::create(Rcpp::_["mean"] = stats.mean,
                        Rcpp::_["percentile_2_5"] = stats.percentile_2_5,
                        Rcpp::_["percentile_97_5"] = stats.percentile_97_5,
                        Rcpp::_["p_value"] = stats.p_value);
}

// [[Rcpp::export]]
IntegerMatrix fastmed_test_two_bootstrap_samples(int n,
                                                 uint64_t base_seed,
                                                 uint64_t global_combination_idx,
                                                 uint64_t rep_idx) {
    if (n <= 0) {
        throw std::invalid_argument("n must be positive");
    }

    std::mt19937_64 gen(derive_seed(base_seed, global_combination_idx, rep_idx));
    std::uniform_int_distribution<int> dis(0, n - 1);

    IntegerMatrix draws(n, 2);
    for (int i = 0; i < n; ++i) {
        draws(i, 0) = dis(gen);
    }
    for (int i = 0; i < n; ++i) {
        draws(i, 1) = dis(gen);
    }
    return draws;
}

// [[Rcpp::export]]
List fastmed_test_ols_sigmas_cpp(NumericVector exposure,
                                NumericVector mediator,
                                NumericVector outcome) {
    const auto n = static_cast<int>(exposure.size());
    if (mediator.size() != n || outcome.size() != n) {
        throw std::invalid_argument(
            "Mismatch in lengths between exposure, mediator, and outcome");
    }
    if (n <= 0) {
        throw std::invalid_argument("Vectors must be non-empty");
    }

    VectorXd exposure_eigen = Rcpp::as<VectorXd>(exposure);
    VectorXd mediator_eigen = Rcpp::as<VectorXd>(mediator);
    VectorXd outcome_eigen = Rcpp::as<VectorXd>(outcome);

    MatrixXd X_med(n, 2);
    X_med.col(0).setOnes();
    X_med.col(1) = exposure_eigen;
    VectorXd beta_med = linear_regression(X_med, mediator_eigen);

    MatrixXd X_out(n, 3);
    X_out.col(0).setOnes();
    X_out.col(1) = mediator_eigen;
    X_out.col(2) = exposure_eigen;
    VectorXd beta_out = linear_regression(X_out, outcome_eigen);

    const int p_med = static_cast<int>(X_med.cols());
    const int p_out = static_cast<int>(X_out.cols());
    if (n <= p_med) {
        throw std::invalid_argument("Insufficient observations for mediator model");
    }
    if (n <= p_out) {
        throw std::invalid_argument("Insufficient observations for outcome model");
    }

    VectorXd resid_med = mediator_eigen - X_med * beta_med;
    double sigma_med =
        std::sqrt(resid_med.squaredNorm() / static_cast<double>(n - p_med));

    VectorXd resid_out = outcome_eigen - X_out * beta_out;
    double sigma_out =
        std::sqrt(resid_out.squaredNorm() / static_cast<double>(n - p_out));

    return List::create(Rcpp::_["sigma_med"] = sigma_med,
                        Rcpp::_["sigma_out"] = sigma_out,
                        Rcpp::_["df_med"] = n - p_med,
                        Rcpp::_["df_out"] = n - p_out);
}

