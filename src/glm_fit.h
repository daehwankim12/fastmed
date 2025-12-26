#pragma once

#include <RcppEigen.h>

#include <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;

enum class GlmFamily { Gaussian, Binomial, Poisson };

struct GlmFit {
    VectorXd coef;      // p
    MatrixXd vcov;      // p x p
    double dispersion;  // gaussian: sigma^2, else 1
    bool converged;
    int iterations;
    double dev;
    int rank;
};

GlmFamily parse_family_or_auto(const std::string& fam_str, const VectorXd& y);

GlmFit glm_fit_irls_qr(const MatrixXd& X,
                       const VectorXd& y,
                       GlmFamily fam,
                       const VectorXd& prior_w,
                       const VectorXd& offset,
                       int maxit,
                       double epsilon,
                       double qr_tol);

