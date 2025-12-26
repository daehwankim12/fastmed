#include "glm_fit.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <vector>

static inline std::string lower_copy(std::string s) {
    std::transform(s.begin(),
                   s.end(),
                   s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static inline double inv_logit(double x) {
    if (x >= 0.0) {
        double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    double z = std::exp(x);
    return z / (1.0 + z);
}

static inline double safe_exp(double x) {
    if (x > 700.0) {
        x = 700.0;
    }
    if (x < -700.0) {
        x = -700.0;
    }
    return std::exp(x);
}

static inline double linkinv(double eta, GlmFamily fam) {
    if (fam == GlmFamily::Gaussian) {
        return eta;
    }
    if (fam == GlmFamily::Binomial) {
        return inv_logit(eta);
    }
    return safe_exp(eta); // Poisson
}

static inline double linkfun(double mu, GlmFamily fam) {
    const double eps = 1e-15;
    if (fam == GlmFamily::Gaussian) {
        return mu;
    }
    if (fam == GlmFamily::Binomial) {
        double m = std::min(std::max(mu, eps), 1.0 - eps);
        return std::log(m / (1.0 - m));
    }
    // Poisson
    double m = std::max(mu, eps);
    return std::log(m);
}

static inline double variance_mu(double mu, GlmFamily fam) {
    if (fam == GlmFamily::Gaussian) {
        return 1.0;
    }
    if (fam == GlmFamily::Binomial) {
        return mu * (1.0 - mu);
    }
    return mu; // Poisson
}

static inline double mu_eta(double mu, GlmFamily fam) {
    if (fam == GlmFamily::Gaussian) {
        return 1.0;
    }
    if (fam == GlmFamily::Binomial) {
        return mu * (1.0 - mu);
    }
    return mu; // Poisson
}

static inline bool valid_mu(double mu, GlmFamily fam) {
    if (!std::isfinite(mu)) {
        return false;
    }
    if (fam == GlmFamily::Gaussian) {
        return true;
    }
    if (fam == GlmFamily::Binomial) {
        return (mu > 0.0 && mu < 1.0);
    }
    return (mu > 0.0); // Poisson
}

static inline void validate_y_for_family(const VectorXd& y, GlmFamily fam) {
    const double eps = 1e-8;
    if (fam == GlmFamily::Binomial) {
        for (int i = 0; i < y.size(); ++i) {
            double v = y[i];
            if (!std::isfinite(v)) {
                throw std::runtime_error("Binomial y contains non-finite");
            }
            if (std::fabs(v) <= eps) {
                continue;
            }
            if (std::fabs(v - 1.0) <= eps) {
                continue;
            }
            throw std::runtime_error("Binomial requires y in {0,1} (found non-binary).");
        }
    } else if (fam == GlmFamily::Poisson) {
        for (int i = 0; i < y.size(); ++i) {
            double v = y[i];
            if (!std::isfinite(v)) {
                throw std::runtime_error("Poisson y contains non-finite");
            }
            if (v < -eps) {
                throw std::runtime_error("Poisson requires y >= 0.");
            }
            double r = std::round(v);
            if (std::fabs(v - r) > eps) {
                throw std::runtime_error("Poisson requires integer y.");
            }
        }
    }
}

static inline GlmFamily detect_family_from_y(const VectorXd& y) {
    const double eps = 1e-8;
    bool all01 = true;
    bool all_nonneg_int = true;
    for (int i = 0; i < y.size(); ++i) {
        double v = y[i];
        if (!std::isfinite(v)) {
            throw std::runtime_error("Auto-detect: non-finite y");
        }
        if (!(std::fabs(v) <= eps || std::fabs(v - 1.0) <= eps)) {
            all01 = false;
        }
        if (v < -eps) {
            all_nonneg_int = false;
        } else {
            double r = std::round(v);
            if (std::fabs(v - r) > eps) {
                all_nonneg_int = false;
            }
        }
    }
    if (all01) {
        return GlmFamily::Binomial;
    }
    if (all_nonneg_int) {
        return GlmFamily::Poisson;
    }
    return GlmFamily::Gaussian;
}

GlmFamily parse_family_or_auto(const std::string& fam_str, const VectorXd& y) {
    std::string f = lower_copy(fam_str);
    if (f == "auto") {
        return detect_family_from_y(y);
    }
    if (f == "gaussian" || f == "normal") {
        return GlmFamily::Gaussian;
    }
    if (f == "binomial" || f == "bernoulli" || f == "logistic") {
        return GlmFamily::Binomial;
    }
    if (f == "poisson") {
        return GlmFamily::Poisson;
    }
    throw std::invalid_argument("Unknown family: " + fam_str);
}

static inline double deviance_glm(const VectorXd& y,
                                 const VectorXd& mu,
                                 const VectorXd& w,
                                 GlmFamily fam) {
    const double eps = 1e-15;
    double dev = 0.0;

    if (fam == GlmFamily::Gaussian) {
        dev = (w.array() * (y - mu).array().square()).sum();
        return dev;
    }

    if (fam == GlmFamily::Binomial) {
        for (int i = 0; i < y.size(); ++i) {
            double yi = y[i];
            double mui = std::min(std::max(mu[i], eps), 1.0 - eps);
            double wi = w[i];

            if (yi <= 0.0) {
                dev += 2.0 * wi * (-std::log(1.0 - mui));
            } else if (yi >= 1.0) {
                dev += 2.0 * wi * (-std::log(mui));
            } else {
                dev += 2.0 * wi * (yi * std::log(yi / mui) +
                                   (1.0 - yi) * std::log((1.0 - yi) / (1.0 - mui)));
            }
        }
        return dev;
    }

    // Poisson
    for (int i = 0; i < y.size(); ++i) {
        double yi = y[i];
        double mui = std::max(mu[i], eps);
        double wi = w[i];
        if (yi == 0.0) {
            dev += 2.0 * wi * mui;
        } else {
            dev += 2.0 * wi * (yi * std::log(yi / mui) - (yi - mui));
        }
    }
    return dev;
}

static inline VectorXd wls_solve_qr(const MatrixXd& X,
                                    const VectorXd& z,
                                    const VectorXd& sqrt_w,
                                    double qr_tol,
                                    int* out_rank = nullptr) {
    if (X.rows() != z.size() || z.size() != sqrt_w.size()) {
        throw std::runtime_error("wls_solve_qr: dimension mismatch");
    }

    MatrixXd Xw = X.array().colwise() * sqrt_w.array();
    VectorXd zw = z.array() * sqrt_w.array();

    Eigen::ColPivHouseholderQR<MatrixXd> qr(Xw);
    qr.setThreshold(qr_tol);
    if (out_rank) {
        *out_rank = qr.rank();
    }
    return qr.solve(zw);
}

GlmFit glm_fit_irls_qr(const MatrixXd& X,
                       const VectorXd& y,
                       GlmFamily fam,
                       const VectorXd& prior_w,
                       const VectorXd& offset,
                       int maxit,
                       double epsilon,
                       double qr_tol) {
    const int n = static_cast<int>(X.rows());
    const int p = static_cast<int>(X.cols());
    if (y.size() != n) {
        throw std::runtime_error("glm_fit_irls_qr: X/y mismatch");
    }
    if (prior_w.size() != n) {
        throw std::runtime_error("glm_fit_irls_qr: weights mismatch");
    }
    if (offset.size() != n) {
        throw std::runtime_error("glm_fit_irls_qr: offset mismatch");
    }
    if (n <= p) {
        throw std::runtime_error("glm_fit_irls_qr: n must be > p");
    }

    if (fam != GlmFamily::Gaussian) {
        validate_y_for_family(y, fam);
    }

    VectorXd mustart(n);
    if (fam == GlmFamily::Gaussian) {
        mustart = y;
    } else if (fam == GlmFamily::Binomial) {
        mustart = (prior_w.array() * y.array() + 0.5) /
                  (prior_w.array() + 1.0);
        mustart = mustart.array().min(1.0 - 1e-8).max(1e-8);
    } else { // Poisson
        mustart = (y.array() + 0.1).max(1e-8);
    }

    VectorXd eta(n);
    VectorXd mu(n);
    for (int i = 0; i < n; ++i) {
        eta[i] = linkfun(mustart[i], fam);
        mu[i] = linkinv(eta[i], fam);
    }

    VectorXd beta(p);
    {
        std::vector<int> good;
        good.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (prior_w[i] > 0 && std::isfinite(eta[i]) && std::isfinite(offset[i])) {
                good.push_back(i);
            }
        }
        if (static_cast<int>(good.size()) <= p) {
            throw std::runtime_error("glm_fit_irls_qr: insufficient good rows at init");
        }

        MatrixXd Xg(static_cast<int>(good.size()), p);
        VectorXd zg(static_cast<int>(good.size()));
        VectorXd sw(static_cast<int>(good.size()));
        for (int r = 0; r < static_cast<int>(good.size()); ++r) {
            int i = good[r];
            Xg.row(r) = X.row(i);
            zg[r] = eta[i] - offset[i];
            sw[r] = std::sqrt(prior_w[i]);
        }

        int rank = 0;
        beta = wls_solve_qr(Xg, zg, sw, qr_tol, &rank);
        if (rank < p) {
            throw std::runtime_error("glm_fit_irls_qr: rank deficient at init");
        }
    }

    double dev_old = deviance_glm(y, mu, prior_w, fam);

    bool converged = false;
    int it = 0;
    int rank = p;

    VectorXd z(n);
    VectorXd w_sqrt(n);
    VectorXd mu_eta_v(n);
    VectorXd var_v(n);

    for (it = 0; it < maxit; ++it) {
        eta = X * beta + offset;

        for (int i = 0; i < n; ++i) {
            double mui = linkinv(eta[i], fam);
            if (fam == GlmFamily::Binomial) {
                mui = std::min(std::max(mui, 1e-12), 1.0 - 1e-12);
            } else if (fam == GlmFamily::Poisson) {
                mui = std::max(mui, 1e-12);
            }
            mu[i] = mui;
            mu_eta_v[i] = mu_eta(mui, fam);
            var_v[i] = variance_mu(mui, fam);
        }

        std::vector<int> good;
        good.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (!(prior_w[i] > 0)) {
                continue;
            }
            if (!std::isfinite(eta[i]) || !std::isfinite(mu[i])) {
                continue;
            }

            double d = std::max(mu_eta_v[i], 1e-12);
            double v = std::max(var_v[i], 1e-12);

            z[i] = (eta[i] - offset[i]) + (y[i] - mu[i]) / d;
            w_sqrt[i] = std::sqrt(prior_w[i] * (d * d) / v);

            if (std::isfinite(z[i]) && std::isfinite(w_sqrt[i]) && w_sqrt[i] > 0) {
                good.push_back(i);
            }
        }

        if (static_cast<int>(good.size()) <= p) {
            throw std::runtime_error("glm_fit_irls_qr: insufficient good rows");
        }

        MatrixXd Xg(static_cast<int>(good.size()), p);
        VectorXd zg(static_cast<int>(good.size()));
        VectorXd sw(static_cast<int>(good.size()));
        for (int r = 0; r < static_cast<int>(good.size()); ++r) {
            int i = good[r];
            Xg.row(r) = X.row(i);
            zg[r] = z[i];
            sw[r] = w_sqrt[i];
        }

        VectorXd beta_new;
        {
            int rrank = 0;
            beta_new = wls_solve_qr(Xg, zg, sw, qr_tol, &rrank);
            rank = rrank;
            if (rank < p) {
                throw std::runtime_error("glm_fit_irls_qr: rank deficient");
            }
            if (!beta_new.allFinite()) {
                throw std::runtime_error("glm_fit_irls_qr: non-finite coefficients");
            }
        }

        auto dev_for_beta = [&](const VectorXd& b, bool* out_valid_mu) -> double {
            VectorXd e = X * b + offset;
            VectorXd m(n);
            bool ok = true;
            for (int i = 0; i < n; ++i) {
                double mui = linkinv(e[i], fam);
                if (fam == GlmFamily::Binomial) {
                    mui = std::min(std::max(mui, 1e-12), 1.0 - 1e-12);
                } else if (fam == GlmFamily::Poisson) {
                    mui = std::max(mui, 1e-12);
                }
                m[i] = mui;
                if (!valid_mu(mui, fam)) {
                    ok = false;
                }
            }
            if (out_valid_mu) {
                *out_valid_mu = ok;
            }
            return deviance_glm(y, m, prior_w, fam);
        };

        bool mu_ok_new = true;
        double dev_new = dev_for_beta(beta_new, &mu_ok_new);

        int half = 0;
        while ((!std::isfinite(dev_new) || dev_new > dev_old || !mu_ok_new) && half < 25) {
            beta_new = 0.5 * (beta_new + beta);
            dev_new = dev_for_beta(beta_new, &mu_ok_new);
            ++half;
        }

        beta = beta_new;

        if (std::isfinite(dev_new)) {
            double denom = 0.1 + std::fabs(dev_new);
            if (std::fabs(dev_new - dev_old) / denom < epsilon) {
                converged = true;
                dev_old = dev_new;
                ++it;
                break;
            }
            dev_old = dev_new;
        } else {
            throw std::runtime_error("glm_fit_irls_qr: deviance is not finite");
        }
    }

    eta = X * beta + offset;
    for (int i = 0; i < n; ++i) {
        double mui = linkinv(eta[i], fam);
        if (fam == GlmFamily::Binomial) {
            mui = std::min(std::max(mui, 1e-12), 1.0 - 1e-12);
        } else if (fam == GlmFamily::Poisson) {
            mui = std::max(mui, 1e-12);
        }
        mu[i] = mui;
    }

    double phi = 1.0;
    if (fam == GlmFamily::Gaussian) {
        int df = n - rank;
        if (df <= 0) {
            throw std::runtime_error("glm_fit_irls_qr: df_resid <= 0");
        }
        phi = deviance_glm(y, mu, prior_w, fam) / static_cast<double>(df);
        phi = std::max(phi, 1e-15);
    }

    VectorXd W(n);
    for (int i = 0; i < n; ++i) {
        double d = std::max(mu_eta(mu[i], fam), 1e-12);
        double v = std::max(variance_mu(mu[i], fam), 1e-12);
        W[i] = prior_w[i] * (d * d) / v;
    }

    MatrixXd XtWX = X.transpose() * W.asDiagonal() * X;
    XtWX = 0.5 * (XtWX + XtWX.transpose());
    Eigen::LDLT<MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("glm_fit_irls_qr: XtWX LDLT failed");
    }

    MatrixXd vcov = ldlt.solve(MatrixXd::Identity(p, p)) * phi;

    GlmFit fit;
    fit.coef = beta;
    fit.vcov = vcov;
    fit.dispersion = phi;
    fit.converged = converged;
    fit.iterations = it;
    fit.dev = dev_old;
    fit.rank = rank;
    return fit;
}

// [[Rcpp::export]]
Rcpp::List glm_fit_cpp(const Eigen::MatrixXd& X,
                       const Eigen::VectorXd& y,
                       std::string family = "auto",
                       Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue,
                       Rcpp::Nullable<Rcpp::NumericVector> offset = R_NilValue,
                       int maxit = 25,
                       double epsilon = 1e-8,
                       double qr_tol = 1e-12) {
    if (X.rows() != y.size()) {
        throw std::runtime_error("glm_fit_cpp: X/y mismatch");
    }
    const int n = static_cast<int>(X.rows());

    Eigen::VectorXd prior_w(n);
    if (weights.isNotNull()) {
        Rcpp::NumericVector w(weights);
        if (w.size() != n) {
            throw std::runtime_error("glm_fit_cpp: weights length mismatch");
        }
        prior_w = Rcpp::as<Eigen::VectorXd>(w);
    } else {
        prior_w.setOnes();
    }

    Eigen::VectorXd off(n);
    if (offset.isNotNull()) {
        Rcpp::NumericVector o(offset);
        if (o.size() != n) {
            throw std::runtime_error("glm_fit_cpp: offset length mismatch");
        }
        off = Rcpp::as<Eigen::VectorXd>(o);
    } else {
        off.setZero();
    }

    GlmFamily fam = parse_family_or_auto(family, y);
    GlmFit fit = glm_fit_irls_qr(X, y, fam, prior_w, off, maxit, epsilon, qr_tol);

    return Rcpp::List::create(
        Rcpp::Named("coef") = fit.coef,
        Rcpp::Named("vcov") = fit.vcov,
        Rcpp::Named("dispersion") = fit.dispersion,
        Rcpp::Named("converged") = fit.converged,
        Rcpp::Named("iterations") = fit.iterations,
        Rcpp::Named("dev") = fit.dev,
        Rcpp::Named("rank") = fit.rank);
}

