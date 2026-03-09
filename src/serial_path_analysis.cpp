// serial_path_analysis.cpp
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

#include "fastmed_csv.h"
#include "fastmed_rng.h"

#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// [[Rcpp::depends(RcppEigen, RcppParallel)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;
using Rcpp::CharacterVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;

namespace {

constexpr double kMaxExactIntDouble = 9007199254740991.0;  // 2^53 - 1

enum class PathModelType { Partial, Semifull, Full };

struct OnlineVariance {
    int n = 0;
    double mean = 0.0;
    double m2 = 0.0;

    void add(double x) {
        ++n;
        const double delta = x - mean;
        mean += delta / static_cast<double>(n);
        const double delta2 = x - mean;
        m2 += delta * delta2;
    }

    double variance() const {
        if (n < 2) {
            return 0.0;
        }
        return m2 / static_cast<double>(n - 1);
    }
};

inline double safe_sqrt(double x) {
    if (!(x >= 0.0) || !std::isfinite(x)) {
        return NA_REAL;
    }
    return std::sqrt(x);
}

inline double normal_pvalue_two_sided(double z) {
    if (!std::isfinite(z)) {
        return NA_REAL;
    }
    const double absz = std::fabs(z);
    const double upper = 1.0 - R::pnorm(absz, 0.0, 1.0, 1, 0);
    double p = 2.0 * upper;
    if (p < 0.0) {
        p = 0.0;
    }
    if (p > 1.0) {
        p = 1.0;
    }
    return p;
}

inline double chisq_pvalue_upper(double chisq, int df) {
    if (!std::isfinite(chisq) || df < 0) {
        return NA_REAL;
    }
    if (df == 0) {
        return (chisq == 0.0 ? 1.0 : 0.0);
    }
    return 1.0 - R::pchisq(chisq, static_cast<double>(df), 1, 0);
}

inline bool logdet_spd(const Eigen::Ref<const MatrixXd>& mat, double* out_logdet) {
    Eigen::LLT<MatrixXd> llt(mat);
    if (llt.info() != Eigen::Success) {
        return false;
    }
    const auto& L = llt.matrixL();
    double sum_log = 0.0;
    for (int i = 0; i < L.rows(); ++i) {
        const double d = L(i, i);
        if (!(d > 0.0) || !std::isfinite(d)) {
            return false;
        }
        sum_log += std::log(d);
    }
    *out_logdet = 2.0 * sum_log;
    return true;
}

inline bool solve_spd(const Eigen::Ref<const MatrixXd>& mat,
                      const Eigen::Ref<const MatrixXd>& rhs,
                      MatrixXd* out) {
    Eigen::LLT<MatrixXd> llt(mat);
    if (llt.info() != Eigen::Success) {
        return false;
    }
    *out = llt.solve(rhs);
    if (llt.info() != Eigen::Success) {
        return false;
    }
    return out->allFinite();
}

inline bool regression_from_cov(const Eigen::Ref<const MatrixXd>& Sxx,
                                const Eigen::Ref<const VectorXd>& Sxy,
                                double Syy,
                                Eigen::Ref<VectorXd> beta_out,
                                double* resid_var_out) {
    Eigen::LDLT<MatrixXd> ldlt(Sxx);
    if (ldlt.info() != Eigen::Success) {
        return false;
    }
    beta_out = ldlt.solve(Sxy);
    if (ldlt.info() != Eigen::Success || !beta_out.allFinite()) {
        return false;
    }
    const double resid = Syy - Sxy.dot(beta_out);
    if (!(resid > 0.0) || !std::isfinite(resid)) {
        return false;
    }
    *resid_var_out = resid;
    return true;
}

inline int model_df(PathModelType model, int m) {
    switch (model) {
    case PathModelType::Partial:
        return 0;
    case PathModelType::Semifull:
        return 1;
    case PathModelType::Full: {
        const int64_t df64 = (static_cast<int64_t>(m) * static_cast<int64_t>(m + 1)) / 2;
        if (df64 > std::numeric_limits<int>::max()) {
            return NA_INTEGER;
        }
        return static_cast<int>(df64);
    }
    }
    return NA_INTEGER;
}

inline int model_param_count(PathModelType model, int m) {
    const int64_t mm = static_cast<int64_t>(m);
    const int64_t p = mm + 2;
    switch (model) {
    case PathModelType::Partial: {
        const int64_t npar = (p * (p + 1)) / 2;
        if (npar > std::numeric_limits<int>::max()) {
            return NA_INTEGER;
        }
        return static_cast<int>(npar);
    }
    case PathModelType::Semifull: {
        const int64_t npar = (p * (p + 1)) / 2 - 1;
        if (npar > std::numeric_limits<int>::max()) {
            return NA_INTEGER;
        }
        return static_cast<int>(npar);
    }
    case PathModelType::Full: {
        const int64_t npar = 2 * mm + 3;
        if (npar > std::numeric_limits<int>::max()) {
            return NA_INTEGER;
        }
        return static_cast<int>(npar);
    }
    }
    return NA_INTEGER;
}

struct DynamicPathCoefs {
    std::vector<double> a;  // size m
    std::vector<double> d;  // lower triangle, size m(m-1)/2
    std::vector<double> b;  // size m
    double cp = 0.0;

    void resize(int m) {
        a.assign(static_cast<size_t>(m), 0.0);
        b.assign(static_cast<size_t>(m), 0.0);
        const size_t d_size = static_cast<size_t>(m) * static_cast<size_t>(m - 1) / 2;
        d.assign(d_size, 0.0);
        cp = 0.0;
    }

    static size_t d_index(int j, int k) {
        // 1-indexed mediators: 1 <= k < j <= m
        return static_cast<size_t>((j - 1) * (j - 2) / 2 + (k - 1));
    }

    double& d_jk(int j, int k) {
        return d[d_index(j, k)];
    }
    double d_jk(int j, int k) const {
        return d[d_index(j, k)];
    }
};

struct DynamicPathFit {
    bool ok = false;
    int n_complete = 0;
    int m = 0;

    DynamicPathCoefs coefs;
    double var_x = NA_REAL;
    std::vector<double> var_e;  // size m+1: M1..Mm, Y

    MatrixXd sigma;  // p x p

    double chisq = NA_REAL;
    int df = NA_INTEGER;
    double pvalue = NA_REAL;
    double cfi = NA_REAL;
    double tli = NA_REAL;
    double rmsea = NA_REAL;
    double aic = NA_REAL;
};

struct WorkerBuffers {
    int m = 0;
    int p = 0;
    int max_pred = 0;

    MatrixXd S;        // p x p covariance
    MatrixXd Spp;      // max_pred x max_pred predictors covariance
    VectorXd Sxy;      // max_pred predictors-outcome cov
    VectorXd beta;     // max_pred regression coefficients

    VectorXd sums;     // p
    VectorXd means;    // p
    VectorXd centered; // p
    VectorXd sd;       // p

    MatrixXd A;      // p x p structural coefficients
    MatrixXd M;      // p x p (I-A)
    MatrixXd Omega;  // p x p diagonal residual cov
    MatrixXd invM;   // p x p (I-A)^{-1}

    // Wishart bootstrap
    MatrixXd L;    // p x p cholesky factor
    MatrixXd W_A;  // p x p Bartlett factor
    MatrixXd Srep; // p x p covariance draw

    std::vector<uint64_t> decoded;  // length m+2
    std::vector<int> var_cols;      // length p
    std::vector<std::string> m_names;

    size_t d_count = 0;
    size_t param_count_out = 0;
    size_t effect_count_out = 0;
    std::vector<OnlineVariance> param_vars;
    std::vector<OnlineVariance> effect_vars;

    void resize(int new_m) {
        m = new_m;
        p = m + 2;
        max_pred = m + 1;

        S.resize(p, p);
        Spp.resize(max_pred, max_pred);
        Sxy.resize(max_pred);
        beta.resize(max_pred);

        sums.resize(p);
        means.resize(p);
        centered.resize(p);
        sd.resize(p);

        A.resize(p, p);
        M.resize(p, p);
        Omega.resize(p, p);
        invM.resize(p, p);

        L.resize(p, p);
        W_A.resize(p, p);
        Srep.resize(p, p);

        decoded.assign(static_cast<size_t>(m + 2), 0);
        var_cols.assign(static_cast<size_t>(p), -1);
        m_names.clear();
        m_names.reserve(static_cast<size_t>(m));

        d_count = static_cast<size_t>(m) * static_cast<size_t>(m - 1) / 2;
        param_count_out = static_cast<size_t>(m) + d_count + static_cast<size_t>(m) + 1;
        effect_count_out = static_cast<size_t>(4 + m);
        param_vars.assign(param_count_out, OnlineVariance{});
        effect_vars.assign(effect_count_out, OnlineVariance{});
    }
};

inline bool build_dynamic_sigma(const DynamicPathCoefs& coefs,
                               double var_x,
                               const std::vector<double>& var_e,
                               int m,
                               WorkerBuffers& buf,
                               MatrixXd* out_sigma) {
    const int p = m + 2;
    if (p <= 0) {
        return false;
    }
    if (!(var_x > 0.0) || !std::isfinite(var_x)) {
        return false;
    }
    if (static_cast<int>(var_e.size()) != m + 1) {
        return false;
    }
    for (double v : var_e) {
        if (!(v > 0.0) || !std::isfinite(v)) {
            return false;
        }
    }

    buf.A.setZero(p, p);
    // ordering: X, M1..Mm, Y
    for (int j = 1; j <= m; ++j) {
        buf.A(j, 0) = coefs.a[static_cast<size_t>(j - 1)];
        for (int k = 1; k < j; ++k) {
            buf.A(j, k) = coefs.d_jk(j, k);
        }
    }
    const int y_idx = p - 1;
    buf.A(y_idx, 0) = coefs.cp;
    for (int j = 1; j <= m; ++j) {
        buf.A(y_idx, j) = coefs.b[static_cast<size_t>(j - 1)];
    }

    buf.M.setIdentity(p, p);
    buf.M.noalias() -= buf.A;

    buf.invM.setIdentity(p, p);
    buf.M.template triangularView<Eigen::Lower>().solveInPlace(buf.invM);
    if (!buf.invM.allFinite()) {
        return false;
    }

    buf.Omega.setZero(p, p);
    buf.Omega(0, 0) = var_x;
    for (int j = 1; j <= m; ++j) {
        buf.Omega(j, j) = var_e[static_cast<size_t>(j - 1)];
    }
    buf.Omega(y_idx, y_idx) = var_e[static_cast<size_t>(m)];

    *out_sigma = buf.invM * buf.Omega * buf.invM.transpose();
    return out_sigma->allFinite();
}

inline bool fit_dynamic_path_model(const Eigen::Ref<const MatrixXd>& S,
                                  int n_complete,
                                  int m,
                                  PathModelType model,
                                  WorkerBuffers& buf,
                                  DynamicPathFit* out) {
    out->ok = false;
    out->n_complete = n_complete;
    out->m = m;
    out->df = model_df(model, m);

    const int p = m + 2;
    if (S.rows() != p || S.cols() != p) {
        return false;
    }
    if (out->df == NA_INTEGER) {
        return false;
    }

    const double var_x = S(0, 0);
    if (!(var_x > 0.0) || !std::isfinite(var_x)) {
        return false;
    }

    out->coefs.resize(m);
    out->var_x = var_x;
    out->var_e.assign(static_cast<size_t>(m + 1), NA_REAL);
    out->sigma.resize(p, p);

    // Fit mediator regressions.
    for (int j = 1; j <= m; ++j) {
        const double Syy = S(j, j);
        if (!(Syy > 0.0) || !std::isfinite(Syy)) {
            return false;
        }

        int q = 0;
        if (model == PathModelType::Full && j > 1) {
            // Mj ~ M(j-1)
            q = 1;
            buf.Spp(0, 0) = S(j - 1, j - 1);
            buf.Sxy(0) = S(j - 1, j);
            double resid = NA_REAL;
            if (!regression_from_cov(buf.Spp.topLeftCorner(q, q),
                                     buf.Sxy.head(q),
                                     Syy,
                                     buf.beta.head(q),
                                     &resid)) {
                return false;
            }
            out->coefs.d_jk(j, j - 1) = buf.beta(0);
            out->var_e[static_cast<size_t>(j - 1)] = resid;
            continue;
        }

        // Partial/Semifull (and Full for j==1): Mj ~ X + M1..M(j-1)
        q = 1 + (j - 1);
        // predictors indices: 0, 1, ..., j-1
        for (int a = 0; a < q; ++a) {
            for (int b = 0; b < q; ++b) {
                buf.Spp(a, b) = S(a, b);
            }
            buf.Sxy(a) = S(a, j);
        }

        double resid = NA_REAL;
        if (!regression_from_cov(buf.Spp.topLeftCorner(q, q),
                                 buf.Sxy.head(q),
                                 Syy,
                                 buf.beta.head(q),
                                 &resid)) {
            return false;
        }

        out->coefs.a[static_cast<size_t>(j - 1)] = buf.beta(0);
        for (int k = 1; k < j; ++k) {
            out->coefs.d_jk(j, k) = buf.beta(k);
        }
        out->var_e[static_cast<size_t>(j - 1)] = resid;
    }

    // Fit outcome regression.
    const int y_idx = p - 1;
    const double var_y = S(y_idx, y_idx);
    if (!(var_y > 0.0) || !std::isfinite(var_y)) {
        return false;
    }

    if (model == PathModelType::Full) {
        // Y ~ Mm
        const int mm_idx = m;
        buf.Spp(0, 0) = S(mm_idx, mm_idx);
        buf.Sxy(0) = S(mm_idx, y_idx);
        double resid = NA_REAL;
        if (!regression_from_cov(buf.Spp.topLeftCorner(1, 1),
                                 buf.Sxy.head(1),
                                 var_y,
                                 buf.beta.head(1),
                                 &resid)) {
            return false;
        }
        out->coefs.b[static_cast<size_t>(m - 1)] = buf.beta(0);
        out->var_e[static_cast<size_t>(m)] = resid;
    } else if (model == PathModelType::Semifull) {
        // Y ~ M1..Mm (no X)
        const int q = m;
        for (int a = 0; a < q; ++a) {
            for (int b = 0; b < q; ++b) {
                buf.Spp(a, b) = S(a + 1, b + 1);
            }
            buf.Sxy(a) = S(a + 1, y_idx);
        }
        double resid = NA_REAL;
        if (!regression_from_cov(buf.Spp.topLeftCorner(q, q),
                                 buf.Sxy.head(q),
                                 var_y,
                                 buf.beta.head(q),
                                 &resid)) {
            return false;
        }
        for (int j = 1; j <= m; ++j) {
            out->coefs.b[static_cast<size_t>(j - 1)] = buf.beta(j - 1);
        }
        out->coefs.cp = 0.0;
        out->var_e[static_cast<size_t>(m)] = resid;
    } else {
        // Partial: Y ~ X + M1..Mm
        const int q = m + 1;
        for (int a = 0; a < q; ++a) {
            const int ia = a;
            for (int b = 0; b < q; ++b) {
                const int ib = b;
                buf.Spp(a, b) = S(ia, ib);
            }
            buf.Sxy(a) = S(ia, y_idx);
        }
        double resid = NA_REAL;
        if (!regression_from_cov(buf.Spp.topLeftCorner(q, q),
                                 buf.Sxy.head(q),
                                 var_y,
                                 buf.beta.head(q),
                                 &resid)) {
            return false;
        }
        out->coefs.cp = buf.beta(0);
        for (int j = 1; j <= m; ++j) {
            out->coefs.b[static_cast<size_t>(j - 1)] = buf.beta(j);
        }
        out->var_e[static_cast<size_t>(m)] = resid;
    }

    if (!build_dynamic_sigma(out->coefs, out->var_x, out->var_e, m, buf, &out->sigma)) {
        return false;
    }

    out->ok = true;
    return true;
}

inline bool compute_baseline_chisq(const Eigen::Ref<const MatrixXd>& S,
                                  int n_complete,
                                  double* out_chisq) {
    if (n_complete <= 1) {
        return false;
    }
    double logdet_S = NA_REAL;
    if (!logdet_spd(S, &logdet_S)) {
        return false;
    }

    double logdet_diag = 0.0;
    for (int i = 0; i < S.rows(); ++i) {
        const double v = S(i, i);
        if (!(v > 0.0) || !std::isfinite(v)) {
            return false;
        }
        logdet_diag += std::log(v);
    }

    const double f0 = logdet_diag - logdet_S;
    double chisq0 = static_cast<double>(n_complete - 1) * f0;
    if (!std::isfinite(chisq0)) {
        return false;
    }
    if (chisq0 < 0.0) {
        chisq0 = 0.0;
    }
    *out_chisq = chisq0;
    return true;
}

inline bool compute_dynamic_fit_measures(const Eigen::Ref<const MatrixXd>& S,
                                        const Eigen::Ref<const MatrixXd>& sigma,
                                        int n_complete,
                                        PathModelType model,
                                        int m,
                                        double chisq_baseline,
                                        int df_baseline,
                                        DynamicPathFit* out) {
    const int p = m + 2;
    const int df = model_df(model, m);
    const int param_count = model_param_count(model, m);
    if (df == NA_INTEGER || param_count == NA_INTEGER || df < 0 || n_complete <= 1) {
        return false;
    }
    if (S.rows() != p || S.cols() != p || sigma.rows() != p || sigma.cols() != p) {
        return false;
    }

    double logdet_S = NA_REAL;
    double logdet_sigma = NA_REAL;
    if (!logdet_spd(S, &logdet_S)) {
        return false;
    }
    if (!logdet_spd(sigma, &logdet_sigma)) {
        return false;
    }

    MatrixXd solved;
    if (!solve_spd(sigma, S, &solved)) {
        return false;
    }
    const double trace_term = solved.trace();
    if (!std::isfinite(trace_term)) {
        return false;
    }

    const double f_ml = logdet_sigma + trace_term - logdet_S - static_cast<double>(p);
    double chisq = static_cast<double>(n_complete - 1) * f_ml;
    if (!std::isfinite(chisq)) {
        return false;
    }
    if (chisq < 0.0) {
        chisq = 0.0;
    }

    const double pvalue = (df == 0 ? 1.0 : chisq_pvalue_upper(chisq, df));

    double rmsea = 0.0;
    if (df > 0) {
        const double denom = static_cast<double>(df) * static_cast<double>(n_complete - 1);
        const double ratio = (chisq - static_cast<double>(df)) / denom;
        rmsea = safe_sqrt(std::max(0.0, ratio));
    }

    double cfi = NA_REAL;
    double tli = NA_REAL;
    if (df == 0) {
        cfi = 1.0;
        tli = 1.0;
    } else {
        const double num = std::max(chisq - static_cast<double>(df), 0.0);
        const double base = std::max(chisq_baseline - static_cast<double>(df_baseline), 0.0);
        const double den = std::max(base, num);
        if (den <= 0.0) {
            cfi = 1.0;
        } else {
            cfi = 1.0 - num / den;
        }

        const double t0 = chisq_baseline / static_cast<double>(df_baseline);
        const double t1 = chisq / static_cast<double>(df);
        const double denom = t0 - 1.0;
        if (std::fabs(denom) < 1e-12) {
            tli = 1.0;
        } else {
            tli = (t0 - t1) / denom;
        }
    }

    const double neg2loglik =
        static_cast<double>(n_complete) *
        (static_cast<double>(p) * std::log(2.0 * std::acos(-1.0)) + logdet_sigma + trace_term);
    const double aic = neg2loglik + 2.0 * static_cast<double>(param_count);

    // df==0 convention (hard contract): chisq=0, p=1, cfi=1, tli=1, rmsea=0
    out->chisq = (df == 0 ? 0.0 : chisq);
    out->df = df;
    out->pvalue = (df == 0 ? 1.0 : pvalue);
    out->cfi = (df == 0 ? 1.0 : cfi);
    out->tli = (df == 0 ? 1.0 : tli);
    out->rmsea = (df == 0 ? 0.0 : rmsea);
    out->aic = aic;

    return true;
}

inline bool draw_wishart_dynamic(std::mt19937_64& rng,
                                const Eigen::Ref<const MatrixXd>& sigma,
                                int df,
                                WorkerBuffers& buf,
                                MatrixXd* out_S) {
    const int p = sigma.rows();
    if (sigma.cols() != p) {
        return false;
    }
    if (df < p) {
        return false;
    }

    Eigen::LLT<MatrixXd> llt(sigma);
    if (llt.info() != Eigen::Success) {
        return false;
    }
    buf.L = llt.matrixL();

    std::normal_distribution<double> norm(0.0, 1.0);
    buf.W_A.setZero(p, p);

    for (int i = 0; i < p; ++i) {
        const int dof = df - i;
        if (dof <= 0) {
            return false;
        }
        std::chi_squared_distribution<double> chisq(static_cast<double>(dof));
        const double diag = std::sqrt(chisq(rng));
        if (!std::isfinite(diag) || diag <= 0.0) {
            return false;
        }
        buf.W_A(i, i) = diag;
        for (int j = 0; j < i; ++j) {
            buf.W_A(i, j) = norm(rng);
        }
    }

    const MatrixXd B = buf.L * buf.W_A;
    MatrixXd W = B * B.transpose();
    if (!W.allFinite()) {
        return false;
    }
    W /= static_cast<double>(df);
    *out_S = std::move(W);
    return out_S->allFinite();
}

inline void decode_combination(uint64_t idx,
                              const std::vector<uint64_t>& dims,
                              std::vector<uint64_t>& out) {
    out.resize(dims.size());
    for (int d = static_cast<int>(dims.size()) - 1; d >= 0; --d) {
        const uint64_t dim = dims[static_cast<size_t>(d)];
        out[static_cast<size_t>(d)] = idx % dim;
        idx /= dim;
    }
}

inline std::string escape_combination_token(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == '%') {
            out += "%25";
        } else if (c == '|') {
            out += "%7C";
        } else if (c == '=') {
            out += "%3D";
        } else {
            out.push_back(c);
        }
    }
    return out;
}

inline std::string build_combination_id(const std::string& x_name,
                                       const std::vector<std::string>& m_names,
                                       const std::string& y_name) {
    std::string out;
    out.reserve(x_name.size() + y_name.size() + 16 + 8 * m_names.size());
    out += "X=";
    out += escape_combination_token(x_name);
    for (size_t j = 0; j < m_names.size(); ++j) {
        out.push_back('|');
        out += "M";
        out += std::to_string(j + 1);
        out.push_back('=');
        out += escape_combination_token(m_names[j]);
    }
    out += "|Y=";
    out += escape_combination_token(y_name);
    return out;
}

inline std::string build_m_vars_field(const std::vector<std::string>& m_names) {
    std::string out;
    for (size_t j = 0; j < m_names.size(); ++j) {
        if (j > 0) {
            out.push_back(';');
        }
        out += "M";
        out += std::to_string(j + 1);
        out.push_back('=');
        out += m_names[j];
    }
    return out;
}

inline void append_csv_int(std::string& out, int value) {
    if (value == NA_INTEGER) {
        out += "NA";
        return;
    }
    out += std::to_string(value);
}

inline void append_est_block(std::string& out,
                            double est,
                            double se,
                            double sd_pred,
                            double sd_out) {
    out.push_back(',');
    append_csv_double_fixed6(out, est);
    out.push_back(',');
    append_csv_double_fixed6(out, se);
    out.push_back(',');

    double z = NA_REAL;
    double p = NA_REAL;
    if (std::isfinite(se)) {
        if (se > 0.0) {
            z = est / se;
            p = normal_pvalue_two_sided(z);
        } else {
            z = 0.0;
            p = 1.0;
        }
    }
    append_csv_double_fixed6(out, z);
    out.push_back(',');
    append_csv_double_fixed6(out, p);
    out.push_back(',');

    const double scale =
        (std::isfinite(sd_pred) && std::isfinite(sd_out) && sd_out > 0.0 ? (sd_pred / sd_out)
                                                                         : NA_REAL);
    const double est_std = (std::isfinite(scale) ? est * scale : NA_REAL);
    append_csv_double_fixed6(out, est_std);
}

inline bool parse_whole_number_index(double x, const char* name, uint64_t* out) {
    if (!std::isfinite(x)) {
        return false;
    }
    if (x < 0.0 || x > kMaxExactIntDouble) {
        return false;
    }
    if (std::floor(x) != x) {
        return false;
    }
    *out = static_cast<uint64_t>(x);
    return true;
}

inline uint64_t checked_mul_u64(uint64_t a, uint64_t b, const char* context) {
    if (a == 0 || b == 0) {
        return 0;
    }
    if (a > std::numeric_limits<uint64_t>::max() / b) {
        throw std::overflow_error(std::string("Combination count overflow: ") + context);
    }
    return a * b;
}

struct SerialPathWorker : public RcppParallel::Worker {
    const Eigen::Map<const Eigen::MatrixXd>& data;
    const std::vector<std::string>& column_names;
    const std::vector<int>& x_cols;
    const std::vector<std::vector<int>>& mediator_cols;  // length m
    const std::vector<int>& y_cols;
    const std::vector<uint64_t>& dims;  // [nx, nm1..nmm, ny]
    const int m;
    const int nrep;
    const uint64_t base_seed;
    const double alpha;
    const bool excel_safe_csv;
    const uint64_t first_global_idx;
    const uint64_t shard_count;
    const uint64_t chunk_offset;
    const bool write_fit;
    const bool write_params;
    const bool write_effects;
    const bool include_failure_reason;

    std::vector<std::string>& fit_lines;
    std::vector<std::string>& param_lines;
    std::vector<std::string>& effect_lines;

    SerialPathWorker(const Eigen::Map<const Eigen::MatrixXd>& data_,
                     const std::vector<std::string>& column_names_,
                     const std::vector<int>& x_cols_,
                     const std::vector<std::vector<int>>& mediator_cols_,
                     const std::vector<int>& y_cols_,
                     const std::vector<uint64_t>& dims_,
                     int m_,
                     int nrep_,
                     uint64_t base_seed_,
                     double alpha_,
                     bool excel_safe_csv_,
                     uint64_t first_global_idx_,
                     uint64_t shard_count_,
                     uint64_t chunk_offset_,
                     bool write_fit_,
                     bool write_params_,
                     bool write_effects_,
                     bool include_failure_reason_,
                     std::vector<std::string>& fit_lines_,
                     std::vector<std::string>& param_lines_,
                     std::vector<std::string>& effect_lines_)
        : data(data_),
          column_names(column_names_),
          x_cols(x_cols_),
          mediator_cols(mediator_cols_),
          y_cols(y_cols_),
          dims(dims_),
          m(m_),
          nrep(nrep_),
          base_seed(base_seed_),
          alpha(alpha_),
          excel_safe_csv(excel_safe_csv_),
          first_global_idx(first_global_idx_),
          shard_count(shard_count_),
          chunk_offset(chunk_offset_),
          write_fit(write_fit_),
          write_params(write_params_),
          write_effects(write_effects_),
          include_failure_reason(include_failure_reason_),
          fit_lines(fit_lines_),
          param_lines(param_lines_),
          effect_lines(effect_lines_) {}

    void operator()(std::size_t begin, std::size_t end) {
        WorkerBuffers buf;
        buf.resize(m);

        for (std::size_t i = begin; i < end; ++i) {
            const uint64_t local_idx = static_cast<uint64_t>(i);
            const uint64_t global_idx = first_global_idx + (chunk_offset + local_idx) * shard_count;
            process_one(global_idx, local_idx, buf);
        }
    }

    void process_one(uint64_t global_idx, uint64_t local_idx, WorkerBuffers& buf) const {
        const int p = m + 2;
        const int y_pos = p - 1;

        decode_combination(global_idx, dims, buf.decoded);
        const size_t x_list = static_cast<size_t>(buf.decoded[0]);
        const size_t y_list = static_cast<size_t>(buf.decoded.back());
        const int x_idx = x_cols[x_list];
        const int y_idx = y_cols[y_list];

        std::vector<int>& cols = buf.var_cols;
        cols[0] = x_idx;
        cols[y_pos] = y_idx;
        std::vector<std::string>& m_names = buf.m_names;
        m_names.clear();
        m_names.reserve(static_cast<size_t>(m));
        for (int j = 0; j < m; ++j) {
            const size_t mj_list = static_cast<size_t>(buf.decoded[static_cast<size_t>(j + 1)]);
            const int mj_idx = mediator_cols[static_cast<size_t>(j)][mj_list];
            cols[j + 1] = mj_idx;
            m_names.push_back(column_names[mj_idx]);
        }

        const std::string& x_name = column_names[x_idx];
        const std::string& y_name = column_names[y_idx];
        const std::string combination_id = build_combination_id(x_name, m_names, y_name);
        const std::string m_vars = build_m_vars_field(m_names);

        auto write_failed_fit_row = [&](int n_complete, const std::string& failure_reason) {
            if (!write_fit) {
                return;
            }
            std::string out;
            out.reserve(combination_id.size() + m_vars.size() + 256);
            out += csv_escape(combination_id, excel_safe_csv);
            out.push_back(',');
            out += csv_escape(x_name, excel_safe_csv);
            out.push_back(',');
            out += csv_escape(y_name, excel_safe_csv);
            out.push_back(',');
            out += csv_escape(m_vars, excel_safe_csv);
            out.push_back(',');
            append_csv_int(out, m);
            out.push_back(',');
            append_csv_int(out, n_complete);
            out.push_back(',');
            out += "NA";  // selected_model
            if (include_failure_reason) {
                out.push_back(',');
                out += csv_escape(failure_reason, excel_safe_csv);
            }
            // Remaining columns: p-values and fit measures
            const int na_cols = 2 /*pvals*/ + 3 * 7 /*fit blocks*/;
            for (int k = 0; k < na_cols; ++k) {
                out += ",NA";
            }
            out.push_back('\n');
            fit_lines[static_cast<size_t>(local_idx)] = std::move(out);
        };

        // Detect duplicate columns: failed combination (no params/effects).
        {
            bool dup = false;
            for (int a = 0; a < p && !dup; ++a) {
                for (int b = a + 1; b < p; ++b) {
                    if (cols[a] == cols[b]) {
                        dup = true;
                        break;
                    }
                }
            }
            if (dup) {
                write_failed_fit_row(NA_INTEGER, "duplicate_columns");
                return;
            }
        }

        const int n_total = static_cast<int>(data.rows());
        if (n_total <= 0) {
            write_failed_fit_row(NA_INTEGER, "empty_data");
            return;
        }

        // Compute complete-case means.
        buf.sums.setZero(p);
        int n_complete = 0;
        for (int i = 0; i < n_total; ++i) {
            bool ok = true;
            for (int j = 0; j < p; ++j) {
                const double v = data(i, cols[j]);
                if (!std::isfinite(v)) {
                    ok = false;
                    break;
                }
                buf.centered(j) = v;
            }
            if (!ok) {
                continue;
            }
            buf.sums.noalias() += buf.centered;
            ++n_complete;
        }

        if (n_complete < p + 1) {
            write_failed_fit_row(n_complete, "insufficient_complete_cases");
            return;
        }

        buf.means = buf.sums / static_cast<double>(n_complete);

        // Compute covariance matrix with denominator (n_complete - 1).
        buf.S.setZero(p, p);
        for (int i = 0; i < n_total; ++i) {
            bool ok = true;
            for (int j = 0; j < p; ++j) {
                const double v = data(i, cols[j]);
                if (!std::isfinite(v)) {
                    ok = false;
                    break;
                }
                buf.centered(j) = v - buf.means(j);
            }
            if (!ok) {
                continue;
            }
            for (int a = 0; a < p; ++a) {
                for (int b = a; b < p; ++b) {
                    buf.S(a, b) += buf.centered(a) * buf.centered(b);
                }
            }
        }
        const double inv_df = 1.0 / static_cast<double>(n_complete - 1);
        for (int a = 0; a < p; ++a) {
            for (int b = a; b < p; ++b) {
                const double v = buf.S(a, b) * inv_df;
                buf.S(a, b) = v;
                buf.S(b, a) = v;
            }
        }

        // Standard deviations for standardization.
        for (int j = 0; j < p; ++j) {
            buf.sd(j) = safe_sqrt(buf.S(j, j));
        }

        double chisq_baseline = NA_REAL;
        const int df_baseline = (p * (p - 1)) / 2;
        if (!compute_baseline_chisq(buf.S, n_complete, &chisq_baseline)) {
            write_failed_fit_row(n_complete, "baseline_fit_failed");
            return;
        }

        DynamicPathFit fit_partial;
        DynamicPathFit fit_semifull;
        DynamicPathFit fit_full;
        if (!fit_dynamic_path_model(buf.S, n_complete, m, PathModelType::Partial, buf, &fit_partial)) {
            write_failed_fit_row(n_complete, "partial_model_fit_failed");
            return;
        }
        if (!fit_dynamic_path_model(buf.S, n_complete, m, PathModelType::Semifull, buf, &fit_semifull)) {
            write_failed_fit_row(n_complete, "semifull_model_fit_failed");
            return;
        }
        if (!fit_dynamic_path_model(buf.S, n_complete, m, PathModelType::Full, buf, &fit_full)) {
            write_failed_fit_row(n_complete, "full_model_fit_failed");
            return;
        }

        if (!compute_dynamic_fit_measures(buf.S,
                                          fit_partial.sigma,
                                          n_complete,
                                          PathModelType::Partial,
                                          m,
                                          chisq_baseline,
                                          df_baseline,
                                          &fit_partial)) {
            write_failed_fit_row(n_complete, "partial_fit_measures_failed");
            return;
        }
        if (!compute_dynamic_fit_measures(buf.S,
                                          fit_semifull.sigma,
                                          n_complete,
                                          PathModelType::Semifull,
                                          m,
                                          chisq_baseline,
                                          df_baseline,
                                          &fit_semifull)) {
            write_failed_fit_row(n_complete, "semifull_fit_measures_failed");
            return;
        }
        if (!compute_dynamic_fit_measures(buf.S,
                                          fit_full.sigma,
                                          n_complete,
                                          PathModelType::Full,
                                          m,
                                          chisq_baseline,
                                          df_baseline,
                                          &fit_full)) {
            write_failed_fit_row(n_complete, "full_fit_measures_failed");
            return;
        }

        const double diff_p_sf = std::max(0.0, fit_semifull.chisq - fit_partial.chisq);
        const int diff_df_sf = fit_semifull.df - fit_partial.df;
        const double p_P_vs_SF =
            (diff_df_sf > 0 ? chisq_pvalue_upper(diff_p_sf, diff_df_sf) : NA_REAL);

        const double diff_p_f = std::max(0.0, fit_full.chisq - fit_semifull.chisq);
        const int diff_df_f = fit_full.df - fit_semifull.df;
        double p_SF_vs_F = NA_REAL;
        if (diff_df_f > 0) {
            p_SF_vs_F = chisq_pvalue_upper(diff_p_f, diff_df_f);
        } else if (diff_df_f == 0 && std::isfinite(diff_p_f)) {
            p_SF_vs_F = (diff_p_f == 0.0 ? 1.0 : 0.0);
        }

        PathModelType selected = PathModelType::Full;
        const char* selected_label = "Full";
        if (!std::isfinite(p_P_vs_SF) || p_P_vs_SF < alpha) {
            selected = PathModelType::Partial;
            selected_label = "Partial";
        } else if (!std::isfinite(p_SF_vs_F) || p_SF_vs_F < alpha) {
            selected = PathModelType::Semifull;
            selected_label = "Semifull";
        }

        const DynamicPathFit* selected_fit = nullptr;
        if (selected == PathModelType::Partial) {
            selected_fit = &fit_partial;
        } else if (selected == PathModelType::Semifull) {
            selected_fit = &fit_semifull;
        } else {
            selected_fit = &fit_full;
        }
        const DynamicPathCoefs& b = selected_fit->coefs;

        // Defined effects (O(m^2) recursion; no path enumeration).
        std::vector<double> te_x_m(static_cast<size_t>(m), 0.0);
        te_x_m[0] = b.a[0];
        for (int j = 2; j <= m; ++j) {
            double te = b.a[static_cast<size_t>(j - 1)];
            for (int k = 1; k < j; ++k) {
                te += b.d_jk(j, k) * te_x_m[static_cast<size_t>(k - 1)];
            }
            te_x_m[static_cast<size_t>(j - 1)] = te;
        }
        double total = b.cp;
        for (int j = 1; j <= m; ++j) {
            total += b.b[static_cast<size_t>(j - 1)] * te_x_m[static_cast<size_t>(j - 1)];
        }
        const double direct = b.cp;
        const double total_indirect = total - direct;
        double adjacent_chain = b.a[0];
        for (int j = 2; j <= m; ++j) {
            adjacent_chain *= b.d_jk(j, j - 1);
        }
        adjacent_chain *= b.b[static_cast<size_t>(m - 1)];

        // Bootstrap SEs (only if needed for params/effects).
        const auto param_index_a = [&](int j) -> size_t {
            return static_cast<size_t>(j - 1);
        };
        const auto param_index_d = [&](int j, int k) -> size_t {
            return static_cast<size_t>(m) + DynamicPathCoefs::d_index(j, k);
        };
        const auto param_index_b = [&](int j) -> size_t {
            return static_cast<size_t>(m) + buf.d_count + static_cast<size_t>(j - 1);
        };
        const size_t param_index_cp =
            static_cast<size_t>(m) + buf.d_count + static_cast<size_t>(m);

        const auto effect_index_direct = [&]() -> size_t { return 0; };
        const auto effect_index_total = [&]() -> size_t { return 1; };
        const auto effect_index_total_indirect = [&]() -> size_t { return 2; };
        const auto effect_index_adjacent_chain = [&]() -> size_t { return 3; };
        const auto effect_index_x_to_m_total = [&](int j) -> size_t {
            return static_cast<size_t>(4 + (j - 1));
        };

        const bool need_se = (nrep > 0) && (write_params || write_effects);
        if (need_se) {
            for (auto& v : buf.param_vars) {
                v = OnlineVariance{};
            }
            for (auto& v : buf.effect_vars) {
                v = OnlineVariance{};
            }

            std::mt19937_64 rng(derive_seed(base_seed, global_idx, 0));
            const int wishart_df = n_complete - 1;

            for (int rep = 0; rep < nrep; ++rep) {
                MatrixXd Srep;
                if (!draw_wishart_dynamic(rng, selected_fit->sigma, wishart_df, buf, &Srep)) {
                    continue;
                }
                DynamicPathFit rep_fit;
                if (!fit_dynamic_path_model(Srep, n_complete, m, selected, buf, &rep_fit)) {
                    continue;
                }
                const DynamicPathCoefs& rb = rep_fit.coefs;

                // params
                for (int j = 1; j <= m; ++j) {
                    buf.param_vars[param_index_a(j)].add(rb.a[static_cast<size_t>(j - 1)]);
                }
                for (int j = 2; j <= m; ++j) {
                    for (int k = 1; k < j; ++k) {
                        buf.param_vars[param_index_d(j, k)].add(rb.d_jk(j, k));
                    }
                }
                for (int j = 1; j <= m; ++j) {
                    buf.param_vars[param_index_b(j)].add(rb.b[static_cast<size_t>(j - 1)]);
                }
                buf.param_vars[param_index_cp].add(rb.cp);

                // effects
                std::vector<double> r_te_x_m(static_cast<size_t>(m), 0.0);
                r_te_x_m[0] = rb.a[0];
                for (int j = 2; j <= m; ++j) {
                    double te = rb.a[static_cast<size_t>(j - 1)];
                    for (int k = 1; k < j; ++k) {
                        te += rb.d_jk(j, k) * r_te_x_m[static_cast<size_t>(k - 1)];
                    }
                    r_te_x_m[static_cast<size_t>(j - 1)] = te;
                }
                double r_total = rb.cp;
                for (int j = 1; j <= m; ++j) {
                    r_total += rb.b[static_cast<size_t>(j - 1)] * r_te_x_m[static_cast<size_t>(j - 1)];
                }
                const double r_direct = rb.cp;
                const double r_total_ind = r_total - r_direct;
                double r_adj_chain = rb.a[0];
                for (int j = 2; j <= m; ++j) {
                    r_adj_chain *= rb.d_jk(j, j - 1);
                }
                r_adj_chain *= rb.b[static_cast<size_t>(m - 1)];

                buf.effect_vars[effect_index_direct()].add(r_direct);
                buf.effect_vars[effect_index_total()].add(r_total);
                buf.effect_vars[effect_index_total_indirect()].add(r_total_ind);
                buf.effect_vars[effect_index_adjacent_chain()].add(r_adj_chain);
                for (int j = 1; j <= m; ++j) {
                    buf.effect_vars[effect_index_x_to_m_total(j)].add(
                        r_te_x_m[static_cast<size_t>(j - 1)]);
                }
            }
        }

        const auto se_or_na = [&](const OnlineVariance& v) -> double {
            if (!need_se) {
                return NA_REAL;
            }
            if (v.n < 2) {
                return NA_REAL;
            }
            return safe_sqrt(v.variance());
        };

        // Write fit row (1 per combination).
        if (write_fit) {
            std::string out;
            out.reserve(combination_id.size() + m_vars.size() + 512);
            out += csv_escape(combination_id, excel_safe_csv);
            out.push_back(',');
            out += csv_escape(x_name, excel_safe_csv);
            out.push_back(',');
            out += csv_escape(y_name, excel_safe_csv);
            out.push_back(',');
            out += csv_escape(m_vars, excel_safe_csv);
            out.push_back(',');
            append_csv_int(out, m);
            out.push_back(',');
            append_csv_int(out, n_complete);
            out.push_back(',');
            out += csv_escape(selected_label, excel_safe_csv);
            if (include_failure_reason) {
                out += ",NA";
            }
            out.push_back(',');
            append_csv_double_fixed6(out, p_P_vs_SF);
            out.push_back(',');
            append_csv_double_fixed6(out, p_SF_vs_F);
            out.push_back(',');

            auto append_fit_block = [&](const DynamicPathFit& fit) {
                append_csv_double_fixed6(out, fit.chisq);
                out.push_back(',');
                append_csv_int(out, fit.df);
                out.push_back(',');
                append_csv_double_fixed6(out, fit.pvalue);
                out.push_back(',');
                append_csv_double_fixed6(out, fit.cfi);
                out.push_back(',');
                append_csv_double_fixed6(out, fit.tli);
                out.push_back(',');
                append_csv_double_fixed6(out, fit.rmsea);
                out.push_back(',');
                append_csv_double_fixed6(out, fit.aic);
            };

            append_fit_block(fit_partial);
            out.push_back(',');
            append_fit_block(fit_semifull);
            out.push_back(',');
            append_fit_block(fit_full);
            out.push_back('\n');
            fit_lines[static_cast<size_t>(local_idx)] = std::move(out);
        }

        // Write parameter rows (selected model only; long format).
        if (write_params) {
            std::string out;
            out.reserve(256 + static_cast<size_t>(m) * 128);

            for (int j = 1; j <= m; ++j) {
                const double est = b.a[static_cast<size_t>(j - 1)];
                const double se = se_or_na(buf.param_vars[param_index_a(j)]);
                out += csv_escape(combination_id, excel_safe_csv);
                out.push_back(',');
                out += csv_escape(selected_label, excel_safe_csv);
                out.push_back(',');
                out += csv_escape(m_names[static_cast<size_t>(j - 1)], excel_safe_csv);
                out.push_back(',');
                out += "~,";
                out += csv_escape(x_name, excel_safe_csv);
                out.push_back(',');
                out += csv_escape("a" + std::to_string(j), excel_safe_csv);
                append_est_block(out, est, se, buf.sd(0), buf.sd(j));
                out.push_back('\n');

                for (int k = 1; k < j; ++k) {
                    const double est_d = b.d_jk(j, k);
                    const double se_d = se_or_na(buf.param_vars[param_index_d(j, k)]);
                    out += csv_escape(combination_id, excel_safe_csv);
                    out.push_back(',');
                    out += csv_escape(selected_label, excel_safe_csv);
                    out.push_back(',');
                    out += csv_escape(m_names[static_cast<size_t>(j - 1)], excel_safe_csv);
                    out.push_back(',');
                    out += "~,";
                    out += csv_escape(m_names[static_cast<size_t>(k - 1)], excel_safe_csv);
                    out.push_back(',');
                    out += csv_escape("d_" + std::to_string(j) + "_" + std::to_string(k), excel_safe_csv);
                    append_est_block(out, est_d, se_d, buf.sd(k), buf.sd(j));
                    out.push_back('\n');
                }
            }

            // Y ~ X (cp)
            {
                const double est = b.cp;
                const double se = se_or_na(buf.param_vars[param_index_cp]);
                out += csv_escape(combination_id, excel_safe_csv);
                out.push_back(',');
                out += csv_escape(selected_label, excel_safe_csv);
                out.push_back(',');
                out += csv_escape(y_name, excel_safe_csv);
                out.push_back(',');
                out += "~,";
                out += csv_escape(x_name, excel_safe_csv);
                out.push_back(',');
                out += "cp";
                append_est_block(out, est, se, buf.sd(0), buf.sd(y_pos));
                out.push_back('\n');
            }

            for (int j = 1; j <= m; ++j) {
                const double est = b.b[static_cast<size_t>(j - 1)];
                const double se = se_or_na(buf.param_vars[param_index_b(j)]);
                out += csv_escape(combination_id, excel_safe_csv);
                out.push_back(',');
                out += csv_escape(selected_label, excel_safe_csv);
                out.push_back(',');
                out += csv_escape(y_name, excel_safe_csv);
                out.push_back(',');
                out += "~,";
                out += csv_escape(m_names[static_cast<size_t>(j - 1)], excel_safe_csv);
                out.push_back(',');
                out += csv_escape("b" + std::to_string(j), excel_safe_csv);
                append_est_block(out, est, se, buf.sd(j), buf.sd(y_pos));
                out.push_back('\n');
            }

            param_lines[static_cast<size_t>(local_idx)] = std::move(out);
        }

        // Write effects rows (selected model only; long format).
        if (write_effects) {
            std::string out;
            out.reserve(256 + static_cast<size_t>(m) * 128);

            auto write_effect = [&](const std::string& name,
                                    double est,
                                    double se,
                                    double sd_out) {
                out += csv_escape(combination_id, excel_safe_csv);
                out.push_back(',');
                out += csv_escape(name, excel_safe_csv);
                append_est_block(out, est, se, buf.sd(0), sd_out);
                out.push_back('\n');
            };

            write_effect("direct",
                         direct,
                         se_or_na(buf.effect_vars[effect_index_direct()]),
                         buf.sd(y_pos));
            write_effect("total",
                         total,
                         se_or_na(buf.effect_vars[effect_index_total()]),
                         buf.sd(y_pos));
            write_effect("total_indirect",
                         total_indirect,
                         se_or_na(buf.effect_vars[effect_index_total_indirect()]),
                         buf.sd(y_pos));
            write_effect("adjacent_chain",
                         adjacent_chain,
                         se_or_na(buf.effect_vars[effect_index_adjacent_chain()]),
                         buf.sd(y_pos));

            for (int j = 1; j <= m; ++j) {
                write_effect("x_to_m" + std::to_string(j) + "_total",
                             te_x_m[static_cast<size_t>(j - 1)],
                             se_or_na(buf.effect_vars[effect_index_x_to_m_total(j)]),
                             buf.sd(j));
            }

            effect_lines[static_cast<size_t>(local_idx)] = std::move(out);
        }
    }
};

}  // namespace

// [[Rcpp::export]]
void serial_path_analysis_cpp(NumericMatrix data,
                              CharacterVector column_names,
                              IntegerVector x_col_idx,
                              List mediator_col_idx_list,
                              IntegerVector y_col_idx,
                              int m,
                              int nrep,
                              std::string output_fit_file,
                              std::string output_params_file,
                              std::string output_effects_file,
                              uint64_t base_seed = 0,
                              double alpha = 0.05,
                              double combination_start = 0.0,
                              double combination_end = NA_REAL,
                              uint64_t shard_id = 0,
                              uint64_t shard_count = 1,
                              int chunk_size = 1024,
                              int grain_size = 1,
                              bool overwrite = true,
                              bool excel_safe_csv = false,
                              bool write_fit = true,
                              bool write_params = true,
                              bool write_effects = true,
                              bool include_failure_reason = false) {
    if (data.nrow() == 0 || data.ncol() == 0) {
        throw std::invalid_argument("Data matrix is empty");
    }
    if (static_cast<int>(column_names.size()) != data.ncol()) {
        throw std::invalid_argument("column_names length does not match data.ncol()");
    }
    if (m < 1) {
        throw std::invalid_argument("m must be >= 1");
    }
    if (mediator_col_idx_list.size() != m) {
        throw std::invalid_argument("mediator_col_idx_list length must equal m");
    }
    if (nrep < 0) {
        throw std::invalid_argument("nrep must be >= 0");
    }
    if (!std::isfinite(alpha) || alpha <= 0.0 || alpha >= 1.0) {
        throw std::invalid_argument("alpha must be in (0,1)");
    }
    if (chunk_size <= 0) {
        throw std::invalid_argument("chunk_size must be positive");
    }
    if (grain_size <= 0) {
        throw std::invalid_argument("grain_size must be positive");
    }
    if (shard_count < 1) {
        throw std::invalid_argument("shard_count must be >= 1");
    }
    if (shard_id >= shard_count) {
        throw std::invalid_argument("shard_id must be in [0, shard_count)");
    }
    if (!write_fit && !write_params && !write_effects) {
        return;
    }
    if (write_fit && output_fit_file.empty()) {
        throw std::invalid_argument("output_fit_file must be non-empty when write_fit=TRUE");
    }
    if (write_params && output_params_file.empty()) {
        throw std::invalid_argument("output_params_file must be non-empty when write_params=TRUE");
    }
    if (write_effects && output_effects_file.empty()) {
        throw std::invalid_argument("output_effects_file must be non-empty when write_effects=TRUE");
    }

    std::vector<std::string> colnames_cpp = as<std::vector<std::string>>(column_names);

    std::vector<int> x_cols = as<std::vector<int>>(x_col_idx);
    std::vector<int> y_cols = as<std::vector<int>>(y_col_idx);
    std::vector<std::vector<int>> mediator_cols;
    mediator_cols.reserve(static_cast<size_t>(m));
    for (int i = 0; i < m; ++i) {
        IntegerVector stage = mediator_col_idx_list[i];
        mediator_cols.push_back(as<std::vector<int>>(stage));
    }

    auto validate_idx = [&](const std::vector<int>& idxs, const char* name) {
        if (idxs.empty()) {
            throw std::invalid_argument(std::string(name) + " must be non-empty");
        }
        for (int v : idxs) {
            if (v < 0 || v >= data.ncol()) {
                throw std::invalid_argument(std::string(name) + " out of bounds");
            }
        }
    };

    validate_idx(x_cols, "x_col_idx");
    for (int i = 0; i < m; ++i) {
        validate_idx(mediator_cols[static_cast<size_t>(i)], "mediator_col_idx_list");
    }
    validate_idx(y_cols, "y_col_idx");

    std::vector<uint64_t> dims;
    dims.reserve(static_cast<size_t>(m + 2));
    dims.push_back(static_cast<uint64_t>(x_cols.size()));
    for (int i = 0; i < m; ++i) {
        dims.push_back(static_cast<uint64_t>(mediator_cols[static_cast<size_t>(i)].size()));
    }
    dims.push_back(static_cast<uint64_t>(y_cols.size()));

    uint64_t total = 1;
    for (size_t i = 0; i < dims.size(); ++i) {
        total = checked_mul_u64(total, dims[i], "total combinations");
    }
    if (total == 0) {
        throw std::invalid_argument("No combinations to process");
    }

    uint64_t start_idx = 0;
    if (!parse_whole_number_index(combination_start, "combination_start", &start_idx)) {
        throw std::invalid_argument(
            "combination_start must be a whole number in [0, 2^53-1]");
    }
    uint64_t end_idx = total;
    if (!Rcpp::NumericVector::is_na(combination_end)) {
        if (!parse_whole_number_index(combination_end, "combination_end", &end_idx)) {
            throw std::invalid_argument(
                "combination_end must be NA or a whole number in [0, 2^53-1]");
        }
    }
    if (start_idx > end_idx) {
        throw std::invalid_argument("combination_start must be <= combination_end");
    }
    if (end_idx > total) {
        throw std::invalid_argument("combination_end exceeds total combinations");
    }

    // First global index in [start,end) with idx % shard_count == shard_id.
    uint64_t first_global_idx = 0;
    bool has_any = false;
    if (start_idx < end_idx) {
        const uint64_t start_mod = start_idx % shard_count;
        const uint64_t delta = (shard_id + shard_count - start_mod) % shard_count;
        first_global_idx = start_idx + delta;
        if (first_global_idx < end_idx) {
            has_any = true;
        }
    }

    uint64_t n_to_process = 0;
    if (has_any) {
        const uint64_t last = end_idx - 1;
        n_to_process = 1 + (last - first_global_idx) / shard_count;
    }

    // Open output files and write headers.
    auto open_for_write = [&](const std::string& path) -> std::ofstream {
        if (!overwrite) {
            std::ifstream existing(path.c_str());
            if (existing.good()) {
                throw std::invalid_argument("output file already exists: " + path);
            }
        }
        std::ofstream out_stream(path.c_str(), std::ios::out | std::ios::trunc);
        if (!out_stream.is_open()) {
            throw std::runtime_error("Unable to open output file: " + path);
        }
        return out_stream;
    };

    std::ofstream fit_stream;
    std::ofstream params_stream;
    std::ofstream effects_stream;

    if (write_fit) {
        fit_stream = open_for_write(output_fit_file);
        fit_stream <<
            "Combination,X_var,Y_var,M_vars,m,N_complete,selected_model,";
        if (include_failure_reason) {
            fit_stream << "failure_reason,";
        }
        fit_stream <<
            "p_P_vs_SF,p_SF_vs_F,"
            "partial_chisq,partial_df,partial_p,partial_cfi,partial_tli,partial_rmsea,partial_aic,"
            "semifull_chisq,semifull_df,semifull_p,semifull_cfi,semifull_tli,semifull_rmsea,semifull_aic,"
            "full_chisq,full_df,full_p,full_cfi,full_tli,full_rmsea,full_aic\n";
    }
    if (write_params) {
        params_stream = open_for_write(output_params_file);
        params_stream << "Combination,selected_model,lhs,op,rhs,label,est,se,z,p,std\n";
    }
    if (write_effects) {
        effects_stream = open_for_write(output_effects_file);
        effects_stream << "Combination,effect,est,se,z,p,std\n";
    }

    if (n_to_process == 0) {
        return;
    }

    Eigen::Map<const Eigen::MatrixXd> data_map(data.begin(), data.nrow(), data.ncol());
    const uint64_t chunk_size_u64 = static_cast<uint64_t>(chunk_size);
    const size_t grain_size_cpp = static_cast<size_t>(grain_size);

    for (uint64_t chunk_offset = 0; chunk_offset < n_to_process; chunk_offset += chunk_size_u64) {
        Rcpp::checkUserInterrupt();
        const uint64_t chunk_end = std::min(chunk_offset + chunk_size_u64, n_to_process);
        const size_t chunk_n = static_cast<size_t>(chunk_end - chunk_offset);

        std::vector<std::string> fit_lines(chunk_n);
        std::vector<std::string> param_lines(chunk_n);
        std::vector<std::string> effect_lines(chunk_n);

        SerialPathWorker worker(data_map,
                                colnames_cpp,
                                x_cols,
                                mediator_cols,
                                y_cols,
                                dims,
                                m,
                                nrep,
                                base_seed,
                                alpha,
                                excel_safe_csv,
                                first_global_idx,
                                shard_count,
                                chunk_offset,
                                write_fit,
                                write_params,
                                write_effects,
                                include_failure_reason,
                                fit_lines,
                                param_lines,
                                effect_lines);

        RcppParallel::parallelFor(static_cast<size_t>(0), chunk_n, worker, grain_size_cpp);

        for (size_t i = 0; i < chunk_n; ++i) {
            if (write_fit) {
                fit_stream << fit_lines[i];
            }
            if (write_params) {
                params_stream << param_lines[i];
            }
            if (write_effects) {
                effects_stream << effect_lines[i];
            }
        }
    }
}
