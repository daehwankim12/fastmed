#pragma once

#include "glm_fit.h"

#include <cmath>

inline double inv_logit(double x) {
    if (x >= 0.0) {
        double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    double z = std::exp(x);
    return z / (1.0 + z);
}

inline double safe_exp(double x) {
    if (x > 700.0) {
        x = 700.0;
    }
    if (x < -700.0) {
        x = -700.0;
    }
    return std::exp(x);
}

inline double linkinv(double eta, GlmFamily fam) {
    if (fam == GlmFamily::Gaussian) {
        return eta;
    }
    if (fam == GlmFamily::Binomial) {
        return inv_logit(eta);
    }
    return safe_exp(eta);  // Poisson
}

