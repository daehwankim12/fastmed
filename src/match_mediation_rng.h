#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

enum class MatchMediationNoiseKind : uint8_t { None = 0, Normals = 1, Uniforms = 2 };

inline constexpr std::size_t kMatchMediationNoOffset =
    std::numeric_limits<std::size_t>::max();

struct MatchMediationRngSlice {
    std::size_t coef_offset = kMatchMediationNoOffset;  // doubles
    MatchMediationNoiseKind noise_kind = MatchMediationNoiseKind::None;
    std::size_t noise_offset = kMatchMediationNoOffset;  // doubles
    std::size_t idx_offset = kMatchMediationNoOffset;    // ints
};

struct MatchMediationRngBlock {
    std::size_t begin = 0;
    std::size_t end = 0;
    int n = 0;
    int nrep = 0;
    bool bootstrap = false;

    std::vector<MatchMediationRngSlice> slices;

    std::vector<double> coef_normals;
    std::vector<double> noise_normals;
    std::vector<double> noise_uniforms;
    std::vector<int> bootstrap_indices;
};

