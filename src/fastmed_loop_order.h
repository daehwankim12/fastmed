#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

enum class LoopDimension : unsigned char { Exposure = 0, Mediator = 1, Outcome = 2 };

using LoopOrder = std::array<LoopDimension, 3>;  // outer, middle, inner

namespace {
inline std::string ascii_tolower_copy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

inline std::vector<std::string> split_tokens(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            if (!cur.empty()) {
                out.push_back(cur);
                cur.clear();
            }
            continue;
        }
        cur.push_back(c);
    }
    if (!cur.empty()) {
        out.push_back(cur);
    }
    return out;
}

inline LoopDimension parse_loop_dim_token(const std::string& token) {
    if (token == "e" || token == "exp" || token == "exposure") {
        return LoopDimension::Exposure;
    }
    if (token == "m" || token == "med" || token == "mediator") {
        return LoopDimension::Mediator;
    }
    if (token == "o" || token == "out" || token == "outcome") {
        return LoopDimension::Outcome;
    }
    throw std::invalid_argument("Invalid loop_order token: " + token);
}
}  // namespace

inline LoopOrder parse_loop_order(const std::string& loop_order) {
    std::string key = ascii_tolower_copy(loop_order);
    for (char& c : key) {
        if (c == '-' || c == ' ' || c == ',' || c == '/' || c == '.') {
            c = '_';
        }
    }

    if (key == "emo") {
        return {LoopDimension::Exposure, LoopDimension::Mediator, LoopDimension::Outcome};
    }
    if (key == "eom") {
        return {LoopDimension::Exposure, LoopDimension::Outcome, LoopDimension::Mediator};
    }
    if (key == "meo") {
        return {LoopDimension::Mediator, LoopDimension::Exposure, LoopDimension::Outcome};
    }
    if (key == "moe") {
        return {LoopDimension::Mediator, LoopDimension::Outcome, LoopDimension::Exposure};
    }
    if (key == "oem") {
        return {LoopDimension::Outcome, LoopDimension::Exposure, LoopDimension::Mediator};
    }
    if (key == "ome") {
        return {LoopDimension::Outcome, LoopDimension::Mediator, LoopDimension::Exposure};
    }

    const std::vector<std::string> tokens = split_tokens(key, '_');
    if (tokens.size() != 3) {
        throw std::invalid_argument(
            "loop_order must be a permutation of {exposure, mediator, outcome} "
            "(e.g. 'exposure_mediator_outcome' or 'EMO')");
    }

    LoopOrder order{
        parse_loop_dim_token(tokens[0]),
        parse_loop_dim_token(tokens[1]),
        parse_loop_dim_token(tokens[2]),
    };

    if (order[0] == order[1] || order[0] == order[2] || order[1] == order[2]) {
        throw std::invalid_argument(
            "loop_order must contain each of exposure, mediator, outcome exactly once");
    }

    return order;
}

inline void decode_combination_indices(std::size_t idx,
                                      std::size_t e,
                                      std::size_t m,
                                      std::size_t o,
                                      const LoopOrder& order,
                                      std::size_t& exp_list_idx,
                                      std::size_t& med_list_idx,
                                      std::size_t& out_list_idx) {
    const std::array<std::size_t, 3> sizes{e, m, o};
    std::array<std::size_t, 3> indices{0, 0, 0};

    const auto dim_index = [](LoopDimension dim) -> std::size_t {
        return static_cast<std::size_t>(dim);
    };

    const std::size_t outer_dim = dim_index(order[0]);
    const std::size_t middle_dim = dim_index(order[1]);
    const std::size_t inner_dim = dim_index(order[2]);

    const std::size_t inner_size = sizes[inner_dim];
    const std::size_t middle_size = sizes[middle_dim];
    if (inner_size == 0 || middle_size == 0 || sizes[outer_dim] == 0) {
        throw std::invalid_argument("decode_combination_indices: empty index list");
    }

    indices[inner_dim] = idx % inner_size;
    idx /= inner_size;
    indices[middle_dim] = idx % middle_size;
    idx /= middle_size;
    indices[outer_dim] = idx;

    exp_list_idx = indices[0];
    med_list_idx = indices[1];
    out_list_idx = indices[2];
}

inline uint64_t canonical_combination_index(std::size_t exp_list_idx,
                                           std::size_t med_list_idx,
                                           std::size_t out_list_idx,
                                           std::size_t m,
                                           std::size_t o) {
    return static_cast<uint64_t>((exp_list_idx * m + med_list_idx) * o + out_list_idx);
}
