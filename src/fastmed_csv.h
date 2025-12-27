#pragma once

#include "fastmed_stats.h"

#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <string>
#include <system_error>
#include <type_traits>

inline std::string csv_escape(const std::string& field) {
    bool needs_escape = false;
    for (char c : field) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            needs_escape = true;
            break;
        }
    }
    if (!needs_escape) {
        return field;
    }

    std::string escaped = "\"";
    for (char c : field) {
        if (c == '"') {
            escaped += "\"\"";
        } else {
            escaped += c;
        }
    }
    escaped += "\"";
    return escaped;
}

inline std::string csv_escape(const std::string& field, bool excel_safe_csv) {
    if (!excel_safe_csv) {
        return csv_escape(field);
    }

    size_t i = 0;
    while (i < field.size() &&
           std::isspace(static_cast<unsigned char>(field[i])) != 0) {
        ++i;
    }
    if (i < field.size()) {
        const char c = field[i];
        if (c == '=' || c == '+' || c == '-' || c == '@') {
            std::string prefixed;
            prefixed.reserve(field.size() + 1);
            prefixed.push_back('\'');
            prefixed += field;
            return csv_escape(prefixed);
        }
    }

    return csv_escape(field);
}

namespace fastmed_detail {
template <typename T, typename = void>
struct has_to_chars_fixed : std::false_type {};

template <typename T>
struct has_to_chars_fixed<
    T,
    std::void_t<decltype(std::to_chars(static_cast<char*>(nullptr),
                                       static_cast<char*>(nullptr),
                                       std::declval<T>(),
                                       std::chars_format::fixed,
                                       6))>> : std::true_type {};
}  // namespace fastmed_detail

inline void append_csv_double_fixed6(std::string& out, double value) {
    if (!std::isfinite(value)) {
        out += "NA";
        return;
    }

    char buf[64];

    if constexpr (fastmed_detail::has_to_chars_fixed<double>::value) {
        auto [ptr, ec] =
            std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::fixed, 6);
        if (ec != std::errc()) {
            out += "NA";
            return;
        }
        out.append(buf, static_cast<size_t>(ptr - buf));
        return;
    }

    const int n = std::snprintf(buf, sizeof(buf), "%.6f", value);
    if (n <= 0 || static_cast<size_t>(n) >= sizeof(buf)) {
        out += "NA";
        return;
    }
    out.append(buf, static_cast<size_t>(n));
}

inline void append_statistics(std::string& out, const StatisticsSummary& stats) {
    append_csv_double_fixed6(out, stats.mean);
    out.push_back(',');
    append_csv_double_fixed6(out, stats.percentile_2_5);
    out.push_back(',');
    append_csv_double_fixed6(out, stats.percentile_97_5);
    out.push_back(',');
    append_csv_double_fixed6(out, stats.p_value);
}
