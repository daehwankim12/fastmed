#pragma once

#include <cstdint>

inline uint64_t derive_seed(uint64_t base_seed,
                            uint64_t global_combination_idx,
                            uint64_t rep_idx) {
    uint64_t x = base_seed ^ (global_combination_idx << 32) ^ rep_idx;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

