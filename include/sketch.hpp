#pragma once

#ifndef _SKETCH_H_
#define _SKETCH_H_

#include <array>
#include <string>
#include <cstdint>
#include <math.h>
#include "constants.hpp"

const uint64_t BOTTOM_BITS = log2(M); /* bits used to first hash */
const uint64_t MASK = M - 1;

const int NONE_HASH_ID = 0;
const int MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID = 1;
const int XXHASH64_ID = 2;

template <size_t K, int HASH_ID>
std::array<uint64_t, M> get_sketch_full_template( std::string input );

template <int HASH_ID>
std::array<uint64_t, M> get_sketch_hash_template( std::string input, uint8_t kmerlen );

std::array<uint64_t, M> get_sketch( std::string input, const uint8_t kmerlen, const int hash_id );

#endif