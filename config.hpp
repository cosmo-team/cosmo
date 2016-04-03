#pragma once
#ifndef CONFIG_HPP 
#define CONFIG_HPP
// This file contains type definitions, constants,
// And any other config stuff

#include <string>
#include <boost/tuple/tuple.hpp>
#include "utility.hpp"

using std::string;

// TODO: add support for wider integers
//#include <boost/multiprecision/cpp_int.hpp>
typedef __uint128_t uint128_t;
// need to define numeric limits since max() returns 0 on some compilers!!!!
namespace std {
  template <>
  struct numeric_limits<uint128_t> {
    static const size_t digits = CHAR_BIT*sizeof(uint128_t);

    static uint128_t min() {
      return 0;
    }

    static uint128_t max() {
      uint64_t max_64 = numeric_limits<uint64_t>::max();
      uint128_t temp = max_64;
      temp <<= 64;
      temp += max_64;
      return temp;
    }
  };
}

namespace cosmo {

// TYPES
// Select kmer type
static_assert(K_LEN <= 64 && K_LEN > 0,
  "k is currently limited to a positive integer <= 64.");

#if K_LEN <= 32
typedef uint64_t kmer_t;
#elif K_LEN <= 64
typedef uint128_t kmer_t;
#endif

typedef uint64_t color_t;
// For measuring the length of a kmer (LCS and dummies)
typedef uint8_t length_t;
typedef boost::tuple<kmer_t, length_t> dummy_t;

const string version = VERSION;
const string banner  = BANNER;
const size_t max_k = (K_LEN <= 32)? 32 : 64;
const size_t max_colors = bitwidth<color_t>::width;
const size_t mb_to_bytes = 1024 * 1024;
const size_t block_size       = 2 * 1024 * 1024; // KB
const size_t default_mem_size = 2 * 1024; // MB

// File extensions
const string graph_ext     = ".dbg";
const string bitvector_ext = ".bits";
const string packed_ext    = ".packed";
const string color_ext     = ".colors";
const string lcs_ext       = ".lcs";
}

#endif
