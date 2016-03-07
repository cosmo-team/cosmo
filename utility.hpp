#pragma once
#ifndef UTILITY_HPP
#define UTILITY_HPP
// Utility functions

// For CHAR_BIT, which might not be 8 on certain systems
#include <climits>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
//#include <gmp.h>
#include <boost/random.hpp>
#include <sdsl/bit_vectors.hpp> // for size_in_bytes

#include "debug.hpp"

template <typename Container>
double bits_per_element(const Container & c) {
  return sdsl::size_in_bytes(c) * 8.0 / c.size();
}

namespace cosmo {

namespace fs = boost::filesystem;
using fs::extension;
using fs::basename;

// TODO: fill this in with comparators, shifters, whatever else to make it compile
// https://gmplib.org/manual/Low_002dlevel-Functions.html#Low_002dlevel-Functions
//template <size_t n_limbs>
//struct gmp_int {
//  mp_limb_t limbs[n_limbs];
//};

template <typename T>
struct bitwidth {
  const static size_t width = sizeof(T) * CHAR_BIT;
};

template <typename T>
T & deconst(const T & x) {
  return const_cast<T&>(x);
}

template <typename value_type, class IndexFunction>
ssize_t function_binary_search(size_t lo, size_t hi, value_type key, IndexFunction f) {
  while ( lo <= hi ) {
    size_t mid = lo + (hi - lo)/2;

    value_type temp = f(mid);
    if (temp < key) lo = mid + 1;
    else if (temp == key) return mid;
    else if (mid == 0) return -1;
    else hi = mid - 1;
  }
  return -1;
}

struct nt_to_int {
  uint64_t operator()(char x) const {
    switch(x) {
      case '$': return 0;
      case 'a': return 1;
      case 'c': return 2;
      case 'g': return 3;
      case 't': return 4;
      case 'A': return 5;
      case 'C': return 6;
      case 'G': return 7;
      case 'T': return 8;
      // Should never reach here
      default:
        COSMO_ASSERT(false);
        return -1;
    }
  }
};

struct int_to_nt {
  char operator()(size_t x) const {
    assert(x <= 8);
    return "$acgtACGT"[x];
  }
};

std::vector<size_t> random_uints(size_t lo, size_t hi, size_t n) {
  typedef boost::mt19937 rng_type;
  rng_type rng(time(0));

  std::vector<size_t> v(n);

  // Define random distributions
  boost::uniform_int<size_t> index_dist(lo, hi);
  boost::variate_generator<rng_type &, boost::uniform_int<size_t>> random_idx(rng, index_dist);

  for (size_t i = 0; i < n; ++i) {
    v[i] = random_idx();
  }
  return v;
}

std::string random_string(std::string alphabet, size_t n) {
  typedef boost::mt19937 rng_type;

  std::string s(n, alphabet[0]);
  rng_type rng(time(0));

  // Define random distributions
  boost::uniform_int<size_t> index_dist(0, alphabet.size()-1); // inclusive range
  boost::variate_generator<rng_type &, boost::uniform_int<size_t>> random_idx(rng, index_dist);

  for (size_t i = 0; i < n; ++i) {
    s[i] = alphabet[random_idx()];
  }
  return s;
}

}

// TODO: This may not work in VS etc... Add code to support other compilers?
#define clz(x) __builtin_clzll((x))

#endif
