#pragma once
#ifndef UTILITY_HPP
#define UTILITY_HPP
// Utility functions

// For CHAR_BIT, which might not be 8 on certain systems
#include <climits>
#include <string>
#include <boost/filesystem.hpp>
#include "debug.hpp"

namespace cosmo {

namespace fs = boost::filesystem;
using fs::extension;
using fs::basename;

template <typename T>
struct bitwidth {
  const static size_t width = sizeof(T) * CHAR_BIT;
};

template <typename T>
T & deconst(const T & x) {
  return const_cast<T&>(x);
}

template <typename T>
struct identity {
  const T& operator()(const T& value) {
    return value;
  }
};

enum class input_format { raw, dsk, kmc };
const std::string input_format_strings[] = { "raw", "DSK1", "KMC2" };

inline input_format get_format(const std::string & filename) {
  std::string ext = extension(filename);
  if (ext == ".solid_kmers_binary") return input_format::dsk;
  if (ext == ".kmc" || ext == ".kmc_pre" || ext == ".kmc_suf") return input_format::kmc;
  return input_format::raw;
}

inline bool probably_list_of_files(const std::string & filename, input_format * fmt = nullptr) {
  size_t lines_left = 100;
  const size_t max_line_length = 100000;
  char buffer[max_line_length];
  std::ifstream in(filename);
  while (lines_left-- && in.getline(buffer, max_line_length)) {
    std::string path(buffer);
    // kmc_pre + kmc_suf are needed, but we only test for one
    auto extensions = { "", ".solid_kmers_binary", ".kmc_pre" };

    bool none_matched = true;
    for (auto ext : extensions) {
      auto p = path + ext;
      if (boost::filesystem::exists(p)) {
        none_matched = false;
        if (fmt) *fmt = get_format(p);
        break;
      }
    }

    if (none_matched) {
      in.close();
      return false;
    }
  }
  in.close();
  return true;
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

}

// TODO: This may not work in VS etc... Add code to support other compilers?
#define clz(x) __builtin_clzll((x))

#endif
