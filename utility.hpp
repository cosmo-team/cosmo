#pragma once
#ifndef _UTILITY_HPP
#define _UTILITY_HPP

#include <climits> // For CHAR_BIT
//#include <ttmath/ttmath.h> // http://www.ttmath.org/

//typedef ttmath::UInt<2> uint128_t;
typedef __uint128_t uint128_t;

//typedef ttmath::UInt<2> uint128_t;
typedef __uint128_t uint128_t;

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

// TODO: This may not work in VS etc... Add code to support other compilers?
#define clz(x) __builtin_clzll((x))

#endif
