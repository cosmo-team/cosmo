#pragma once
#ifndef _UTILITY_HPP
#define _UTILITY_HPP

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

template <typename T>
T & deconst(const T & x) {
  return const_cast<T&>(x);
}

// This may not work in VS etc... Add code to support other compilers?
#define clz(x) __builtin_clzll((x))

#endif
