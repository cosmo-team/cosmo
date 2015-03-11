#pragma once
#ifndef _UTILITY_HPP
#define _UTILITY_HPP

#include <boost/range/iterator_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/range.hpp>
#include <climits> // For CHAR_BIT
//#include <ttmath/ttmath.h> // http://www.ttmath.org/

const uint64_t MB_TO_BYTES = 1024 * 1024;

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

// http://stackoverflow.com/questions/8511035/sequence-zip-function-for-c11
template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

// TODO: This may not work in VS etc... Add code to support other compilers?
#define clz(x) __builtin_clzll((x))

#endif
