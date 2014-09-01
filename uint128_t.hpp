#pragma once
#ifndef UINT128_T_HPP
#define UINT128_T_HPP

#include <climits> // For CHAR_BIT

// Possibly handy function for users of this type
template <typename T>
struct bitwidth {
  const static size_t width = sizeof(T) * CHAR_BIT;
};

/* UINT128_T
 * Useful for wider blocks (e.g. when 32 < k <= 64 in a DNA kmer program)
 * Based off https://github.com/calccrypto/uint128_t/ but with fewer functions (not all operations are required:
 * arithmetic and printing functions are removed)
 * and rewritten following guidelines from http://courses.cms.caltech.edu/cs11/material/cpp/donnie/cpp-ops.html
 * Also, all inline for easy inclusion in a project, and the member uint64_ts are left public (which makes sense
 * for a simple binary block type, imo).
 *
 * The above was distributed under the MIT license, and the message has been included below.
 *
 * Note: if extended to support wider blocks, the template functions take a RHS assumed to be smaller than 128
 * bits, so we would need to check the width of the incoming type
 * Also, shifting negatively isn't supported. The compiler gives a warning but some people might want to do it...
 * Sorry! Maybe in another version.
 */

// ORIGINAL UINT128_T LICENSE MESSAGE:

/*
uint128_t.h
An unsigned 128 bit integer type for C++
Copyright (c) 2014 Jason Lee @ calccrypto at gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

With much help from Auston Sterling

Thanks to Stefan Deigmüller for finding
a bug in operator*.

Thanks to François Dessenne for convincing me
to do a general rewrite of this class.
*/

struct uint128_t {
  private:
  // These are included to make operations easier to read
  // And may assist with generalizing later to support uint256_t, etc...
  typedef uint64_t block_t;
  const static size_t block_width = sizeof(block_t) * CHAR_BIT;
  const static size_t num_blocks = 2;
  const static size_t total_width = num_blocks * block_width;

  // These remain public because this class is mainly for bit manipulation
  // and some algorithms that use it might require the speed of direct member access (e.g. reversal)
  public:
  block_t _upper, _lower;

  // Constructors
  uint128_t() : _upper(0), _lower(0) {}
  uint128_t(const uint128_t & rhs) : _upper(rhs._upper), _lower(rhs._lower) {}
  template <typename T>
  uint128_t(const T & rhs) : _upper(0), _lower(rhs) {}
  template <typename S, typename T>
  uint128_t(const S & upper_rhs, const T & lower_rhs) : _upper(upper_rhs), _lower(lower_rhs) {}

  // Assignment Operator
  uint128_t & operator=(const uint128_t & rhs) {
    _upper = rhs._upper;
    _lower = rhs._lower;
    return *this;
  }

  template <typename T>
  uint128_t & operator=(const T & rhs) {
    _upper = 0;
    _lower = rhs;
    return *this;
  }

  // Bitwise Operators
  uint128_t & operator&=(const uint128_t & rhs) {
    _upper &= rhs._upper;
    _lower &= rhs._lower;
    return *this;
  }

  uint128_t & operator|=(const uint128_t & rhs) {
    _upper |= rhs._upper;
    _lower |= rhs._lower;
    return *this;
  }

  uint128_t & operator^=(const uint128_t & rhs) {
    _upper ^= rhs._upper;
    _lower ^= rhs._lower;
    return *this;
  }

  const uint128_t operator&(const uint128_t & rhs) const { return uint128_t(*this) &= rhs; }
  const uint128_t operator|(const uint128_t & rhs) const { return uint128_t(*this) |= rhs; }
  const uint128_t operator^(const uint128_t & rhs) const { return uint128_t(*this) ^= rhs; }
  const uint128_t operator~() const { return uint128_t(~_upper, ~_lower); }

  template <typename T>
  uint128_t & operator&=(const T & rhs) {
    _upper = 0;
    _lower &= rhs;
    return *this;
  }

  template <typename T>
  uint128_t & operator|=(const T & rhs){
    _lower |= (block_t) rhs;
    return *this;
  }

  template <typename T>
  uint128_t & operator^=(const T & rhs){
    _lower ^= (block_t) rhs;
    return *this;
  }

  template <typename T>
  uint128_t operator&(const T & rhs) const { return uint128_t(*this) &= rhs; }
  template <typename T>
  uint128_t operator|(const T & rhs) const { return uint128_t(*this) |= rhs; }
  template <typename T>
  uint128_t operator^(const T & rhs) const { return uint128_t(*this) ^= rhs; }

  // Shift Operators
  template <typename T>
  uint128_t operator<<(const T & rhs) const {
    uint64_t shift = uint64_t(rhs);
    if (shift >= total_width) return uint128_t(0);
    else if (shift == block_width) return uint128_t(_lower, 0);
    else if (shift == 0) return *this;
    else if (shift < block_width) return uint128_t((_upper << shift) + (_lower >> (block_width - shift)), _lower << shift);
    else if ((total_width > shift) && (shift > block_width)) return uint128_t(_lower << (shift - block_width), 0);
    else return uint128_t(0);
  }

  template <typename T>
  uint128_t operator>>(const T & rhs) const {
    uint64_t shift = uint64_t(rhs);
    if (shift >= total_width) return uint128_t(0);
    else if (shift == block_width) return uint128_t(0, _upper);
    else if (shift == 0) return *this;
    else if (shift < block_width) return uint128_t(_upper >> shift, (_upper << (block_width - shift)) | (_lower >> shift));
    else if ((total_width > shift) && (shift > block_width)) return uint128_t(0, (_upper >> (shift - block_width)));
    else return uint128_t(0);
  }

  uint128_t operator<<(const uint128_t & rhs) const {
    if ((bool)rhs._upper) return uint128_t(0);
    else return *this << (uint64_t) rhs;
  }

  uint128_t operator>>(const uint128_t & rhs) const {
    if ((bool)rhs._upper) return uint128_t(0);
    else return *this >> (uint64_t) rhs;
  }

  uint128_t & operator<<=(const uint128_t & rhs) {
    *this = *this << rhs;
    return *this;
  }
 
  template <typename T>
  uint128_t & operator<<=(const T & rhs) {
    *this = *this << (uint64_t) rhs;
    return *this;
  }

  uint128_t & operator>>=(const uint128_t & rhs) {
    *this = *this >> rhs;
    return *this;
  }
  
  template <typename T>
  uint128_t & operator>>=(const T & rhs) {
    *this = *this >> (uint64_t) rhs;
    return *this;
  }

  // Comparison Functions
  bool operator!=(const uint128_t & rhs) const { return (_upper != rhs._upper || _lower != rhs._lower); }
  bool operator==(const uint128_t & rhs) const { return !(*this != rhs); }
  bool operator>(const uint128_t & rhs) const { return (_upper == rhs._upper)? (_lower > rhs._lower) : (_upper > rhs._upper); }
  bool operator<(const uint128_t & rhs) const { return (_upper == rhs._upper)? (_lower < rhs._lower) : (_upper < rhs._upper); }
  bool operator>=(const uint128_t & rhs) const { return *this > rhs || *this == rhs; }
  bool operator<=(const uint128_t & rhs) const { return *this < rhs || *this == rhs; }

  template <typename T> bool operator==(const T & rhs) const {
    return (!_upper && (_lower == (uint64_t) rhs));
  }

  template <typename T> bool operator!=(const T & rhs) const {
    return (_upper || (_lower != (uint64_t) rhs));
  }

  template <typename T> bool operator>(const T & rhs) const {
    return (_upper || (_lower > (uint64_t) rhs));
  }

  template <typename T> bool operator<(const T & rhs) const {
    return (!_upper)? (_lower < (uint64_t) rhs) : false;
  }

  template <typename T> bool operator>=(const T & rhs) const {
    return ((*this > rhs) | (*this == rhs));
  }

  template <typename T> bool operator<=(const T & rhs) const {
    return ((*this < rhs) | (*this == rhs));
  }

  // Casting
  operator bool() const { return (bool) _lower; }
  operator char() const { return (char) _lower; }
  operator int() const { return (int) _lower; }
  operator uint8_t() const { return (uint8_t) _lower; }
  operator uint16_t() const { return (uint16_t) _lower; }
  operator uint32_t() const { return (uint32_t) _lower; }
  operator uint64_t() const { return (uint64_t) _lower; }
};

// lhs type T as first arguemnt
// If the output is not a bool, casts to type T

// Bitwise Operators
// I get errors if 
/*
template <typename T> T operator&(const T & lhs, const uint128_t & rhs){
  return (T) (lhs & (T) rhs._lower);
}

template <typename T> T operator|(const T & lhs, const uint128_t & rhs){
  return (T) (lhs | (T) rhs._lower);
}

template <typename T> T operator^(const T & lhs, const uint128_t & rhs){
  return (T) (lhs ^ (T) rhs._lower);
}

template <typename T> T operator&=(T & lhs, const uint128_t & rhs){
  lhs &= (T) rhs._lower; return lhs;
}

template <typename T> T operator|=(T & lhs, const uint128_t & rhs){
  lhs |= (T) rhs._lower; return lhs;
}

template <typename T> T operator^=(T & lhs, const uint128_t & rhs){
  lhs ^= (T) rhs._lower; return lhs;
}

// Comparison Operators
template <typename T> bool operator==(const T & lhs, const uint128_t & rhs){
  return (!rhs._upper && ((uint64_t) lhs == rhs._lower));
}

template <typename T> bool operator!=(const T & lhs, const uint128_t & rhs){
  return (rhs._upper || ((uint64_t) lhs != rhs._lower));
}

template <typename T> bool operator>(const T & lhs, const uint128_t & rhs){
  return (!rhs._upper && ((uint64_t) lhs > rhs._lower));
}

template <typename T> bool operator<(const T & lhs, const uint128_t & rhs){
  return (rhs._upper)? true : (uint64_t) lhs < rhs._lower;
}

template <typename T> bool operator>=(const T & lhs, const uint128_t & rhs){
  return (rhs._upper)? false : (uint64_t) lhs >= rhs._lower;
}

template <typename T> bool operator<=(const T & lhs, const uint128_t & rhs){
  return (rhs._upper)? true : (uint64_t) lhs <= rhs._lower;
}
*/

#endif
