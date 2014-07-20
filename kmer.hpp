#pragma once
#ifndef KMER_HPP
#define KMER_HPP

#include <functional>
#include <string>
#include "debug.h"
#include "lut.hpp"
#include "utility.hpp"
#include "uint128_t.hpp"

#define BLOCK_WIDTH 64
#define NT_WIDTH 2
#define DNA "acgt"
#define DUMMY '$'

// Swaps G (11 -> 10) and T (10 -> 11) representation so radix ordering is lexical
// (needed because some kmer counters like DSK swap this representation, but we assume G < T
// in our de bruijn graph implementation)

struct swap_gt_f : std::unary_function<uint64_t, uint64_t> {
  inline uint64_t operator() (const uint64_t & x) const { return (x ^ ((x & 0xAAAAAAAAAAAAAAAA) >> 1)); }
} swap_gt;

// This could be used for uint128_t conversion, but it isn't needed, and including it requires messy templating
// inline uint128_t swap_gt(uint128_t x) { return uint128_t(swap_gt(x._upper), swap_gt(x._lower)); }

inline uint8_t get_nt(uint64_t block, uint8_t i) {
  // Assumes the nts are numbered from the left (which allows them to be compared as integers)
  return (block << (i * NT_WIDTH) >> (BLOCK_WIDTH - NT_WIDTH));
}

inline uint8_t get_nt(const uint128_t & block, uint8_t i) {
  const uint8_t nts_per_block = BLOCK_WIDTH/NT_WIDTH;
  uint8_t block_idx = i/nts_per_block;
  // this assumes that the block type is a multiple of 64 bits and has no other data before
  // but should be reusable for larger block types
  uint64_t block_64 = ((uint64_t*)(&block))[block_idx];
  return get_nt(block_64, i%nts_per_block);
}

template <typename T>
uint8_t get_edge_label(const T & x) {
  return get_nt(x, 0);
}

// return kmer[lo, hi) - like pythons x[lo:hi] indexing
template <typename T>
T get_range(T x, uint8_t lo = 0, uint8_t hi = -1) {
  if (hi <= lo) return T(0);
  const size_t NUM_NTS = bitwidth<T>()/NT_WIDTH;
  uint8_t shift = (hi < NUM_NTS)? (NUM_NTS - hi) * NT_WIDTH : 0;
  return (x >> shift) << (shift + lo * NT_WIDTH);
}

// It might look like the implementations for get_start_node and get_end_node
// are around the wrong way, but remember that the kmers are stored in reverse
// to enable integer comparison => colex ordering of the string representation
template <typename T>
T get_start_node(const T & x, uint8_t k) {
  return get_range(x, 1, k);
}

template <typename T>
T get_end_node(const T & x, uint8_t k) {
  return get_range(x, 0, k-1);
}

// Doesn't reverse on bit level, reverses at the two-bit level
inline uint64_t reverse_block(uint64_t x) {
  uint64_t output;

  unsigned char * p = (unsigned char *) &x;
  unsigned char * q = (unsigned char *) &output;
  q[7] = reverse_8(p[0]);
  q[6] = reverse_8(p[1]);
  q[5] = reverse_8(p[2]);
  q[4] = reverse_8(p[3]);
  q[3] = reverse_8(p[4]);
  q[2] = reverse_8(p[5]);
  q[1] = reverse_8(p[6]);
  q[0] = reverse_8(p[7]);
  return output;
}

template <class T>
struct reverse_nt : std::unary_function<T, T> {
  T operator() (const T& x) const {return reverse_block((uint64_t)x);}
};

template <>
struct reverse_nt<uint128_t> : std::unary_function<uint128_t, uint128_t> {
  inline uint128_t operator() (const uint128_t & x) const {
    return uint128_t(reverse_block(x._lower), reverse_block(x._upper));
  }
};

inline static uint64_t revcomp_block(uint64_t x) {
  uint64_t output;

  unsigned char * p = (unsigned char *) &x;
  unsigned char * q = (unsigned char *) &output;
  q[7] = revcomp_8(p[0]);
  q[6] = revcomp_8(p[1]);
  q[5] = revcomp_8(p[2]);
  q[4] = revcomp_8(p[3]);
  q[3] = revcomp_8(p[4]);
  q[2] = revcomp_8(p[5]);
  q[1] = revcomp_8(p[6]);
  q[0] = revcomp_8(p[7]);
  return output;
}

template <class T>
struct reverse_complement : std::unary_function<T, T> {
  const uint8_t _k;
  reverse_complement(uint8_t k) : _k(k) {}
  T operator() (const T& x) const {
    return revcomp_block((uint64_t)x) << (BLOCK_WIDTH - _k * NT_WIDTH);
  }
};

template <>
struct reverse_complement<uint128_t> : std::unary_function<uint128_t, uint128_t> {
  const uint8_t _k;
  reverse_complement(uint8_t k) : _k(k) {}
  uint128_t operator() (const uint128_t & x) const {
    return uint128_t(revcomp_block(x._lower), revcomp_block(x._upper)) << (bitwidth<uint128_t>() - _k * NT_WIDTH);
  }
};

template <typename T>
std::string kmer_to_string(const T & kmer_block, uint8_t max_k, uint8_t this_k = -1) {
  std::string buf(max_k, DUMMY);

  // To enable not giving a this_k value -> we can print full kmers or dummy kmers
  if (this_k > max_k) this_k = max_k;

  for (uint32_t j = 0; j < this_k; j++) {
    buf[max_k-j-1] = DNA[get_nt(kmer_block, j)];
  }
  return buf;
}

#endif
