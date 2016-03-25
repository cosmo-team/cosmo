#pragma once
#ifndef KMER_HPP
#define KMER_HPP

#include <ostream>
#include <functional>
#include <algorithm>
#include <string>
#include "utility.hpp"
#include "debug.h"
#include "lut.hpp"
#include "uint128_t.hpp"
#include <bitset>

#define BLOCK_WIDTH 64
#define NT_WIDTH 2
#define DNA_RADIX 4
#define DNA_ALPHA "acgt"
#define DUMMY_SYM '$'


typedef std::bitset<NUM_COLS> color_bv;
//FIXME: find a way to make NUM_COLS not be compile time.
// http://en.cppreference.com/w/cpp/utility/bitset says:
// "Notes: If the size of the bitset is not known at compile time, std::vector<bool> or boost::dynamic_bitset may be used."
    
void clear_bv(color_bv &bv);
void set_bit(color_bv &bv, uint32_t j);
void serialize_color_bv(std::ofstream &cfs, std::vector<color_bv>::iterator &colors, uint64_t index);
void deserialize_color_bv(std::ifstream &colorfile, color_bv &value);
// Swaps G (11 -> 10) and T (10 -> 11) representation so radix ordering is lexical
// (needed because some kmer counters like DSK swap this representation, but we assume G < T
// in our de bruijn graph implementation)

// TODO: Could probably make this faster using a static k param

static struct swap_gt_f : std::unary_function<uint64_t, uint64_t> {
  inline uint64_t operator() (const uint64_t & x) const { return (x ^ ((x & 0xAAAAAAAAAAAAAAAA) >> 1)); }
} swap_gt;

inline uint8_t get_nt(uint64_t block, uint8_t i) {
  // Assumes the nts are numbered from the left (which allows them to be compared as integers)
  return (block << (i * NT_WIDTH) >> (BLOCK_WIDTH - NT_WIDTH));
}

inline uint8_t get_nt(const uint128_t & block, uint8_t i) {
  const uint8_t nts_per_block = BLOCK_WIDTH/NT_WIDTH;
  uint8_t block_idx = i/nts_per_block;
  // this assumes that the block type is a multiple of 64 bits and has no other data before
  // but should be reusable for larger block types
  uint64_t block_64 = ((const uint64_t*const)(&block))[block_idx];
  return get_nt(block_64, i%nts_per_block);
}

template <typename kmer_t>
kmer_t clear_nt(const kmer_t & x, uint8_t i) {
  // Keep in mind that 0 is the leftmost nt, but that the leftmost is the edge (so rightmost in our diagrams)
  kmer_t mask = ~((kmer_t(3) << (bitwidth<kmer_t>::width - (1+i)*NT_WIDTH)));
  return x & mask;
}

// TODO: remove or test/fix
template <typename kmer_t>
kmer_t set_nt(const kmer_t & x, uint8_t i, uint8_t v) {
  assert(v < DNA_RADIX);
  kmer_t cleared = clear_nt(x, i);
  return cleared | (kmer_t(v) << (bitwidth<kmer_t>::width - (1+i)*NT_WIDTH));
}

template <typename T>
struct get_nt_functor : std::binary_function<T, uint8_t, uint8_t> {
  uint8_t operator() (const T & x, uint8_t i) {
    return get_nt(x, i);
  }
};

template <typename T>
uint8_t get_edge_label(const T & x) {
  return get_nt(x, 0);
}

// return kmer[lo, hi) - like pythons x[lo:hi] indexing
// interface ignores the fact that these might be stored in reverse...
// so the 0th element is still the last nucleotide
template <typename T>
T get_range(const T & x, uint8_t lo, uint8_t hi) {
  const size_t NUM_NTS = bitwidth<T>::width/NT_WIDTH;
  uint8_t shift = (hi < NUM_NTS)? (NUM_NTS - hi) * NT_WIDTH : 0;
  return (x >> shift) << (shift + lo * NT_WIDTH);
}

// It might look like the implementations for get_start_node and get_end_node
// are around the wrong way, but remember that the kmers are stored in reverse
// to enable integer comparison => colex ordering of the string representation
template <typename T>
T get_start_node(const T & x, uint8_t) {
  return get_start_node(x);
}

template <typename T>
T get_start_node(const T & x) {
  return x << NT_WIDTH;
  //return get_range(x, 1, k);
}

template <typename T>
T get_start_node_suffix(const T & x, uint8_t k) {
  return get_range(x, 1, k-1);
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
    return revcomp_block((uint64_t)x) << (64 - _k * NT_WIDTH);
  }
};

template <>
struct reverse_complement<uint128_t> : std::unary_function<uint128_t, uint128_t> {
  const uint8_t _k;
  reverse_complement(uint8_t k) : _k(k) {}
  uint128_t operator() (const uint128_t & x) const {
    return uint128_t(revcomp_block(x._lower), revcomp_block(x._upper)) << (128 - _k * NT_WIDTH);
  }
};

template <typename kmer_t>
bool is_palindrome(const kmer_t & x, uint8_t k) {
  return (k%2==0 && x == reverse_complement<kmer_t>(k)(x));
}

template <typename kmer_t>
kmer_t representative(const kmer_t & x, uint8_t k) {
  kmer_t twin = reverse_complement<kmer_t>(k)(x);
  return (x < twin)? x : twin;
}

// Shift and attach a value v on to the end (so the 0th element) of the kmer,
// and delete the last (kth) symbol (like following a path in a
// de bruijn graph)
template <typename kmer_t>
kmer_t follow_edge(const kmer_t & x, uint8_t k, uint8_t v) {
  kmer_t y(x);
  // unset character
  y = set_nt(y, k-1, 0);
  y >>= NT_WIDTH;
  // set character
  y = set_nt(y, 0, v);
  return y;
}

template <typename T>
std::string kmer_to_string(const T & kmer_block, uint8_t max_k, uint8_t this_k = -1) {
  std::string buf(max_k, DUMMY_SYM);

  // To enable not giving a this_k value -> we can print full kmers or dummy kmers
  if (this_k > max_k) this_k = max_k;

  for (uint32_t j = 0; j < this_k; j++) {
    buf[max_k-j-1] = DNA_ALPHA[get_nt(kmer_block, j)];
  }
  return buf;
}

template <typename kmer_t>
void convert_representation(const kmer_t * kmers_in, kmer_t * kmers_out, size_t num_kmers, bool swap) {
  // Swap G and T and reverses the nucleotides so that they can
  // be compared and sorted as integers to give colexicographical ordering
  if (swap)
    std::transform((const uint64_t*const )kmers_in, (const uint64_t*const)(kmers_in + num_kmers), (uint64_t*)kmers_out, swap_gt);
  std::transform(kmers_in, kmers_in + num_kmers, kmers_out, reverse_nt<kmer_t>());
}

// Convenience function to print array of kmers
template <typename kmer_t>
void print_kmers(std::ostream & out, kmer_t * kmers, size_t num_kmers, uint32_t k, uint8_t * dummy_lengths = 0) {
  for (size_t i = 0; i<num_kmers; i++) {
    out << kmer_to_string(kmers[i], k, ((dummy_lengths)? dummy_lengths[i]:k)) << std::endl;
  }
}

// Longest common suffix
template <typename kmer_t>
size_t lcs(const kmer_t & a, const kmer_t & b, size_t k) {
  static const size_t num_blocks = bitwidth<kmer_t>::width/BLOCK_WIDTH;
  //kmer_t x = a ^ b; // Do this in the loop below instead to potentially save some cycles
  uint64_t * p = (uint64_t*)&a;
  uint64_t * q = (uint64_t*)&b;
  size_t total = 0;
  // This should unroll. for 128 bits its only 2 iters
  for (size_t i = 0; i < num_blocks; i++) {
    if (p[i] == q[i]) {
      total += BLOCK_WIDTH;
      continue;
    }
    else {
      // COUNT *LEADING* ZEROS - we store kmers backwards
      total += clz(p[i] ^ q[i]);
      break;
    }
  }
  total /= NT_WIDTH;
  return std::min(total, k);
}

template <typename kmer_t>
size_t node_lcs(const kmer_t & a, const kmer_t & b, size_t k) {
  //assert(k>0); // Shouldnt be called this way, but we also minus 1 down below...
  if (k == 0) return 0;
  kmer_t x(a);
  kmer_t y(b);
  // TODO: make position-templated set_nt()
  ((uint64_t*)&x)[0] &= 0x3FFFFFFFFFFFFFFF;
  ((uint64_t*)&y)[0] &= 0x3FFFFFFFFFFFFFFF;
  return lcs(x,y,k) - 1;
}

#endif
