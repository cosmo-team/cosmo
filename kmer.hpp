#pragma once
#ifndef KMER_HPP
#define KMER_HPP

//#include <ostream>
//#include <functional>
//#include <algorithm>
#include <string>
#include "utility.hpp"
#include "lut.hpp"

#define BLOCK_WIDTH 64
#define NT_WIDTH 2
#define DNA_RADIX 4
#define DNA_ALPHA "acgt"
#define DUMMY_SYM '$'

// Swaps G (11 -> 10) and T (10 -> 11) representation so radix ordering is lexical
inline uint64_t _swap_gt_64(const uint64_t x) {
  return (x ^ ((x & 0xAAAAAAAAAAAAAAAA) >> 1));
}

template <class T>
struct swap_gt : std::unary_function<T, T> {
  inline T operator() (const T& x) const {
    T temp(x);
    uint64_t * blocks = reinterpret_cast<uint64_t*>(&temp);
    for (int i = 0; i < sizeof(T)/sizeof(uint64_t); i++) {
      blocks[i] = _swap_gt_64(blocks[i]);
    }
    return temp;
  }
};

// Reverses block at two-bit (nt) level
inline uint64_t _reverse_nt_64(const uint64_t x) {
  uint64_t output;
  unsigned char * p = (unsigned char *) &x;
  unsigned char * q = (unsigned char *) &output;

  for (size_t i = 0; i < 8; i++) {
    q[8 - i - 1] = reverse_8(p[i]);
  }

  return output;
}

template <class T>
struct reverse_nt : std::unary_function<T, T> {
  T operator() (const T& x) const {
    size_t num_blocks = sizeof(T)/sizeof(uint64_t);

    T result(x);
    for (size_t i = 0; i < num_blocks; i++) {
      uint64_t block = ((uint64_t*)&x)[i];
      ((uint64_t*)&result)[num_blocks - i - 1] = _reverse_nt_64(block);
    }

    return result;
  }
};

template <typename T>
inline uint8_t get_nt(const T & x, size_t i) {
  size_t num_blocks = sizeof(T)/sizeof(uint64_t);
  size_t nts_per_block = 32;

  // Nucleotides are numbered starting from 0 on the right of kmer,
  // and are numbered from MSB as to be sorted in colex order.
  // 0.....k-1
  // MSB...LSB
  size_t block_idx = num_blocks - i/nts_per_block - 1;
  size_t nt_idx    = i%nts_per_block;
  uint64_t block   = ((uint64_t*)&x)[block_idx];
  return (block << (nt_idx * 2)) >> 62;
}

template <typename T>
std::string kmer_to_string(const T & kmer, size_t max_k, size_t this_k = -1) {
  std::string buf(max_k, DUMMY_SYM);

  // To enable not giving a this_k value -> we can print full kmers or dummy kmers
  if (this_k > max_k) this_k = max_k;

  for (size_t nt_pos = 0; nt_pos < this_k; nt_pos++) {
    buf[max_k-nt_pos-1] = DNA_ALPHA[get_nt(kmer, nt_pos)];
  }
  return buf;
}

inline static uint64_t _reverse_complement_64(uint64_t x) {
  uint64_t output;

  unsigned char * p = (unsigned char *) &x;
  unsigned char * q = (unsigned char *) &output;

  for (size_t i = 0; i < 8; i++) {
    q[8 - i - 1] = revcomp_8(p[i]);
  }

  return output;
}

template <class T>
struct reverse_complement : std::unary_function<T, T> {
  const size_t _k;
  reverse_complement(size_t k) : _k(k) {}
  T operator() (const T& x) const {
    size_t num_blocks = sizeof(T)/sizeof(uint64_t);

    T result(x);
    for (size_t i = 0; i < num_blocks; i++) {
      uint64_t block = ((uint64_t*)&x)[i];
      ((uint64_t*)&result)[num_blocks - i - 1] = _reverse_complement_64(block);
    }
    return result << (sizeof(T)*8 - _k * NT_WIDTH);
  }
};

template <typename T>
inline T get_start_node(const T & x) {
  return x << NT_WIDTH;
}

template <typename T>
inline uint8_t get_edge_label(const T & x) {
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

template <typename T>
T get_end_node(const T & x, uint8_t k) {
  return get_range(x, 0, k-1);
}

/*
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
*/

#endif
