#pragma once
#ifndef KMER_HPP
#define KMER_HPP

#include <string>
#include "lut.hpp"
#include "utility.hpp"
#include "uint128_t.hpp"

#define BLOCK_WIDTH 64
#define NT_WIDTH 2
#define DNA "acgt"
#define DUMMY '$'

// Swaps G (11 -> 11) and T (10 -> 11) representation so radix ordering is lexical
inline uint64_t swap_gt(uint64_t x) {
  return (x ^ ((x & 0xAAAAAAAAAAAAAAAA) >> 1));
}

inline uint128_t swap_gt(const uint128_t & x) {
  return uint128_t(swap_gt(x._upper), swap_gt(x._lower));
}

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

// Doesn't reverse on bit level, reverses at the two-bit level
inline uint64_t reverse_nt(uint64_t x) {
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

inline uint128_t reverse_nt(const uint128_t & x) {
  return uint128_t(reverse_nt(x._lower), reverse_nt(x._upper));
}

inline static uint64_t block_revcomp(uint64_t x) {
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

inline uint64_t reverse_complement(uint64_t x, uint8_t k) {
  return block_revcomp(x) << (BLOCK_WIDTH - k * NT_WIDTH);
}

inline uint128_t reverse_complement(const uint128_t & x, uint8_t k) {
  // should be rewritable for larger block types using map and reverse for arrays
  return uint128_t(block_revcomp(x._lower), block_revcomp(x._upper)) << (bitwidth<uint128_t>() - k * NT_WIDTH);
}

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
