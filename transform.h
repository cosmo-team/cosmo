#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "common.h"
#include "uint128_t.h"
#include "lut.h"

// Swaps G (11 -> 11) and T (10 -> 11) representation so radix ordering is lexical
#define swap_gt_64(x) ((x) ^ (((x) & 0xAAAAAAAAAAAAAAAA) >> 1))

uint64_t block_revcomp_64(uint64_t x);
uint64_t block_reverse_64(uint64_t x);
uint64_t reverse_complement_64(uint64_t x, uint32_t k);
uint128_t reverse_complement_128(uint128_t x, uint32_t k);

#endif
