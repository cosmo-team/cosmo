#pragma once
#ifndef UINT128_T_H
#define UINT128_T_H

#include <inttypes.h>

#define BLOCK_WIDTH 64
//uint128_t x = { 99, 20 };
//printf("%016llx%016llx\n", x.upper, x.lower);

typedef struct {
  uint64_t upper;
  uint64_t lower;
} uint128_t;

// Prototypes
uint128_t right_shift_128(uint128_t, unsigned int);
uint128_t left_shift_128(uint128_t, unsigned int);
uint128_t and_128(uint128_t, uint128_t);
uint128_t or_128(uint128_t, uint128_t);
uint128_t not_128(uint128_t);

inline uint128_t right_shift_128(uint128_t x, unsigned int distance) {
  x.lower >>= distance;
  x.lower |= (x.upper << (BLOCK_WIDTH - distance));
  x.upper >>= distance;
  return x;
}

inline uint128_t left_shift_128(uint128_t x, unsigned int distance) {
  x.upper <<= distance;
  x.upper |= (x.lower >> (BLOCK_WIDTH - distance));
  x.lower <<= distance;
  return x;
}

inline uint128_t and_128(uint128_t x, uint128_t y) {
  x.upper &= y.upper;
  x.lower &= y.lower;
  return x;
}

inline uint128_t or_128(uint128_t x, uint128_t y) {
  x.upper |= y.upper;
  x.lower |= y.lower;
  return x;
}

inline uint128_t not_128(uint128_t x) {
  x.upper = !x.upper;
  x.lower = !x.lower;
  return x;
}

#endif
