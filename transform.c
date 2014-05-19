#include "transform.h"

// Doesn't reverse on bit level, reverses at the two-bit level
uint64_t block_reverse_64(uint64_t x) {
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

uint64_t block_revcomp_64(uint64_t x) {
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

// Different to block_revcomp_64 because it shifts the correct amount after
inline uint64_t reverse_complement_64(uint64_t x, uint32_t k) {
  return block_revcomp_64(x) >> (64 - k * 2);
}

uint128_t reverse_complement_128(uint128_t x, uint32_t k) {
  uint64_t temp = block_revcomp_64(x.parts.upper);
  x.parts.upper = block_revcomp_64(x.parts.lower);
  x.parts.lower = temp;
  x.v = (x.v >> (128 - k*2));
  return x;
}

