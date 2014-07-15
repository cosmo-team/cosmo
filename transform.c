
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
  uint64_t temp = x.upper;
  x.upper = block_revcomp_64(x.lower);
  x.lower = block_revcomp_64(temp);
  //x.v >>= (128 - k*2);
  x = right_shift_128(x, 128 - k*2);
  return x;
}

void add_reverse_complements(const uint64_t * kmers_in, uint64_t * kmers_out, size_t num_kmers, uint32_t k) {
  assert(k <= 64);
  if (k <= 32) {
    for (size_t i = 0; i < num_kmers; i++) {
      kmers_out[i] = reverse_complement_64(kmers_in[i], k);
    }
  }
  else {
    for (size_t i = 0; i < num_kmers; i++) {
      uint128_t * in = (uint128_t*) kmers_in;
      uint128_t * out = (uint128_t*) kmers_out;
      out[i] = reverse_complement_128(in[i], k);
    }
  }
}
