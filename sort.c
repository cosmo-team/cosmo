#include "sort.h"

#define BASE 4
#define DIGIT_WIDTH 2
#define DIGIT_MASK 0x3
#define BLOCK_WIDTH 64

#define get_digit_64(x,i) (((x) >> i * DIGIT_WIDTH) & DIGIT_MASK)
#define get_block_idx(i) (1 - (i)*DIGIT_WIDTH/BLOCK_WIDTH)
#define get_digit_128(x,i) (((uint64_t*)(x))[get_block_idx(i)] >> (((i)*DIGIT_WIDTH)%BLOCK_WIDTH))

void lsd_radix_sort_64(uint64_t * a, uint64_t * b, size_t num_records, uint32_t k) {
  // Init counts
  size_t bases[BASE];

  // For each 2-bit-digit from the right
  for (ssize_t digit_pos = k-1; digit_pos >= 0; digit_pos--) {
    // Reset counts
    for (int c = 0; c < BASE; c++) bases[c] = 0;
    // Count how many of each digit there are in each position
    for (size_t i = 0; i < num_records; i++) bases[get_digit_64(a[i], digit_pos)]++;

    // prefix sum the counts to give us starting positions
    for (int c = 1; c < BASE; c++) bases[c] += bases[c-1];

    // Stably copy each element into its corresponding location
    // Done in reverse to simplify sub-array calculations
    for (ssize_t i = num_records - 1; i >= 0; i--) {
      int x = get_digit_64(a[i], digit_pos);
      b[--bases[x]] = a[i];
    }

    // swap array ptrs
    uint64_t * temp = a;
    a = b;
    b = temp;
  }
  // TODO: could maybe return an int to show which one the result is in, and swap the values instead of copy (above)
  // so we can use it to sort the other array when we need to
  // Copy data from temp array if we need to
  // If we didnt have an even number of iterations, then a will be in b still
  // but the ptrs will have swapped
  if (k%2!=0) {
    memcpy(b, a, num_records * sizeof(uint64_t));
  }
}
