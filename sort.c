#include "sort.h"
#include "transform.h"

#define BASE 4
#define DIGIT_WIDTH 2
#define DIGIT_MASK 0x3
#define BLOCK_WIDTH 64

#define get_digit_64(x,i) (((x) >> i * DIGIT_WIDTH) & DIGIT_MASK)
#define get_block_idx(i) (1 - (i)*DIGIT_WIDTH/BLOCK_WIDTH)
#define get_digit_128(x,i) (((uint64_t*)(x))[get_block_idx(i)] >> (((i)*DIGIT_WIDTH)%BLOCK_WIDTH))

// Base 4 radix sorts colex(range (k, 0]) of each record table_a
// using table_b as the temporary table, and writing the new ptrs to
// new_a and new_b. new_a will point to the final result, while
// new_b will be the results of the second-last iteration.
void colex_partial_radix_sort_64(uint64_t * a, uint64_t * b, size_t num_records, uint32_t k, uint32_t j, uint64_t ** new_a, uint64_t** new_b) {
  // TODO: add optional size check -> each position that is in a > size position is a $
  // TODO: support lo and hi swapping to reverse order
  assert(k > j);
  // Init counts
  size_t bases[BASE];

  // For each 2-bit-digit from the right
  for (ssize_t digit_pos = k-1; digit_pos >= j; digit_pos--) {
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
  // Want a to be the final result, b to be the second-last iteration.
  // If we didnt have an even number of iterations, then desired contents of a will be in array b still,
  // but the ptrs will have swapped.
  // Could use memcpy(b, a, num_records * sizeof(uint64_t)); instead, but if we use out-ptrs we save a memcpy
  // Besides, we have a use for the second-last iteration, which will be in b, so the client code can just pass in
  // new ptrs to be overwritten.
  *new_a = a;
  *new_b = b;
}

void colex_varlen_partial_radix_sort_64(uint64_t * a, uint64_t * b, unsigned char * lengths_a, unsigned char * lengths_b, size_t num_records, uint32_t k, uint32_t j, uint64_t ** new_a, uint64_t** new_b, unsigned char ** new_lengths_a, unsigned char ** new_lengths_b) {
  // TODO: add optional size check -> each position that is in a > size position is a $
  // TODO: support lo and hi swapping to reverse order
  assert(k > j);
  // Init counts (+1 for $ sign)
  size_t bases[BASE + 1];

  // For each 2-bit-digit from the right
  for (ssize_t digit_pos = k-1; digit_pos >= j; digit_pos--) {
    // Reset counts
    for (int c = 0; c < BASE + 1; c++) bases[c] = 0;
    // Count how many of each digit there are in each position
    for (size_t i = 0; i < num_records; i++) 
      if (lengths_a[i] < digit_pos) {
        bases[0]++;
      } else {
        bases[get_digit_64(block_reverse_64(a[i]), digit_pos)+1]++;
      }

    // prefix sum the counts to give us starting positions
    for (int c = 1; c < BASE+1; c++) bases[c] += bases[c-1];

    // Stably copy each element into its corresponding location
    // Done in reverse to simplify sub-array calculations
    for (ssize_t i = num_records - 1; i >= 0; i--) {
      int x = (lengths_a[i] < digit_pos)? 0 : get_digit_64(block_reverse_64(a[i]), digit_pos) + 1;
      b[--bases[x]] = a[i];
      lengths_b[bases[x]] = lengths_a[i];
    }

    // swap array ptrs
    uint64_t * temp = a;
    a = b;
    b = temp;
    unsigned char * temp_lengths = lengths_a;
    lengths_a = lengths_b;
    lengths_b = temp_lengths;
  }
  // Want a to be the final result, b to be the second-last iteration.
  // If we didnt have an even number of iterations, then desired contents of a will be in array b still,
  // but the ptrs will have swapped.
  // Could use memcpy(b, a, num_records * sizeof(uint64_t)); instead, but if we use out-ptrs we save a memcpy
  // Besides, we have a use for the second-last iteration, which will be in b, so the client code can just pass in
  // new ptrs to be overwritten.
  *new_a = a;
  *new_b = b;
  *new_lengths_a = lengths_a;
  *new_lengths_b = lengths_b;
}


