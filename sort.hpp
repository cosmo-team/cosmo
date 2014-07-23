#ifndef SORT_HPP
#define SORT_HPP

// Radix sorts colex(range (hi, lo]) of each record table_a
// using table_b as the temporary table, and writing the new ptrs to
// new_a and new_b. new_a will point to the final result, while
// new_b will be the results of the second-last iteration.
template <int base, typename T, typename F>
void colex_partial_radix_sort(T * a, T * b, size_t num_records, uint32_t lo, uint32_t hi, T ** new_a, T ** new_b, F get_digit) {
  assert(hi > lo);
  // Init counts
  size_t bases[base];

  // For each 2-bit-digit from the right
  for (ssize_t digit_pos = hi-1; digit_pos >= lo; digit_pos--) {
    // Reset counts
    for (int c = 0; c < base; c++) bases[c] = 0;
    // Count how many of each digit there are in each position
    for (size_t i = 0; i < num_records; i++) bases[get_digit(a[i], digit_pos)]++;

    // prefix sum the counts to give us starting positions
    for (int c = 1; c < base; c++) bases[c] += bases[c-1];

    // Stably copy each element into its corresponding location
    // Done in reverse to simplify sub-array calculations
    for (ssize_t i = num_records - 1; i >= 0; i--) {
      int x = get_digit(a[i], digit_pos);
      b[--bases[x]] = a[i];
    }

    // swap array ptrs
    T * temp = a;
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

#endif
