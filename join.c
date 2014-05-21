#include "join.h"

// Remove rightmost nt (2 bits)
#define get_left_64(x)  ((x) >> 2)
// Remove (leftmost) kth nt (2 bits)
#define get_right_64(x, k) (((x) << (62 - (k)*2)) >> (62 - (k)*2))

// NOTE: This part might be possible with a counting approach (since we need the result to have the same number of symbols in each column)

// left join
size_t count_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k) {
  size_t count = 0;
  size_t a_idx = 0, b_idx = 0;
  uint64_t x, y;
  // Iterate over table a (first k-1 nts per record) as outer loop as x
  for (;a_idx < num_records; a_idx++) {
    x = get_left_64(table_a[a_idx]);
    // Iterate over table b (last k-1 nts per record) as inner loop as y
    for(;b_idx < num_records; b_idx++) {
      y = get_right_64(table_b[b_idx], k);
      // If x < y then count++ or output $x, and finally progress x
      // NOTE: comparison done colexically - need to reverse nts in integers
      // TODO: Reverse nts in integers
      // TODO: Move comparator function from sort into a non-dereferencing version
      if (x < y) {
        count++;
        break; // progress outer loop (x in table_a)
      }
      // else progress inner loop (y in table_b)
    }
    // if y gets to end, then the rest of xs require dummy edges
    // may need to check that they arent the same
    if (b_idx >= num_records) {
      count += (num_records - a_idx);
      break;
    }
  }
  // if x gets to end, dont worry about it - all incoming dummy edges handled
  return count;
}

// right join
size_t count_outgoing_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k) {
  // Iterate over table b (last k-1 nts per record) as outer loop as y
  // Iterate over table a (first k-1 nts per record) as inner loop as x
  // NOTE: comparison done colexically - need to reverse it
  // If x > y output y$ and progress x
  // else progress y
  return 0;
}
