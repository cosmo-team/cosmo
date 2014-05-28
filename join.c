#include "transform.h"
#include "join.h"

// Remove rightmost nt (2 bits)
#define get_left_64(x)  ((x) >> 2)
// Remove (leftmost) kth nt (2 bits)
#define get_right_64(x, k) (((x) << (66 - (k)*2)) >> (66 - (k)*2))

// NOTE: This part might be possible with a counting approach (since we need the result to have the same number of symbols in each column)
// similar to radix sorting, but maintain the counts for each column
// then use an in place radix sort for the first part (where we make tables a and b)
// Already sorted on previous character, so use count for 0th col on symbol from 1th col to find range to match to
// Haven't worked it out, but seems like id have to end up sorting it again, so probably faster to use two tables instead
// especially since they are both a product of a single radix sort

// left join
size_t count_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k) {
  size_t count = 0;
  size_t a_idx = 0, b_idx = 0;
  uint64_t x = 0, y = 0;

  // Iterate over table a (first k-1 nts per record) as outer loop as x
  for (;a_idx < num_records; a_idx++) {
    // NOTE: comparison done colexically - need to reverse nts in integers
    x = block_reverse_64(get_left_64(table_a[a_idx]));
    fprintf(stderr, "x: %016llx\n", x);
    // Iterate over table b (last k-1 nts per record) as inner loop as y
    // Can safely skip non-matching values for y < x, as they correspond to outgoing dummy edges
    for(;b_idx < num_records; b_idx++) {
      y = block_reverse_64(get_right_64(table_b[b_idx], k));
      fprintf(stderr, "y: %016llx\n", y);
      // If x < y then count++ or output $x, and progress x
      // Try to see how it works if we test x<=y instead (might be easier)
      if (x < y) {
        fprintf(stderr, "x < y\n");
        count++;
        break; // progress outer loop (x in table_a)
      }
      // else progress inner loop (y in table_b)
    }
    // if y gets to end, then the rest of xs require dummy edges
    // may need to check that they arent the same
    // TODO: add checks for if this x is the same as the previous one
    // (Same y check not needed since y only progresses if x >= y)
    //if (b_idx >= num_records) {
      //count += (num_records - a_idx);
      //break;
  }
  // if x gets to end, dont worry about it - all incoming dummy edges handled
  return count;
}

/*
// right join
// Don't think I actually need this function? Just pass in the tables as opposite pointers
size_t count_outgoing_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k) {
  // Iterate over table b (last k-1 nts per record) as outer loop as y
  // Iterate over table a (first k-1 nts per record) as inner loop as x
  // NOTE: comparison done colexically - need to reverse it
  // If x > y output y$ and progress x
  // else progress y
  return 0;
}
*/
