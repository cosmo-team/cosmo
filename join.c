#include "transform.h"
#include "io.h"
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

// A - B
size_t count_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k) {
  #define get_a(i) (block_reverse_64(get_left_64(table_a[(i)])))
  #define get_b(i) (block_reverse_64(get_right_64(table_b[(i)], k)))
  size_t count = 0;
  size_t a_idx = 0, b_idx = 0;
  uint64_t x = 0, y = 0;
  char buf[k+1];

  while (a_idx < num_records && b_idx < num_records) {
    x = get_a(a_idx);
    y = get_b(b_idx);

    // TODO: handle duplicates in both tables
    // These y-nodes from table B will require outgoing dummy edges
    while (y < x) {
      // add y to result
      /*
      uint64_t temp = get_right_64(table_b[b_idx],k) << 2;
      sprint_kmer_acgt(buf, &temp, k);
      buf[k-1] = '$';
      fprintf(stderr, "%s\n", buf);
      // TODO: add checks that we aren't going past the end of the arrays
      */
      y = get_b(++b_idx);
    }

    // These x-nodes from table A will require incoming dummy edges
    while (y > x) {
      // add x to result
      /*
      uint64_t temp = table_a[a_idx] >> 2;
      sprint_kmer_acgt(buf, &temp, k);
      buf[0] = '$';
      fprintf(stderr, "%s\n", buf);
      */
      x = get_a(++a_idx);
    }

    // These are the nodes that don't need dummy edges
    // TODO: understand why these are printing when I have output a dummy - probably duplication?
    while (y == x) {
      // progress to end
      while (y == x) {
        uint64_t temp = table_a[a_idx];
        sprint_kmer_acgt(buf, &temp, k);
        fprintf(stderr, "%s\n", buf);
        // if outputting, add x to result
        x = get_a(++a_idx);
      }
      y = get_b(++b_idx);
    }
  }
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
