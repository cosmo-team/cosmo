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
  uint64_t x = 0, y = 0, x_prev = 0, y_prev = 0;
  char buf[k+1];

  if (num_records == 0) return 0;

  x = get_a(a_idx);
  y = get_b(b_idx);
  // nothing special about the NOT operation here, just need a guaranteed different value to start with
  x_prev = ~x;
  y_prev = ~y;

  while (a_idx < num_records && b_idx < num_records) {
    // B - A: These y-nodes from table B will require outgoing dummy edges
    while (b_idx < num_records && y < x) {
      // add y to result
      if (y != y_prev) {
        // should probably use the already fetched y
        uint64_t temp = get_right_64(table_b[b_idx],k) << 2;
        sprint_kmer_acgt(buf, &temp, k);
        buf[k-1] = '$';
        fprintf(stderr, "%s\n", buf);
        count++;
      }
      if (++b_idx >= num_records) break;
      y_prev = y;
      y = get_b(b_idx);
    }

    // A - B: These x-nodes from table A will require incoming dummy edges
    // TODO: work out how to merge the results in, as they wont be in the correct position (but will be relatively sorted)
    while (a_idx < num_records && y > x) {
      // add x to result
      if (x != x_prev) {
        uint64_t temp = table_a[a_idx] >> 2;
        sprint_kmer_acgt(buf, &temp, k);
        buf[0] = '$';
        fprintf(stderr, "%s\n", buf);
        count++;
      }
      if (++a_idx >= num_records) break;
      x_prev = x;
      x = get_a(a_idx);
    }

    // These are the nodes that don't need dummy edges
    while (a_idx < num_records && b_idx < num_records && y == x) {
      x_prev = x;
      y_prev = y;

      // Scan to next non-equal (outputing all table entries?)
      while (a_idx < num_records && x == x_prev) {
        uint64_t temp = table_a[a_idx];
        sprint_kmer_acgt(buf, &temp, k);
        fprintf(stderr, "%s\n", buf);
        if (++a_idx >= num_records) break;
        x_prev = x;
        x = get_a(a_idx);
      }

      // Scan to next non-equal
      while (b_idx < num_records && y == y_prev) {
        if (++b_idx >= num_records) break;
        y_prev = y;
        y = get_b(b_idx);
      }
    }
  }
  return count;
}

