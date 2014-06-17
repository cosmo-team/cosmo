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
  size_t a_idx = 0, b_idx = 0, g_idx = 0;
  uint64_t x = 0, y, g = 0, x_prev = 0, y_prev = 0;
  char buf[k+1];

  x = get_a(a_idx);
  g = x;
  y = get_b(b_idx);
  // nothing special about the operation here, just need a different value
  x_prev = ~x;
  y_prev = ~y;

  while (a_idx < num_records && b_idx < num_records) {

    // TODO: handle duplicates in both tables
    // These y-nodes from table B will require outgoing dummy edges
    // B - A
    // if y = y_prev: set skip flag
    while (b_idx < num_records && y < g) {
      // add y to result
      // should probably use the already fetched y
      if (y != y_prev) {
        uint64_t temp = get_right_64(table_b[b_idx],k) << 2;
        sprint_kmer_acgt(buf, &temp, k);
        buf[k-1] = '$';
        fprintf(stderr, "%s\n", buf);
        count++;
      }
      if (++b_idx >= num_records) {
        fprintf(stderr, ">> BREAK\n");
        break;
      }
      y_prev = y;
      y = get_b(b_idx);
    }

    // These x-nodes from table A will require incoming dummy edges
    // A - B
    // TODO: work out how to merge the results in, as they wont be in the correct position (but will be relatively sorted)
    while (g_idx < num_records && y > g) {
      // add x to result
      if (g != x_prev) {
        uint64_t temp = table_a[g_idx] >> 2;
        sprint_kmer_acgt(buf, &temp, k);
        buf[0] = '$';
        fprintf(stderr, "%s\n", buf);
        count++;
      }
      if (++g_idx >= num_records) {
        fprintf(stderr, ">> BREAK\n");
        break;
      }
      x_prev = g;
      g = get_a(g_idx);
    }
    fprintf(stderr, ">> WELP...\n");

    a_idx = g_idx;
    x_prev = x;
    x = g;

    // These are the nodes that don't need dummy edges
    // TODO: add a partition reset for y as well - or is there a better way? since I don't care... just skip to the next change in the inner loop...
    while (a_idx < num_records && b_idx < num_records && g_idx < num_records && y == g) {
      fprintf(stderr, ">> PARTITION: %016llx\n", g);
      fprintf(stderr, ">> a_idx, b_idx, g_idx = %zu, %zu, %zu\n", a_idx, b_idx, g_idx);
      a_idx = g_idx;
      x_prev = x;
      x = g;
      // progress to end
      while (a_idx < num_records && b_idx < num_records && y == x) {
        uint64_t temp = table_a[a_idx];
        sprint_kmer_acgt(buf, &temp, k);
        fprintf(stderr, "%s\n", buf);
        if (++a_idx >= num_records) {
          fprintf(stderr, ">> BREAK\n");
          break;
        }
        x_prev = x;
        x = get_a(a_idx);
      }
      if (++b_idx >= num_records) {
        fprintf(stderr, ">> BREAK\n");
        break;
      }
      y_prev = y;
      y = get_b(b_idx);
    }
    g_idx = a_idx;
    g = x;
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
