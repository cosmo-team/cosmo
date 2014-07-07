#include "transform.h"
#include "io.h"
#include "sort.h"
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

size_t count_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k) {
  #define get_a(i) (block_reverse_64(get_left_64(table_a[(i)])))
  #define get_b(i) (block_reverse_64(get_right_64(table_b[(i)], k)))
  size_t count = 0;
  size_t a_idx = 0, b_idx = 0;
  uint64_t x = 0, y = 0, x_prev = 0, y_prev = 0;
  //char buf[k+1];

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
      if (++b_idx >= num_records) break;
      y = get_b(b_idx);
    }

    // A - B: These x-nodes from table A will require incoming dummy edges
    while (a_idx < num_records && y > x) {
      // add x to result
      if (x != x_prev) {
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

void get_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k, uint64_t * incoming_dummies, size_t num_incoming_dummies) {
  #define get_a(i) (block_reverse_64(get_left_64(table_a[(i)])))
  #define get_b(i) (block_reverse_64(get_right_64(table_b[(i)], k)))
  size_t next = 0;
  size_t a_idx = 0, b_idx = 0;
  uint64_t x = 0, y = 0, x_prev = 0, y_prev = 0;

  if (num_records == 0 || num_incoming_dummies == 0) return;

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
      }
      if (++b_idx >= num_records) break;
      y_prev = y;
      y = get_b(b_idx);
    }

    // A - B: These x-nodes from table A will require incoming dummy edges
    while (a_idx < num_records && y > x) {
      // add x to result
      if (x != x_prev) {
        // make this already fetched, and x calculated from it
        uint64_t temp = table_a[a_idx];
        temp >>= 2;
        incoming_dummies[next++] = temp;
        if (next == num_incoming_dummies) return;
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
  return;
}

void prepare_k_values(unsigned char * k_values, size_t num_dummies, uint32_t k);
void prepare_k_values(unsigned char * k_values, size_t num_dummies, uint32_t k) {
  // first num_dummies are k, then it is k-1 down to 1 num_dummies times
  memset(k_values, k, num_dummies);
  unsigned char * new_k_values = k_values + num_dummies;
  // TODO: check this loop
  for (size_t i = 0; i < num_dummies * (k-1); i++) {
    new_k_values[i] = k - (i % (k-1)) - 1;
  }
}

void generate_dummies_64(uint64_t dummy_node, uint64_t * output, uint32_t k);
void generate_dummies_64(uint64_t dummy_node, uint64_t * output, uint32_t k) {
  for (size_t i = 0; i < k-1; i++) {
    // shift the width of 1 nucleotide. At this stage the kmers should be represented as reversed
    output[i] = (dummy_node >>= 2);
  }
}

void prepare_incoming_dummy_edges_64(uint64_t * dummy_nodes, unsigned char * k_values, size_t num_dummies, uint32_t k) {
  // should be k * num_dummies space allocated for dummy_nodes and k_values
  // loop over each dummy edge, and have a pointer for next, loop over k, write k and edge >> 2
  uint64_t * output = dummy_nodes + num_dummies;
  for (size_t i = 0; i < num_dummies; i++) {
    generate_dummies_64(dummy_nodes[i], output + i * (k-1), k);
  }
  prepare_k_values(k_values, num_dummies, k);
}

void merge_dummies(FILE * outfile, uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k, uint64_t * incoming_dummies, size_t num_incoming_dummies, unsigned char * dummy_lengths) {
  #define get_a(i) (block_reverse_64(get_left_64(table_a[(i)])))
  #define get_b(i) (block_reverse_64(get_right_64(table_b[(i)], k)))
  #define get_dummy(i) (block_reverse_64(incoming_dummies[(i)]))
  #define get_edge(x) (((x) & 0xC000000000000000) >> 62)
  #define get_node_suffix(d,k) (((d)<<2) >> 2*(32 - (k) + 1))
  size_t next = 0;
  size_t a_idx = 0, b_idx = 0, d_idx = 0;
  unsigned char d_len = 0, d_len_prev = 0, prev_out_len = 0;
  uint64_t x = 0, y = 0, x_prev = 0, y_prev = 0, d = 0, d_prev =0, prev_out = 0;
  char buf[k+1];
  char edge_flags[4] = {0, 0, 0, 0};
  int first=1;
  int needs_edge_flag=0;

  if (num_records == 0 && num_incoming_dummies == 0) return;

  if (num_records > 0) {
    x = get_a(a_idx);
    y = get_b(b_idx);
    // nothing special about the NOT operation here, just need a guaranteed different value to start with
    x_prev = ~x;
    y_prev = ~y;
  }

  if (num_incoming_dummies > 0) {
    d = get_dummy(d_idx);
    d_len = dummy_lengths[d_idx];
    d_prev = ~d;
  }

  while (a_idx < num_records && b_idx < num_records) {
    // B - A: These y-nodes from table B will require outgoing dummy edges
    while (b_idx < num_records && y < x) {
      // add y to result
      if (y != y_prev) {
        // Check if we need to output an incoming dummy edge
        while (d_idx < num_incoming_dummies && (d << 2) <= y) {
          if (d != d_prev || d_len != d_len_prev) {
            // print d
            first = (d<<2 != prev_out || d_len_prev != prev_out_len);
            // if the suffix is equal, the length might be different
            // but d_len = k-1
            // If change then reset flags, output 0
            int edge = get_edge(d);
            // do we need a minus flag on our edge?
            // Is this a new node suffix group?
            if (get_node_suffix(d, d_len) != get_node_suffix(prev_out >> 2,prev_out_len)
                || d_len != prev_out_len) {
              needs_edge_flag = 0;
              // These are incoming dummies so have outgoing edges that arent $
              // reset edge_flags
              memset(edge_flags, 0, 4);
              // set flag for this edge (unless $, but wont happen here)
            } else {
              needs_edge_flag = (edge_flags[edge]);
            }
            // update edge flag since we saw this edge, regardless of if we are at a new suffix group or not
            edge_flags[edge] = 1;
            sprint_dummy_acgt(buf, block_reverse_64(d), k, d_len);
            fprintf(outfile, "a %d %s %d\n", first, buf, needs_edge_flag);
            prev_out = d << 2;
            prev_out_len = d_len;
          }
          if (++d_idx >= num_incoming_dummies) break;
          d_prev = d;
          d = get_dummy(d_idx);
          d_len_prev = d_len;
          d_len = dummy_lengths[d_idx];
        }
        // Should always print this after the incoming dummies
        uint64_t temp = get_right_64(table_b[b_idx],k);
        if (get_node_suffix(temp, d_len) != get_node_suffix(prev_out >> 2,prev_out_len)
                || k != prev_out_len) {
          // reset edge_flags
          memset(edge_flags, 0, 4);
        }
        prev_out = temp << 2;
        prev_out_len = k;
        sprint_kmer_acgt(buf, &temp, k);
        buf[k-1] = '$';
        // Outgoing dummies by definition are the last edge of their node
        // And by definition we don't care about the minus flags for dummies (it is marked as a $ instead)
        fprintf(outfile, "b 1 %s 0\n", buf);
      }
      if (++b_idx >= num_records) break;
      y_prev = y;
      y = get_b(b_idx);
    }

    // A - B: These x-nodes from table A will require incoming dummy edges
    while (a_idx < num_records && y > x) {
      if (x != x_prev) {
        while (d_idx < num_incoming_dummies && (d << 2) <= x) {
          if (d != d_prev || d_len != d_len_prev) {
            // print d
            first = (d << 2 != prev_out || d_len != prev_out_len);
            int edge = get_edge(d);
            if (get_node_suffix(d, d_len) != get_node_suffix(prev_out >> 2,prev_out_len)
                || d_len != prev_out_len) {
              needs_edge_flag = 0;
              // These are incoming dummies so have outgoing edges that arent $
              // reset edge_flags
              memset(edge_flags, 0, 4);
              // set flag for this edge (unless $, but wont happen here)
            } else {
              needs_edge_flag = (edge_flags[edge]);
            }
            edge_flags[edge] = 1;
            sprint_dummy_acgt(buf, block_reverse_64(d), k, d_len);
            fprintf(outfile, "c %d %s %d\n", first, buf, needs_edge_flag);
            prev_out = d<<2;
            prev_out_len = d_len;
          }
          if (++d_idx >= num_incoming_dummies) break;
          d_prev = d;
          d = get_dummy(d_idx);
          d_len_prev = d_len;
          d_len = dummy_lengths[d_idx];
        }
        uint64_t temp = table_a[a_idx];
        first = (temp >> 2 != prev_out || k != prev_out_len);
        int edge = get_edge(temp);
        if (get_node_suffix(temp, k) != get_node_suffix(prev_out >> 2,prev_out_len)
          || k != prev_out_len) {
          needs_edge_flag = 0;
          // These are incoming dummies so have outgoing edges that arent $
          // reset edge_flags
          memset(edge_flags, 0, 4);
          // set flag for this edge (unless $, but wont happen here)
        } else {
          needs_edge_flag = (edge_flags[edge]);
        }
        edge_flags[edge] = 1;
 
        prev_out = temp <<2;
        prev_out_len = k;
        sprint_kmer_acgt(buf, &temp, k);
        fprintf(outfile, "d %d %s %d\n", first, buf, needs_edge_flag);
        //incoming_dummies[next++] = block_reverse_64(temp);
        if (next == num_incoming_dummies) return;
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
        while (d_idx < num_incoming_dummies && (d << 2) <= x) {
          if (d != d_prev || d_len != d_len_prev) {
            // print d
            first = (d << 2 != prev_out || d_len != prev_out_len);
            int edge = get_edge(d);
            if (get_node_suffix(d, d_len) != get_node_suffix(prev_out >> 2,prev_out_len)
              || d_len != prev_out_len) {
              needs_edge_flag = 0;
              // These are incoming dummies so have outgoing edges that arent $
              // reset edge_flags
              memset(edge_flags, 0, 4);
              // set flag for this edge (unless $, but wont happen here)
            } else {
              needs_edge_flag = (edge_flags[edge]);
            }
            edge_flags[edge] = 1;

            prev_out = d <<2;
            prev_out_len = d_len;
            sprint_dummy_acgt(buf, block_reverse_64(d), k, d_len);
            fprintf(outfile, "e %d %s %d\n", first, buf, needs_edge_flag);
          }
          if (++d_idx >= num_incoming_dummies) break;
          d_prev = d;
          d = get_dummy(d_idx);
          d_len_prev = d_len;
          d_len = dummy_lengths[d_idx];
        }
        uint64_t temp = table_a[a_idx];
        first = (temp << 2 != prev_out || k != prev_out_len);
        int edge = get_edge(temp);
        if (get_node_suffix(temp, k) != get_node_suffix(prev_out >> 2,prev_out_len)
          || k != prev_out_len) {
          //TRACE("k, prev_out_len = %d, %d..... ", k, prev_out_len)
          needs_edge_flag = 0;
          // These are incoming dummies so have outgoing edges that arent $
          // reset edge_flags
          memset(edge_flags, 0, 4);
          // set flag for this edge (unless $, but wont happen here)
        } else {
          needs_edge_flag = (edge_flags[edge]);
        }
        edge_flags[edge] = 1;

        prev_out = temp <<2;
        prev_out_len = k;
        sprint_kmer_acgt(buf, &temp, k);
        fprintf(outfile, "f %d %s %d\n", first, buf, needs_edge_flag);
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
  while (d_idx < num_incoming_dummies) {
    if (d != d_prev || d_len != d_len_prev) {
      // print d
      first = (d << 2 != prev_out || d_len != prev_out_len);
      int edge = get_edge(d);
      if (get_node_suffix(d, d_len) != get_node_suffix(prev_out >> 2,prev_out_len)
        || d_len != prev_out_len) {
        needs_edge_flag = 0;
        // These are incoming dummies so have outgoing edges that arent $
        // reset edge_flags
        memset(edge_flags, 0, 4);
        // set flag for this edge (unless $, but wont happen here)
      } else {
        needs_edge_flag = (edge_flags[edge]);
      }
      edge_flags[edge] = 1;

      prev_out = d <<2;
      prev_out_len = d_len;
      sprint_dummy_acgt(buf, block_reverse_64(d), k, d_len);
      fprintf(outfile, "g %d %s _\n", first, buf);
    }
    if (++d_idx >= num_incoming_dummies) break;
    d_prev = d;
    d = get_dummy(d_idx);
    d_len_prev = d_len;
    d_len = dummy_lengths[d_idx];
  }
  return;
}
