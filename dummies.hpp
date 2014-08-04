#pragma once
#ifndef DUMMIES_HPP
#define DUMMIES_HPP

#include <boost/range/adaptor/transformed.hpp>     // Map function to inputs
#include <boost/range/adaptor/uniqued.hpp>         // Uniquify
#include <boost/range/algorithm/set_algorithm.hpp> // set_difference
#include <boost/function_output_iterator.hpp>      // for capturing output of set_algorithms
#include <utility>                                 // make_pair (for boost ranges)
#include <functional>                              // function (to avoid errors with the lambdas)
#include <cstring>                                 // memset

#include "kmer.hpp"

using namespace boost::adaptors;

template <typename kmer_t, typename OutputIterator>
void find_incoming_dummy_edges(const kmer_t * table_a, const kmer_t * table_b, size_t num_kmers, uint32_t k, OutputIterator out) {
  auto a_range = std::make_pair(table_a, table_a + num_kmers);
  auto b_range = std::make_pair(table_b, table_b + num_kmers);
  auto a_lam   = std::function<kmer_t(kmer_t)>([](kmer_t x) -> kmer_t {return get_start_node(x);});
  auto b_lam   = std::function<kmer_t(kmer_t)>([k](kmer_t x) -> kmer_t {return get_end_node(x,k);});
  auto a = a_range | transformed(a_lam) | uniqued;
  auto b = b_range | transformed(b_lam) | uniqued;

  //auto out_count = boost::make_function_output_iterator(inc_count);
  //boost::set_difference(a, b, out_count);
  boost::set_difference(a, b, out);
}

template <typename kmer_t>
size_t count_incoming_dummy_edges(kmer_t * table_a, kmer_t * table_b, size_t num_kmers, uint32_t k) {
  size_t count = 0;

  auto inc_count = [&count](kmer_t) {count++;};
  // This is required because set_difference requires an output iterator :<
  // Would be more self descriptive with a better pipeline lib
  auto out_count = boost::make_function_output_iterator(inc_count);
  find_incoming_dummy_edges(table_a, table_b, num_kmers, k, out_count);
  return count;
}

inline void prepare_k_values(uint8_t * k_values, size_t num_dummies, uint32_t k) {
  // first num_dummies are k, then it is k-1 down to 1 num_dummies times
  memset(k_values, k, num_dummies);
  uint8_t * new_k_values = k_values + num_dummies;
  for (size_t i = 0; i < num_dummies * (k-1); i++) {
    new_k_values[i] = k - (i % (k-1)) - 1;
  }
}

template <typename kmer_t>
void generate_dummies(kmer_t dummy_node, kmer_t * output, uint32_t k) {
  // until k-1 because we need at least one symbol left
  for (size_t i = 0; i < k-1; i++) {
    // shift the width of 1 nucleotide. (We are storing the kmers in reverse order, so a shift left
    // is a shift right when printed out)
    // Yeah, maybe this should be in kmer's code and produce an iterator instead...
    // but who has time to well-abstracted code when papers need publishing?
    output[i] = dummy_node <<= NT_WIDTH;
  }
}

template <typename kmer_t>
void prepare_incoming_dummy_edges(kmer_t * dummy_nodes, uint8_t * k_values, size_t num_dummies, uint32_t k) {
  // should be k * num_dummies space allocated for dummy_nodes and k_values
  // loop over each dummy edge, and have a pointer for next, loop over k, write k and edge >> 2
  kmer_t * output = dummy_nodes + num_dummies;
  for (size_t i = 0; i < num_dummies; i++) {
    generate_dummies(dummy_nodes[i], output + i * (k-1), k);
  }
  prepare_k_values(k_values, num_dummies, k);
}

// Could be done cleaner: set_difference iterator as outgoing dummies, transform to have tuple with k value, merge + merge again iterator with comp functor.
//
// TODO: set_difference_iterator, tested against set_difference algorithm
// TODO: find_outgoing_dummy_edges iterator
// TODO: merge_iterator, tested with merge algorithm
//
// but I think an extra set_difference indirection might make things slower (since we have to access outgoing dummies multiple times),
// and this was already written and tested from an earlier C implementation.
// Visitor functor takes 4 params: kmer, size, first flag, edge flag (could also just take kmer and size)
// planned Visitor functors: ascii_full_edge, ascii_edge_only, binary (5 bits per row, x12 per 64 bit block, 4 bits waste per 12, or just per 8 bits at
// first to make parsing easy)
/*
template <typename kmer_t, typename Visitor>
void merge_dummies(kmer_t * table_a, kmer_t * table_b, size_t num_records, uint32_t k,
                   kmer_t * incoming_dummies, kmer_t num_incoming_dummies, uint8_t * dummy_lengths,
                   Visitor visit) {
  #define get_a(i) (block_reverse_64(get_left_64(table_a[(i)])))
  #define get_b(i) (block_reverse_64(get_right_64(table_b[(i)], k)))
  #define get_dummy(i) (block_reverse_64(incoming_dummies[(i)]))
  #define get_edge(x) (((x) & 0xC000000000000000) >> 62)
  #define get_node_suffix(d) (((d)<<2) & ~(0xC000000000000000 >> (k-2)*2))
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
            first = (d<<2 != prev_out || d_len != prev_out_len);
            // if the suffix is equal, the length might be different
            // but d_len = k-1
            // If change then reset flags, output 0
            int edge = get_edge(d);
            // do we need a minus flag on our edge?
            // Is this a new node suffix group?
            if (get_node_suffix(d) != get_node_suffix(prev_out >> 2)
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
            fprintf(outfile, "%d %s %d\n", first, buf, needs_edge_flag);
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
        if (get_node_suffix(block_reverse_64(temp)) != get_node_suffix(prev_out >> 2)
                || k != prev_out_len) {
          // reset edge_flags
          memset(edge_flags, 0, 4);
        }
        prev_out = block_reverse_64(temp) << 2;
        prev_out_len = k;
        sprint_kmer_acgt(buf, &temp, k);
        buf[k-1] = '$';
        // Outgoing dummies by definition are the last edge of their node
        // And by definition we don't care about the minus flags for dummies (it is marked as a $ instead)
        fprintf(outfile, "1 %s 0\n", buf);
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
            if (get_node_suffix(d) != get_node_suffix(prev_out >> 2)
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
            fprintf(outfile, "%d %s %d\n", first, buf, needs_edge_flag);
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
        first = (block_reverse_64(temp) << 2 != prev_out || k != prev_out_len);
        int edge = get_edge(block_reverse_64(temp));
        if (get_node_suffix(block_reverse_64(temp)) != get_node_suffix(prev_out >> 2)
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
 
        prev_out = block_reverse_64(temp) <<2;
        prev_out_len = k;
        sprint_kmer_acgt(buf, &temp, k);
        fprintf(outfile, "%d %s %d\n", first, buf, needs_edge_flag);
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
            if (get_node_suffix(d) != get_node_suffix(prev_out >> 2)
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
            fprintf(outfile, "%d %s %d\n", first, buf, needs_edge_flag);
          }
          if (++d_idx >= num_incoming_dummies) break;
          d_prev = d;
          d = get_dummy(d_idx);
          d_len_prev = d_len;
          d_len = dummy_lengths[d_idx];
        }
        uint64_t temp = table_a[a_idx];
        first = (block_reverse_64(temp) << 2 != prev_out || k != prev_out_len);
        int edge = get_edge(block_reverse_64(temp));
        if (get_node_suffix(block_reverse_64(temp)) != get_node_suffix(prev_out >> 2)
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

        prev_out = block_reverse_64(temp) <<2;
        prev_out_len = k;
        sprint_kmer_acgt(buf, &temp, k);
        fprintf(outfile, "%d %s %d\n", first, buf, needs_edge_flag);
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
      if (get_node_suffix(d) != get_node_suffix(prev_out >> 2)
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
      fprintf(outfile, "%d %s %d\n", first, buf, needs_edge_flag);
    }
    if (++d_idx >= num_incoming_dummies) break;
    d_prev = d;
    d = get_dummy(d_idx);
    d_len_prev = d_len;
    d_len = dummy_lengths[d_idx];
  }
  return;
}
*/

#endif
