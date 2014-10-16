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

enum edge_tag { standard, in_dummy, out_dummy };

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

template <class Visitor>
class Unique {
  Visitor _v;
  bool first_iter = true;
  edge_tag last_tag;
  uint32_t last_k;

  public:
    Unique(Visitor v) : _v(v) {}

    template <typename kmer_t>
    void operator()(edge_tag tag, const kmer_t & x, uint32_t k) {
      static kmer_t last_kmer;

      if (first_iter || tag != last_tag || x != last_kmer || k != last_k) {
        _v(tag, x, k);
      }
      first_iter = false;
      last_tag = tag;
      last_kmer = x;
      last_k = k;
    };
};

template <class Visitor>
class FirstStartNodeFlagger {
  Visitor _v;
  uint32_t _graph_k;

  public:
    FirstStartNodeFlagger(Visitor v, uint32_t k) : _v(v), _graph_k(k) {}
    template <typename kmer_t>
    void operator()(edge_tag tag, const kmer_t & x, const uint32_t k) {
      // Note: for multithreading, this might be dangerous? could make class members instead
      static uint32_t last_k = 0;

#ifdef VAR_ORDER
      static kmer_t last_edge = 0;
      size_t result = node_lcs(x, last_edge, std::min(k, last_k));
      // If dummy edge length is equal, and the LCS length is the length of the non-dummy suffix,
      // include the $ signs as well
      if (last_k == k && k - 1 == result)
        result = _graph_k - 1;
      last_edge = x;
      // This should bump to K when it is == this_k - 1
      //std::cerr << "kmer: " << kmer_to_string(x,k) << std::endl;
      //std::cerr << "lcs: " << result << std::endl;
#else
      static kmer_t last_start_node = 0;
      static bool first_iter = true;
      kmer_t this_start_node = get_start_node(x);
      size_t result = (first_iter || this_start_node != last_start_node || k != last_k);
      first_iter = false;
      last_start_node  = this_start_node;
#endif

      _v(tag, x, k, result);
      last_k     = k;
    }
};

template <class Visitor>
class FirstEndNodeFlagger {
  Visitor  _v;
  uint32_t _graph_k;
  bool first_iter = true;

  public:
    FirstEndNodeFlagger(Visitor v, uint32_t k) : _v(v), _graph_k(k) {}
    template <typename kmer_t>
    void operator()(edge_tag tag, const kmer_t & x, const uint32_t k, uint32_t start_node_flag) {
      static kmer_t last_suffix;
      static uint32_t last_k;
      static bool edge_seen[DNA_RADIX];
      #define reset_flags() memset(edge_seen, 0, DNA_RADIX)

      bool edge_flag = true;
      kmer_t this_suffix = get_start_node_suffix(x, _graph_k);

      // reset "edge seen" flags
      if (this_suffix != last_suffix || k != last_k || first_iter) {
        reset_flags();
        edge_flag = true;
        first_iter = false;
      }
      // Only get the edge label if not a dummy out edge
      if (tag != out_dummy){
        uint8_t edge = get_edge_label(x);
        edge_flag = !edge_seen[edge];
        edge_seen[edge] = true;
      }

      _v(tag, x, k, start_node_flag, edge_flag);

      last_suffix = this_suffix;
      last_k      = k;
    }
};

template <class Visitor>
auto uniquify(Visitor v) -> Unique<decltype(v)> {
  return Unique<decltype(v)>(v);
}

template <class Visitor>
auto add_first_start_node_flag(Visitor v, uint32_t k) -> FirstStartNodeFlagger<decltype(v)> {
  return FirstStartNodeFlagger<decltype(v)>(v, k);
}

template <class Visitor>
auto add_first_end_node_flag(Visitor v, uint32_t k) -> FirstEndNodeFlagger<decltype(v)> {
  return FirstEndNodeFlagger<decltype(v)>(v, k);
}

// Could be done cleaner: set_difference iterator as outgoing dummies, transform to have tuple with k value, merge + merge again iterator with comp functor.
// but I think an extra set_difference indirection might make things slower (since we have to access outgoing dummies multiple times),
// and this was already written and tested from an earlier C implementation.
// Visitor functor takes 4 params: kmer, size, first flag, edge flag (could also just take kmer and size)
// planned Visitor functors: ascii_full_edge, ascii_edge_only, binary (5 bits per row, x12 per 64 bit block, 4 bits waste per 12, or just per 8 bits at
// first to make parsing easy)
template <typename kmer_t, class Visitor>
void merge_dummies(kmer_t * table_a, kmer_t * table_b, const size_t num_records, const uint32_t k,
                   kmer_t * in_dummies, size_t num_incoming_dummies, uint8_t * dummy_lengths,
                   Visitor visitor_f) {
  // runtime speed: O(num_records) (since num_records >= num_incoming_dummies)
  auto visit = uniquify(add_first_start_node_flag(add_first_end_node_flag(visitor_f, k),k));

  #define get_a(i) (get_start_node(table_a[(i)]) >> 2)
  #define get_b(i) (get_end_node(table_b[(i)], k) >> 2) // shifting to give dummy check call consistency
  #define inc_b() while (b_idx < num_records && get_b(++b_idx) == b) {}

  // **Standard edges**: Table a (already sorted by colex(node), then edge).
  // Table a May not be unique (if k is odd and had "palindromic" [in DNA sense] kmer in input)
  // **Outgoing dummies**: {get_b} - {get_a} (the difference of the sets constructed by mapping get_a and get_b)
  // (in the correct order, as B is sorted by last colex(row)).
  // **Incoming Dummies**: in_dummies calculated previously {get_a} - {get_b}, sorted after adding shifts to produce all required $-prefixes
  // e.g. {$acgt, $ccag} -> {$$$$a, $$$$a, $$$ac, $ccag, $$$cc, $$acg, $$cca, $acgt } [note $$$$a is repeated -> not unique]
  size_t a_idx = 0, b_idx = 0, d_idx = 0;

  // Ew macros! I know, I know...
  // (d<=s) because dummies should always sort before anything that is equal to them
  // << 2 to compare node instead. Don't need to compare edge since already sorted
  #define check_for_in_dummies(s) while (d_idx < num_incoming_dummies){ \
    kmer_t d = in_dummies[d_idx]; \
    uint8_t len  = dummy_lengths[d_idx]; \
    if (d<<2 <= ((s)<<2)) { visit(in_dummy, d, len); ++d_idx; }\
    else break; \
  }

  // set_difference loop with mixed-in merge (for in_dummies)
  // visit each a in A, otherwise B if b < a and b not in A (so A U B-A)
  // then print all remaining if either one is depleted
  // at each print, visit all in_dummies < this
  // visit(standard, table_a[a_idx++], k);
  while (a_idx < num_records && b_idx < num_records) {
    kmer_t x = table_a[a_idx];
    kmer_t a = get_a(a_idx);
    kmer_t b = get_b(b_idx);
    // B - A -> dummy out
    if (b < a) {
      check_for_in_dummies(b);
      visit(out_dummy, b, k);
      inc_b();
    }
    // incoming dummy detection, but not the correct position, so just inc a-ptr
    else if (a < b) {
      check_for_in_dummies(x);
      visit(standard, x, k);
      ++a_idx;
    }
    else {
      check_for_in_dummies(x);
      visit(standard, x, k);
      ++a_idx;
      inc_b();
    }
  }

  // Might have entries in a even if b is depleted (e.g. if all b < a)
  while (a_idx < num_records) {
    kmer_t x = table_a[a_idx++];
    check_for_in_dummies(x);
    visit(standard, x, k);
  }

  // Might have entries in b even if a is depleted
  while (b_idx < num_records) {
    kmer_t b = get_b(b_idx++);
    check_for_in_dummies(b);
    visit(out_dummy, b, k);
  }

  // Might have in-dummies remaining
  while (d_idx < num_incoming_dummies) {
    visit(in_dummy, in_dummies[d_idx], dummy_lengths[d_idx]);
    ++d_idx;
  }
}

template <typename kmer_t>
uint8_t get_w(edge_tag tag, const kmer_t & x) {
  if (tag == out_dummy) return 0;
  else return get_edge_label(x) + 1;
}

template <typename kmer_t>
uint8_t get_f(edge_tag tag, const kmer_t & x, const uint32_t k) {
  uint8_t sym;
  if (tag == in_dummy && k == 1) {
    return 0; // in_dummies might have $ if only an edge label
  }
  else if (tag == out_dummy) {
    sym = get_nt(x, 1); // since out_dummies are shifted
  }
  else {
    sym = get_nt(x, 1); // 2nd last symbol
  }
  return sym+1;
}


#endif
