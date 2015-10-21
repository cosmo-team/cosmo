#pragma once
#ifndef DUMMIES_HPP
#define DUMMIES_HPP

#include <boost/range/adaptor/transformed.hpp>     // Map function to inputs
#include <boost/range/adaptor/uniqued.hpp>         // Uniquify
#include <boost/range/algorithm/set_algorithm.hpp> // set_difference
#include <boost/function_output_iterator.hpp>      // for capturing output of set_algorithms
#include <algorithm>
#include <utility>                                 // make_pair (for boost ranges)
#include <functional>                              // function (to avoid errors with the lambdas)
#include <cstring>                                 // memset

#include "kmer.hpp"

using namespace boost::adaptors;

enum edge_tag { standard, in_dummy, out_dummy };

template <typename kmer_t, typename InputRange1, typename InputRange2, typename Func>
void find_incoming_dummy_nodes(const InputRange1 a_range, const InputRange2 b_range, uint32_t k, Func out) {
  //typedef decltype(*a_range.begin()) kmer_t;
  //typedef typename OutputIterator::value_type pair_t;
  auto a_lam   = std::function<kmer_t(kmer_t)>([](kmer_t x) -> kmer_t {return get_start_node(x);});
  auto b_lam   = std::function<kmer_t(kmer_t)>([k](kmer_t x) -> kmer_t {return get_end_node(x,k);});
  // TODO: remove the uniqued adaptors. Kinda slow.
  auto a = a_range | transformed(a_lam) | uniqued;
  auto b = b_range | transformed(b_lam) | uniqued;

  //auto pairer  = std::function<pair_t(kmer_t)>([&](kmer_t x) -> )
  //auto pairer  = [&](kmer_t x) { *out++ = std::make_pair(x, k-1); };
  auto paired_out = boost::make_function_output_iterator(out);
  // TODO: check out std::experimental::parallel?
  // TODO: check order
  boost::set_difference(a, b, paired_out);
  //std::experimental::parallel(a.begin(), a.end(), b.begin(), b.end(), paired_out);
}

template <typename kmer_t, typename OutputIterator>
void generate_dummy_edges(const kmer_t & dummy_node, OutputIterator & output, size_t k) {
  // until k-1 because we need at least one symbol left
  kmer_t temp(dummy_node);
  for (size_t i = 0; i < k-1; i++) {
    *output++ = std::make_pair(temp <<= NT_WIDTH, k - i - 1);
  }
}

// Change to take zip iter
template <typename InputRange, typename OutputIterator1>
void prepare_incoming_dummy_edge_shifts(InputRange dummy_nodes, OutputIterator1 dummy_edges, size_t k) {//OutputIterator2 k_values, size_t k) {
  for (auto & dummy : dummy_nodes) {
    generate_dummy_edges(dummy.first, dummy_edges, k);
  }
  //prepare_k_values(k_values, num_dummies, k);
}

/*
template <typename kmer_t, typename InputRange1, typename InputRange2, typename OutputContainer>
void find_incoming_dummy_edges(const InputRange1 a_range, const InputRange2 b_range, uint32_t k, OutputContainer & incoming_dummies) {
  //typedef typename InputRange1::value_type kmer_t;
  typedef typename OutputContainer::value_type pair_t;

  // Find unique nodes that need dummy edges
  find_incoming_dummy_nodes<kmer_t>(a_range, b_range, k, std::back_inserter(incoming_dummies));

  // Generate dummy edges (all shifts prepended with $)
  //size_t num_dummy_edges = incoming_dummies.size() * (k-1); // non-unique
  //incoming_dummies.reserve(num_dummy_edges);

  prepare_incoming_dummy_edge_shifts(incoming_dummies, std::back_inserter(incoming_dummies), k-1);
}
*/

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
    }
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
// first to make parsing easy)
//
// This really needs a refactor...
// TODO: rewrite using http://www.boost.org/doc/libs/1_48_0/libs/range/doc/html/range/reference/algorithms/mutating/merge.html
template <typename InputRange1, typename InputRange2, typename InputRange3, class Visitor>
void merge_dummies(InputRange1 table_a, InputRange2 table_b,
                   InputRange3 in_dummies, size_t k, Visitor visitor_f) {
  typedef decltype(*table_a.begin()) kmer_t;

  // runtime speed: O(num_records) (since num_records >= num_incoming_dummies)
  auto visit = uniquify(add_first_start_node_flag(add_first_end_node_flag(visitor_f, k),k));

  size_t num_records = table_a.end() - table_a.begin();
  size_t num_incoming_dummies = in_dummies.end() - in_dummies.begin();

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
    kmer_t d = in_dummies[d_idx].first; \
    uint8_t len  = in_dummies[d_idx].second; \
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
    visit(in_dummy, in_dummies[d_idx].first, in_dummies[d_idx].second);
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
