#ifndef DUMMIES_HPP
#define DUMMIES_HPP

#include <boost/range/adaptor/transformed.hpp>     // Map function to inputs
#include <boost/range/adaptor/uniqued.hpp>         // Uniquify
#include <boost/range/algorithm/set_algorithm.hpp> // set_difference
#include <boost/function_output_iterator.hpp>      // for capturing output of set_algorithms
#include <utility>                                 // make_pair (for boost ranges)
#include <functional>                              // function (to avoid errors)
#include <cstring>                                 // memset

#include "kmer.hpp"

#define BOOST_RESULT_OF_USE_DECLTYPE

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

/*
// Write own merging iterator (to output)
// This should go in io.hpp, but we need a function to get regular edges and dummy-out edges iterator (with type and k)
// and a function to provide an iterator for dummy in edges (with type and k)
// OUTPUT must be UNIQUED to avoid COPALINDROMIC edges (i.e. rev_comp(x) maps to an element already in the set of kmers... only to itself if even k, because of our original set)
void merge_dummies(FILE * outfile, uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k, uint64_t * incoming_dummies, size_t num_incoming_dummies, unsigned char * dummy_lengths);
*/

#endif
