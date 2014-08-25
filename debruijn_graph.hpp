#pragma once
#ifndef _DEBRUIJN_GRAPH_H
#define _DEBRUIJN_GRAPH_H

#include <algorithm>
#include <fstream>
#include <vector>
#include <array>
#include <string>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debug.h"

using namespace std;
using namespace sdsl;

// TODO: convert asserts into exceptions? (Copy Boost)
template <size_t t_sigma            = 4, // default: DNA, TODO: change to 0 and make dynamic
          class  t_bit_vector_type  = sd_vector<>,
          class  t_bv_rank_type     = typename t_bit_vector_type::rank_0_type, // We invert the bits so it is 0 instead
          class  t_bv_select_type   = typename t_bit_vector_type::select_0_type,
          class  t_edge_vector_type = wt_huff<rrr_vector<67>>,
          class  t_symbol_type      = typename t_edge_vector_type::value_type,
          class  t_label_type       = string> // can define basic_string<t_symbol_type>, but need to use correct print func
class debruijn_graph {
  static_assert(t_sigma == 4, "Alphabet sizes other than 4 are not yet supported.");

  public:
  const static size_t sigma = t_sigma;
  typedef t_symbol_type symbol_type;
  typedef t_label_type  label_type;

  const size_t           k{};

  private:
  const t_bit_vector_type        m_node_flags{};
  const t_bv_rank_type           m_node_rank{};
  const t_bv_select_type         m_node_select{};
  const t_edge_vector_type       m_edges{};
  // This is the "F table" in the blog/paper. It stores the starting positions of the sorted runs
  // of the k-1th symbol of the edge (or, last symbol of the node)
  // Could be implemented as the start positions, but here we store the cumulative sum (i.e. run ends)
  const array<size_t, 1+sigma> m_symbol_ends{};
  const label_type             m_alphabet{};

  private:
  debruijn_graph(size_t k, const t_bit_vector_type & node_flags, const t_edge_vector_type & edges, const array<size_t, 1+sigma>& symbol_ends, const label_type& alphabet)
    : k(k), m_node_flags(node_flags), m_node_rank(&m_node_flags), m_node_select(&m_node_flags), m_edges(edges), m_symbol_ends(symbol_ends), m_alphabet(alphabet){
  }

  public:
  static debruijn_graph load_from_packed_edges(istream & input, label_type alphabet=label_type{}) {
    // ifstream input(filename, ios::in|ios::binary|ios::ate);
    // check length
    streampos size = input.tellg();
    // should be exceptions...
    assert(size/sizeof(uint64_t) >= (sigma+2)); // space for footer info
    assert(((size_t)size - (sigma+2) * sizeof(uint64_t))%sizeof(uint64_t) == 0); // sequence of uint64_ts

    // read footer
    input.seekg(-(sigma+2) * sizeof(uint64_t), ios::end);
    array<size_t,1+sigma> counts{};
    uint64_t k = 0;
    input.read((char*)&counts[0], (sigma+1) * sizeof(uint64_t));
    input.read((char*)&k, sizeof(uint64_t));
    size_t num_edges = counts[sigma];

    size_t num_blocks = size_t(size)/sizeof(uint64_t) - (sigma+2);
    input.seekg(0, ios::beg); // rewind

    // TODO: Loop through this while reading from file instead (to avoid memory waste)
    vector<uint64_t> blocks(num_blocks,0);
    input.read((char*)&blocks[0], sizeof(uint64_t) * num_blocks);

    int_vector<1> first(num_edges,0);
    // would be nice to fix wavelet trees so the constructor
    // can accept a int_vector<4> instead (which is all we need for DNA)
    int_vector<8> edges(num_edges);

    for (size_t i = 0; i < num_edges; i++) {
      auto x = get_edge(blocks.begin(), i);
      first[i] = 1-get<1>(x); // convert 0s to 1s so we can have a sparse bit vector
      // For branchy graphs it might be better to change this and use RRR
      edges[i] = (get<0>(x) << 1) | !get<2>(x);
    }

    t_bit_vector_type bv(first);
    t_edge_vector_type wt;
    construct_im(wt, edges);
    return debruijn_graph(k, bv, wt, counts, alphabet);
  }

  // Loaders/Writers for sdsl-serialized and JSON
  // size_in_bytes:

  // API
  /*
  size_t outdegree(size_t v) {
    auto range = _node_range(v);
    size_t first = get<0>(range);
    size_t last  = get<1>(range);
    return last - first + 1;
  }
  */
  // outgoing
  // indegree
  // incoming
  // successors
  // predecessors
  /*
  label_type node_label(size_t v) {
  }
  */
  // edge_label <- index
  // node <- string
  size_t num_edges() const { return m_symbol_ends[sigma]; /*_node_flags.size();*/ }
  size_t num_nodes() const { return m_node_rank(num_edges()); }

  private:
  size_t _symbol_start(symbol_type x) {
    assert(x < sigma + 1);
    return (x==0)? 0 : m_symbol_ends[x - 1];
  }

  symbol_type _strip_edge_flag(symbol_type x) {
    return x >> 1;
  }

  symbol_type _with_edge_flag(symbol_type x, bool edge_flag = true) {
    return (x << 1) | edge_flag;
  }

  symbol_type _symbol_access(size_t i) {
    assert(i < num_edges());
    // I assume std binary search is optimised for small ranges (if sigma is small, e.g. DNA)
    return upper_bound(m_symbol_ends.begin(), m_symbol_ends.end(), i) - m_symbol_ends.begin();
  }

  /*
  size_t _succ(size_t i) {
    // TODO: bounds checkS (return value shouldn't be too far out either)
    return _node_select(_node_rank(i) + 1);
  }

  size_t _pred(size_t i) {
    return _node_select(_node_rank(i) - 1);
  }
  */

  size_t _rank_distance(size_t a, size_t b) {
    TRACE("l[%zu] = %llu, r[%zu] = %llu\n", b, m_node_flags[b], b, m_node_rank(b));
    TRACE("l[%zu] = %llu, r[%zu] = %llu\n", a, m_node_flags[a], a, m_node_rank(a));
    return m_node_rank(b) - m_node_rank(a);
  }

  /*
  // Return index of first possible edge obtained by following edge i
  size_t _forward(size_t i) {
    assert(i < num_edges());
    symbol_t x = _edges[i] >> 1; // shift to remove flag
    // TODO: if x == 0 ($) then we can't follow the edge, so exception?
    size_t start = _symbol_start(x);
    return _succ(start);
  }
  */

  public:
  size_t _backward(size_t i) {
    assert(i < num_edges());
    symbol_type x  = _symbol_access(i);
    // This handles x = $ so that we have the all-$ edge at position 0
    // but probably shouldn't be called this way in the normal case
    size_t x_start = _symbol_start(x);
    TRACE("x_start  = %zu\n", x_start);
    // rank is over [0,i) and select is 1-based
    size_t nth = _rank_distance(x_start+1, i+1);
    if (x == 0) return 0;
    // no minus flag because we want the FIRST
    return m_edges.select(nth+1, _with_edge_flag(x, false));
  }

  symbol_type _map_symbol(symbol_type x) {
    return (m_alphabet.size() > 0)? m_alphabet[x] : x;
  }

  size_t _first_edge(size_t v) {
    assert(v < num_nodes());
    // select is 1-based, but nodes are 0-based
    return m_node_select(v+1);
  }

  size_t _last_edge(size_t v) {
    // find the *next* node's first edge and decrement
    // as long as a next node exists!
    assert(v + 1 <= num_nodes());
    if (v+1 == num_nodes()) return num_edges() - 1;
    else return _first_edge(v+1) - 1;
  }

  pair<size_t, size_t> _node_range(size_t v) {
    return make_pair(_first_edge(v), _last_edge(v));
  }
};



#endif
