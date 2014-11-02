#pragma once
#ifndef _DEBRUIJN_HYPERGRAPH_H
#define _DEBRUIJN_HYPERGRAPH_H

#include <vector>
#include <algorithm>
#include <boost/optional.hpp>
//#include <sdsl/wt_algorithm.hpp>
#include "wt_algorithm.hpp"
#include "debruijn_graph.hpp"

using namespace boost;

template <class t_debruijn_graph = debruijn_graph<>,
          class t_lcs_vector     = wt_int<rrr_vector<63>> >
class debruijn_hypergraph : t_debruijn_graph {
  // Nodes are now ranges of edges (could be done for standard de bruijn graph too...
  // would save time by saving state?)
  // A node is a contiguous range of edges, a hypernode is a contiguous range of nodes (hence a range of edges)
  // A hypernode then replaces the defn of node without loss of generality
  public:
  typedef size_t                     edge_type;
  typedef tuple<edge_type, edge_type, size_t> node_type;
  //typedef optional<node_type> hypernode_type; // optional return type

  typedef typename t_debruijn_graph::symbol_type symbol_type;
  const t_debruijn_graph & m_dbg;
  const t_lcs_vector       m_lcs;

  debruijn_hypergraph(const t_debruijn_graph & dbg, const t_lcs_vector & lcs) : m_dbg(dbg), m_lcs(lcs) {}

  // shorter(v, k) - returns the hypernode whose label is the last k characters of v's label (reduce context)
  node_type shorter(const node_type & v, size_t k) {
    if (k <= 1) return node_type(0, m_dbg.num_edges()-1, k);

    // search backward on WT to find the first occurence of a number less than k
    size_t i = get<0>(v);
    size_t j = get<1>(v);

    // find largest i' <= i with L*[i' - 1] < k
    size_t i_prime = prev_lte(m_lcs, i+1, k-1)-1;
    // find smallest j' >= j with L*[j'] < k
    size_t j_prime = next_lte(m_lcs, j+2, k-1)-2;

    return node_type(i_prime, std::max(i_prime, j_prime), k);
  }

  // longer(v, k) - list nodes (new "node") whose labels have length k <= K and end with v's label
  // could make iterator instead
  vector<node_type> longer(const node_type & v, size_t k) {
    vector<node_type> longer_nodes;
    size_t i = get<0>(v);
    size_t j = get<1>(v);
    auto starts = range_lte(m_lcs, i, j, k-1);
    // add code to add first edge?
    for (size_t idx = 1; idx < starts.size()-1; idx++) {
      size_t start = starts[idx-1];
      size_t end   = starts[idx]-1;
      longer_nodes.push_back(node_type(start, end, k));
    }
    size_t last = starts[starts.size()-1];
    longer_nodes.push_back(node_type(last-1, j, k));
    return longer_nodes;
  }

  symbol_type lastchar(const node_type& v) const {
    return m_dbg._symbol_access(get<0>(v));
  }

  optional<node_type> outgoing(const node_type & v, symbol_type x) {
    assert(x < m_dbg.sigma + 1);
    if (x == 0) return optional<node_type>();
    auto max = maxlen(v, x);
    if (!max) return optional<node_type>();
    size_t i = m_dbg._outgoing_edge_pair(get<0>(*max), get<1>(*max), x);
    size_t j = m_dbg._last_edge_of_node(m_dbg._edge_to_node(i));
    return optional<node_type>(shorter(node_type(i,j,m_dbg.k-1), get<2>(v)));
  }

  vector<node_type> backward(const node_type & v) {
    // This could be done lazily, and searched over...
    //if (get<2>(v) == m_dbg.k-1) {} // do standard version - not needed?
    auto l = longer(v, get<2>(v)+1);
    for (auto u : l) { cout << "1. " << get<0>(u) << ", " << get<1>(u) << endl; } cout<<endl;
    // map maxlen to each element in l
    transform(l.begin(), l.end(), l.begin(), [&](const node_type & u){ return maxlen(u); } );
    for (auto u : l) { cout << "2. " << get<0>(u) << ", " << get<1>(u) << endl; }cout<<endl;
    // map standard dbg backward to each element
    transform(l.begin(), l.end(), l.begin(), [&](const node_type & u){
      size_t start = m_dbg._backward(get<0>(u));
      size_t end   = m_dbg._last_edge_of_node(m_dbg._edge_to_node(start));
      return node_type(start, end, get<2>(u));
    });
    for (auto u : l) { cout << "3. " << get<0>(u) << ", " << get<1>(u) << endl; }cout<<endl;
    // map shorter to each element
    transform(l.begin(), l.end(), l.begin(), [&](const node_type & u){ return shorter(u, get<2>(v)); });
    for (auto u : l) { cout << "4. " << get<0>(u) << ", " << get<1>(u) << endl; }cout<<endl;
    return l;
  }

  // maxlen(v, x) - returns some node in the *original* (kmax) graph whose label ends with v's
  // label, and that has an outgoing edge labelled x, or NULL otherwise
  // should return dbg::node_type
  node_type maxlen(const node_type & v) const {
    size_t start = get<0>(v); // have to update this to select to the next node
    // Range_start must also be a start of a node at the top level context, so find the next top level node to find the range
    // TODO: make (or check if it exists) uniform interface for rank/select, so I can easily make next(), prev(), etc
    size_t end = m_dbg._last_edge_of_node(m_dbg._edge_to_node(start));
    return node_type(start, end, m_dbg.k-1);
  }

  optional<node_type> maxlen(const node_type & v, const symbol_type x) const {
    assert(x < m_dbg.sigma + 1);
    // For both flagged and nonflagged symbol in W
    for (symbol_type flag : {0,1}) {
      symbol_type x_with_flag = m_dbg._with_edge_flag(x, flag);
      size_t prev_count = m_dbg.m_edges.rank(get<0>(v), x_with_flag);
      size_t next = m_dbg.m_edges.select(prev_count+1, x_with_flag);
      // check if edge falls outside our range...
      if (next > get<1>(v)) break;
      // Find node range
      size_t node_rank = m_dbg._edge_to_node(next);
      auto n_range = m_dbg._node_range(node_rank);
      return optional<node_type>(node_type(get<0>(n_range), get<1>(n_range), m_dbg.k-1));
    }
    return optional<node_type>();
  }

  // function to get standard node by rank (then use shorter)
  // fwd, back, lastchar
};

#endif
