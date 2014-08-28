#pragma once
#ifndef _ALGORITHM_HPP
#define _ALGORITHM_HPP

#include <sdsl/bit_vectors.hpp>
#include "debruijn_graph.hpp"

using namespace sdsl;

// Calling it with a vector of first minus positions is sliiiightly faster
// But would mean we have to save another vector with our graph if we ever want to traverse
// it again.
// If a vector isnt provided, another method is used which is *almost as fast*
template <bool has_branch = 1, class V=sd_vector<>, class G>
V make_branch_vector(const G& g, const vector<size_t> * v = nullptr) {
  bit_vector branching_flags(g.num_edges(),~has_branch);
  // INDEGREE
  for (size_t edge = 0; edge < g.num_edges()-1; edge++) {
    // This will only unset the bits that are 0,0 (which means no 11, 10, or 01, which is a sibling edge)
    branching_flags[edge] = (g.m_node_flags[edge] || g.m_node_flags[edge+1]) ^ has_branch;
  }
  // The last edge is whatever it is in the graph
  branching_flags[g.num_edges()-1] = g.m_node_flags[g.num_edges()-1] ^ has_branch;

  // OUTDEGREE
  if (v) {
    for (auto edge : *v) {
      ssize_t next = g._forward(edge);
      branching_flags[next] = has_branch;
    }
  }
  else {
    // This is close to the same speed as loading in vector created in the constructor
    // But easier for people to use
    typedef typename G::symbol_type symbol_type;
    for (symbol_type sym = 1; sym < g.sigma + 1; sym++) {
      symbol_type x = g._with_edge_flag(sym, false);
      symbol_type x_minus = g._with_edge_flag(sym, true);
      size_t bound = g.m_edges.rank(g.num_edges(), x_minus);
      for(size_t next_sel = 1; next_sel < bound; ) {
        size_t edge = g.m_edges.select(next_sel, x_minus);
        branching_flags[g._forward(edge)] = has_branch;
        next_sel = 1 + g.m_edges.rank(g.m_edges.select(1+g.m_edges.rank(edge, x),x), x_minus);
      }
    }
  }
  return V(branching_flags);
}

/*
// Should make iterator instead?
template <class G, class F>
void visit_unipaths(const G& g, const F & f) {
  sd_vector<> b = make_branch_vector(g);
  cout << size_in_bytes(b) << endl;
  f(' ');
  //for (auto x : b) {
  //  cout << x << endl;
  //}
}
*/

#endif
