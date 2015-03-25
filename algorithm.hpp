#pragma once
#ifndef _ALGORITHM_HPP
#define _ALGORITHM_HPP

#include <vector>
#include <stack>

#include <sdsl/bit_vectors.hpp>
#include "debruijn_graph.hpp"
#include "kmer.hpp"

using namespace std;
//using namespace sdsl;
/*
// Calling it with a vector of first minus positions is sliiiightly faster
// But would mean we have to save another vector with our graph if we ever want to traverse
// it again.
// If a vector isnt provided, another method is used which is *almost as fast*
template <bool has_branch = 1, class G, class V=sdsl::sd_vector<> >
V make_branch_vector(const G& g, const vector<size_t> * v = nullptr) {
  sdsl::bit_vector branching_flags(g.num_edges(),!has_branch);
  // OUTDEGREE
  for (size_t edge = 0; edge < g.num_edges()-1; edge++) {
    // This will only unset the bits that are 0,0 (which means no 11, 10, or 01, which is a sibling edge)
    // 0,0 -> !has_branch, 0,1 -> has_branch,
    // 1,0 ->  has_branch, 1,1 -> has_branch,
    branching_flags[edge] = ((g.m_node_flags[edge] || g.m_node_flags[edge+1]) == has_branch);
  }
  // The last edge is whatever it is in the graph
  // 0 -> !has_branch, 1-> has_branch
  branching_flags[g.num_edges()-1] = (g.m_node_flags[g.num_edges()-1] == has_branch);

  // INDEGREE
  // Seems like we should iterate over all possible edges, but these will have been captured
  // by the outdegree loop above
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
    // Ignore $ edges, which shouldn't be there anyway, but follow every other symbol
    for (symbol_type sym = 1; sym < g.sigma + 1; sym++) {
      symbol_type x = g._with_edge_flag(sym, false);      // edges without flag
      symbol_type x_minus = g._with_edge_flag(sym, true); // edges with flag
      size_t non_minus_bound = g.m_edges.rank(g.num_edges(), x); // how many edges
      size_t minus_bound = g.m_edges.rank(g.num_edges(), x_minus); // how many minus flags do we have?
      for(size_t next_sel = 1; next_sel <= minus_bound; ) { // bound is inclusive
        size_t edge = g.m_edges.select(next_sel, x_minus); // locate the appropriate pre-edge (must have minus otherwise outdegree <= 1)
        branching_flags[g._forward(edge)] = has_branch; // follow edge and set the flag
        size_t prev_non_flags = g.m_edges.rank(edge, x); // edge is a minus
        size_t next_non_flag = prev_non_flags + 1; // Check that this isn't out of bounds
        size_t next_non_flag_edge = (next_non_flag > non_minus_bound)? g.num_edges() - 1 : g.m_edges.select(next_non_flag, x);
        // Count flags before and select to next one
        if (next_sel == minus_bound) break;
        next_sel = 1 + g.m_edges.rank(next_non_flag_edge, x_minus);
      }
    }
  }
  return V(branching_flags);
}

template <bool has_branch, class V>
struct _type_getter{ };

template <class V> struct _type_getter<0, V> {
  typedef typename V::rank_0_type rank_type;
  typedef typename V::select_0_type select_type;
};
template <class V> struct _type_getter<1, V> {
  typedef typename V::rank_1_type rank_type;
  typedef typename V::select_1_type select_type;
};


// Not technically a visitor... but a visitor pattern
// Should change this to an edge iterator sometime
// Or refactor it to unipath_debruijn_graph
template <class G, const bool has_branch = 1, class B=sdsl::sd_vector<>>
class unipath_visitor {
  typedef sdsl::bit_vector bit_vector;
  typedef typename _type_getter<has_branch, B>::rank_type rank_type;
  typedef typename _type_getter<has_branch, B>::select_type select_type;
  typedef typename G::symbol_type symbol_type;
  typedef typename G::label_type  label_type;

  public: // change this later
  const G & g;
  const B   b;
  const rank_type   branch_rank;
  const select_type branch_select;

  // CTORs
  unipath_visitor(const G& dbg) : g(dbg), b(make_branch_vector(g)), branch_rank(&b), branch_select(&b) {}
  // The below may be slightly faster (if a branch vector is prebuilt at graph construction)
  unipath_visitor(const G& dbg, B & branches) : g(dbg), b(branches), branch_rank(&b), branch_select(&b) {}
  unipath_visitor(const G& dbg, const vector<size_t> & v) : g(dbg), b(make_branch_vector(g, v)), branch_rank(&b), branch_select(&b) {}

  template <class F>
  void operator()(const F & f) const {
    // Need to keep track of visited incase we have isolated cycles (e.g. plasmids)
    bit_vector visited(g.num_edges(),0);
    // We increment this for a faster popcount at the end (to check if we need to traverse cycles)
    size_t     num_visited = 0;

    _traverse_branches(visited, num_visited, f);
    _traverse_cycles(visited, num_visited, f);
  }

  template <class F>
  void _traverse_branches(bit_vector & visited, size_t &, const F &) const {
    // This way was just to avoid a big stack...
    // Select to each contig split
    for (size_t sel_idx = 1; sel_idx <= g.num_edges(); sel_idx++) {
      size_t base = branch_select(sel_idx);
      size_t last = g._last_sibling(base); // first not needed in this calculation
      size_t num_siblings = last - base + 1;

      // Check siblings to see if already visited (then we can skip this)
      size_t non_visited_count = 0;
      for (size_t sibling = 0; sibling < num_siblings; sibling++) {
        if (!visited[base+sibling]) non_visited_count++;
      }
      // Skip the siblings if all are visited
      if (non_visited_count > 0) {
        // Loop backwards to get starting k-1 (if $ before that then we skip... and only do this if the siblings arent visited)

        // loop across sibling branches to avoid selecting (access is faster) and backtracking multiple times
        for (size_t sibling = 0; sibling < num_siblings; sibling++) {
          // If an element is already visited, just skip it (but keep iterating)
          if (visited[base+sibling]) continue;
          // do unipath path traversal
        }
      }
      sel_idx += num_siblings;
    }
  }

  template <class F>
  void _traverse_cycles(bit_vector &, size_t num_visited, const F &) const {
    // If we have visited every edge, then we have nothing left to do!
    // Anything else is a cycle (if we handled the branches before)
    if (num_visited == g.num_edges()) return;

    // Just iterate over all edges that arent visited, then follow them until we get to a visited
  }

  // Assumes starting from edge straight after branch == has_branch
  template <class F>
  void _traverse_unipath(bit_vector &, size_t &, const F & f, ssize_t edge) {
    //auto visit = [&](size_t i) -> void { visited[i] = 1; num_visited++; };
    // Check if first kmer is in set (yes)> return

    //if (edge == -1 || (size_t)edge >= dbg.num_edges()) return;
    // visit elements in buffer to output first kmer

    symbol_type x;
    while(edge != -1 && b[edge] != has_branch) { // and edge not palindrome
      edge = g._forward(edge, x);
      f(x);
    }
    // Add last kmer to set
    // Call this recursively if branch
  }
};

template <class G>
unipath_visitor<G> make_unipath_visitor(const G & g) {
  return unipath_visitor<G>(g);
}
*/

#endif
