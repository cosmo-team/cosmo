#pragma once
#ifndef _ALGORITHM_HPP
#define _ALGORITHM_HPP

#include <iostream>

#include <sdsl/bit_vectors.hpp>
#include "debruijn_graph.hpp"
#include "kmer.hpp"
#include "uint128_t.hpp"

using namespace std;
//using namespace sdsl;

// Calling it with a vector of first minus positions is sliiiightly faster
// But would mean we have to save another vector with our graph if we ever want to traverse
// it again.
// If a vector isnt provided, another method is used which is *almost as fast*
template <bool has_branch = 1, class G, class V>
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


// Should change this to an edge iterator sometime
// Or refactor it to unipath_debruijn_graph
template <const bool has_branch = 1, class G, class B==sdsl::sd_vector<>>
class unipath_visit_strategy {
  typedef sdsl::bit_vector bit_vector;
  typedef typename _type_getter<has_branch, B>::rank_type rank_type;
  typedef typename _type_getter<has_branch, B>::select_type select_type;
  typedef typename G::symbol_type symbol_type;
  typedef typename G::label_type  label_type;

  public: // change this later
  const G & dbg;
  const B   branches;
  const branch_rank;
  const branch_select;

  // CTORs
  unipath_visit_strategy(const G& g) : dbg(g), branches(make_branch_vector(g)),
                                       branch_rank(&branches), branch_select(&branches) {}
  // Slightly faster version if user provides graph constructor with a size_t vector reference
  unipath_visit_strategy(const G& g, const vector<size_t> & v) : dbg(g), branches(make_branch_vector(g, v)),
                                                                 branch_rank(&branches), branch_select(&branches) {}

  template <class F>
  void operator()(F & f) {
    bit_vector visited(g.num_edges());
  }

  // _traverse_roots
  // traverse the dummy nodes depth first until we get to a real contig
  // then follow the contig until a $ out edge
  // finish each contig with $
  // when all the branches are done, pop
}

template <const bool has_branch = 1, class G, class B, class F, typename kmer_t>
void _kmer_t_visit_unipaths(const G& g, const B& b, const F &) {
  rank_type   branch_rank(&b);
  select_type branch_select(&b);
  size_t      num_branches = branch_rank(b.size());

  // Need to keep track of visited incase we have isolated cycles
  // This might sound unusual, but I think it might be frequent with plasmids
  // in bacteria
  sdsl::bit_vector  visited(g.num_edges(),0);
  size_t            num_visited = 0;
  auto              visit = [&](size_t i) -> void { visited[i] = 1; num_visited++; };

  // Since we are keeping track of visited edges, for all non-isolated cycles
  // we can traverse from the root (since all edges have incoming edges, and the
  // only ones without an entry in branching vector are isolated cycles, meaning
  // every other edge is accessible by traversing from the root)
  // The root is the node that ends with $
  // Set up root range [root_start, root_end):
  size_t root_start = 0;
  // If this is 0, we only have isolated cycles or nothing at all
  size_t root_end   = m.m_symbol_ends[0];

  /*
  for (size_t branch_idx = 1; branch_idx <= num_branches; branch_idx++) {
    // Find first branch edge
    size_t root = branch_select(branch_idx);
    // Get previous symbols
    // NOTE: If we used a stack we could avoid this call and start from the $ signs only...
    //       but it might take a lot of memory for very branchy graphs, so leave it for now
    //       ACTUALLY: instead, save previous node label only and do the edge select way of traversing
    //       Then, when we select to the next edge, we check if it has been visited already
    //       If it hasnt, the previous string will be the prefix of it (maybe save the edge too)
    // TODO: change this to kmer_t for easy palindrome check/set?
    auto root_node_label = g.node_label_from_edge(root);
    // TODO: if $ in label, skip it (but check to k+1, because if it is

    // iterate over sibling linearly, incrementing sibling + branch (cuts out some selects)
    // TODO: check that this wont make the above ++ off by one

    size_t num_siblings = g.num_siblings(root);
    for (size_t sibling = 0; sibling < num_siblings; sibling++, branch_idx++ ) {
      // TODO: print only last k-1 characters unless we are at k depth from dummy node start
      std::cerr << root_node_label << std::endl;

      // follow forward until we can't anymore (we meet a branch or $)
      symbol_type x;
      for (ssize_t edge = g._forward(sibling, x);
          // Check for palindromicness, make fn
          x != 0 && b[edge] != has_branch && edge != -1;
          edge = g._forward(edge,x)) {
        // Check if in set
        //visit(edge);
        x = g._strip_edge_flag(x);
        std::cerr << "acgt"[x-1];
      }
      std::cerr << endl;
    }
  }
  */

  // Anything unvisited has no incoming or outgoing branches -> floating cycles (potential plasmids)
  // (not a linear de bruijn graph, because we require every edge have predecessors -> dummies to the root)
  // handle_branchless_contigs(g, left_to_visit, f);
  // NOTE: this means we can handle every other case with handle_branchy_contigs and start from $$$$$Xs
  // (since we need to keep track of visited anyway)

  // TODO: remember representative{lastnode}, check representative{firstnode} in set, then don't traverse
  // Can do this with postprocessing though

  /*
  for (size_t sel_idx = 1; sel_idx <= rank(v.size()); sel_idx++) {
    kmer_t buf{}; // reset buffer
    size_t edge = select(sel_idx);
    assert(v[edge] == has_branch);
    std::cerr << "Edge Label: " << g.edge_label(edge) << std::endl;
    // Traverse back k edges to build buffer
    for (size_t e = edge, i = 0; i < g.k; i++) {
      buf = kmer_t{}; // clear buffer
      typename G::symbol_type x = g._symbol_access(e);
      if (x == 0) continue; // if $ skip this
      // x - 1 because we have to correct for $ being 0
      //buf = follow_edge(buf, x-1, g.k);
      e = g._backward(e);
    }
    std::cerr << kmer_to_string(buf, g.k) << std::endl;

    typename G::symbol_type x = g._strip_edge_flag(g.m_edges[edge]);
    if (x == 0) continue;
    buf = follow_edge(buf, x-1, g.k);
    if (is_palindrome(buf, g.k)) {
      continue;
    }

    for (size_t i = 0; i < g.k; i++) {
      uint8_t y = get_nt(buf, g.k - 1 - i);
      f(g._map_symbol(y));
    }

    edge = g._forward(edge);
    // For each edge... each one should be a has_branch
    // Due to the definition of our branch vector
    // so it will be selected to... so we just follow _forward
    while(v[edge] != has_branch) {
      x = g._strip_edge_flag(g.m_edges[edge]);
      if (x != 0) {
        // update buffer
        buf = follow_edge(buf, x-1, g.k);
        if (is_palindrome(buf, g.k)) {
          f(g._map_symbol(0));
          break;
        }
      }

      f(g._map_symbol(x));
      if (x == 0) break;
      edge = g._forward(edge);
      // for pairs: check here and at the if x == 0, and return node index of this + first edge
    }
  }*/
}

template <const bool has_branch = 1, class G, class V, class F>
void visit_unipaths(const G& g, const V& v, const F & f) {
  size_t k = g.k;
  if (k <= 32) return _kmer_t_visit_unipaths<has_branch, G, V, F, uint64_t>(g, v, f);
  if (k <= 64) return _kmer_t_visit_unipaths<has_branch, G, V, F, uint128_t>(g, v, f);
  assert(false); // no support for larger k yet
}

#endif
