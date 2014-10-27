#pragma once
#ifndef _DEBRUIJN_HYPERGRAPH_H
#define _DEBRUIJN_HYPERGRAPH_H

#include <boost/optional.hpp>
#include "debruijn_graph.hpp"

using namespace boost;

template <class t_debruijn_graph = debruijn_graph<>,
          class t_lcs_vector     = wt_huff<rrr_vector<63>> >
class debruijn_hypergraph : t_debruijn_graph {
  // Nodes are now ranges of edges (could be done for standard de bruijn graph too...
  // would save time by saving state?)
  // A node is a contiguous range of edges, a hypernode is a contiguous range of nodes (hence a range of edges)
  // A hypernode then replaces the defn of node without loss of generality
  public:
  typedef size_t                     edge_type;
  typedef pair<edge_type, edge_type> node_type; // static assert 1st <= 2nd
  //typedef optional<node_type> hypernode_type; // optional return type

  typedef typename t_debruijn_graph::symbol_type symbol_type;
  const t_debruijn_graph & m_dbg;
  const t_lcs_vector       m_lcs;

  debruijn_hypergraph(const t_debruijn_graph & dbg, const t_lcs_vector & lcs) : m_dbg(dbg), m_lcs(lcs) {}

  // shorter(v, k) - returns the hypernode whose label is the last k characters of v's label (reduce context)
  // longer(v, k) - list nodes (new "node") whose labels have length k <= K and end with v's label

  // maxlen(v, x) - returns some node in the *original* (kmax) graph whose label ends with v's
  // label, and that has an outgoing edge labelled x, or NULL otherwise
  node_type maxlen(const node_type & v) const {
    size_t start = get<0>(v); // have to update this to select to the next node
    // Range_start must also be a start of a node at the top level context, so find the next top level node to find the range
    // TODO: make (or check if it exists) uniform interface for rank/select, so I can easily make next(), prev(), etc
    size_t end = m_dbg._last_edge_of_node(m_dbg._edge_to_node(start));
    return node_type(start, end);
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
      return optional<node_type>(m_dbg._node_range(node_rank));
    }
    return optional<node_type>();
  }

  // function to get node by id
  // fwd, back, lastchar
};

#endif
