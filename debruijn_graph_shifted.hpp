#pragma once
#ifndef _DEBRUIJN_GRAPH_SHIFTED_H
#define _DEBRUIJN_GRAPH_SHIFTED_H

#include <algorithm>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <iostream>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "algorithm.hpp"
#include "utility.hpp"
#include "io.hpp"
//#include "debug.h"

using namespace std;
using namespace sdsl;

// TODO: convert asserts into exceptions? (Copy Boost)
template <size_t t_sigma            = 4, // default: DNA, TODO: change to 0 and make dynamic
          class  t_bit_vector_type  = sd_vector<>,
          class  t_bv_rank_type     = typename t_bit_vector_type::rank_0_type, // We invert the bits so it is 0 instead
          class  t_bv_select_type   = typename t_bit_vector_type::select_0_type,
          class  t_edge_vector_type = wt_huff<rrr_vector<63>>,
          class  t_symbol_type      = typename t_edge_vector_type::value_type,
          class  t_label_type       = string> // can define basic_string<t_symbol_type>, but need to use correct print func
class debruijn_graph_shifted {
  static_assert(t_sigma == 4, "Alphabet sizes other than 4 are not yet supported.");

  public:
  const static size_t sigma = t_sigma;
  typedef t_symbol_type symbol_type;
  typedef t_label_type  label_type;
  typedef typename t_bit_vector_type::size_type size_type;
  typedef size_t edge_type;
  typedef pair<edge_type, edge_type> node_type;

  const size_t           k{};

  const t_bit_vector_type      m_node_flags{};
  const t_bv_rank_type         m_node_rank{};
  const t_bv_select_type       m_node_select{};
  const t_edge_vector_type     m_edges{};
  // This is the "F table" in the blog/paper. It stores the starting positions of the sorted runs
  // of the k-1th symbol of the edge (or, last symbol of the node)
  // Could be implemented as the start positions, but here we store the cumulative sum (i.e. run ends)
  const array<size_t, 1+sigma> m_symbol_ends{};
  const array<size_t, 1+sigma> m_edge_max_ranks{};
  const label_type             m_alphabet{};
  const size_t                 m_num_nodes{};

  public:
  debruijn_graph_shifted() {}


  debruijn_graph_shifted(size_t in_k, const t_bit_vector_type & node_flags, const t_edge_vector_type & edges, const array<size_t, 1+sigma>& symbol_ends, const label_type& alphabet)
    : k(in_k), m_node_flags(node_flags), m_node_rank(&m_node_flags), m_node_select(&m_node_flags), m_edges(edges),
      m_symbol_ends(symbol_ends),
      m_edge_max_ranks(_init_max_ranks(edges)),
      m_alphabet(alphabet),
      m_num_nodes(m_node_rank(m_symbol_ends[sigma])) {
  }
  private:
  array<size_t, 1+sigma> _init_max_ranks(const t_edge_vector_type & edges) {
    array<size_t, 1+sigma> max_ranks;
    size_t num_edges = edges.size();
    for (symbol_type x = 0; x<sigma+1;x++) {
      max_ranks[x] = edges.rank(num_edges, _with_edge_flag(x, false));
    }
    return max_ranks;
  }

  public:
  static debruijn_graph_shifted load_from_packed_edges(istream & input, label_type alphabet=label_type{}, vector<size_t> * v=nullptr) {
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

    // TODO: sanity check the inputs (e.g. tally things, convert the above asserts)
    // So we avoid a huge malloc if someone gives us a bad file
    int_vector<1> first(num_edges,0);
    // would be nice to fix wavelet trees so the constructor
    // can accept a int_vector<4> instead (which is all we need for DNA)
    int_vector<8> edges(num_edges);

    bool prev_was_minus = false;
    for (size_t i = 0; i < num_edges; i++) {
      auto x = get_edge(blocks.begin(), i);
      first[i] = 1-get<1>(x); // convert 0s to 1s so we can have a sparse bit vector
      // For branchy graphs it might be better to change this and use RRR
      edges[i] = (get<0>(x) << 1) | !get<2>(x);
      if (v && get<0>(x) && !get<2>(x) && !prev_was_minus) {
        v->push_back(i);
        prev_was_minus = true;
      }
      else if (v && get<0>(x) && get<2>(x)) prev_was_minus = false;
    }

    t_bit_vector_type bv(first);
    t_edge_vector_type wt;
    construct_im(wt, edges);
    return debruijn_graph_shifted(k, bv, wt, counts, alphabet);
  }

  // Loaders/Writers for sdsl-serialized and JSON
  // size_in_bytes:

  // API
  size_t outdegree(size_t v) const {
    assert(v < num_nodes());
    auto range = _node_range(v);
    size_t first = get<0>(range);
    size_t last  = get<1>(range);
    size_t count = last - first + 1;
    // This edge access is kinda annoying, but if the ONE edge is $, we don't really have an outgoing edge
    return count - (count == 1 && _strip_edge_flag(m_edges[first]) == 0);
  }

  vector<node_type> all_preds(const node_type & v) const {
    // node u -> v : edge i -> j
    size_t j = get<0>(v);
    symbol_type y = _symbol_access(j);
    if (y == 0) return vector<node_type>(0);
    size_t i_first = _backward(j);
    size_t i_last  = _next_edge(i_first, y);
    size_t base_rank = m_edges.rank(i_first, _with_edge_flag(y, true));
    size_t last_rank = m_edges.rank(i_last, _with_edge_flag(y, true));
    size_t num_predecessors = last_rank - base_rank + 1;
    // binary search over first -> first + count;
    auto selector = [&](size_t i) -> size_t {
      return (i == 0)? i_first : m_edges.select(base_rank+i, _with_edge_flag(y,true));
    };

    vector<node_type> result(num_predecessors);
    for (size_t i = 0; i<num_predecessors; i++) {
      edge_type e_i = selector(i);
      edge_type e_j = _last_edge_of_node(_edge_to_node(e_i));
      result.push_back(node_type(e_i, e_j));
    }
    return result;
  }

  size_t indegree(size_t v) const {
    // find first predecessor edge i->j
    size_t j = _node_to_edge(v);
    // edge label has to be the last node symbol of v
    symbol_type x = _symbol_access(j);
    if (x == 0) return 0;
    size_t i_first = _backward(j);
    size_t i_last = _next_edge(i_first, x);
    return m_edges.rank(i_last, _with_edge_flag(x, true)) -
           m_edges.rank(i_first, _with_edge_flag(x, true)) + 1;
  }

  // The signed return type is worrying, since it halves the possible answers...
  // but it will only face problems with graphs that have over 2^63 ~= 4^31 edges,
  // which is a saturated debruijn graph of k=31. Which is unlikely.
  // Hopefully an iterator-based API (like BGL) will fix this issue
  ssize_t outgoing(size_t u, symbol_type x) const {
    ssize_t edge = outgoing_edge(u, x);
    if (edge == -1)
      return -1;
    else
      return _edge_to_node(edge);
  }

  // The signed return type is worrying, since it halves the possible answers...
  // but it will only face problems with graphs that have over 2^63 ~= 4^31 edges,
  // which is a saturated debruijn graph of k=31. Which is unlikely.
  // Hopefully an iterator-based API (like BGL) will fix this issue
  ssize_t outgoing_edge(size_t u, symbol_type x) const {
    assert(u < num_nodes());
    assert(x < sigma + 1);
    if (x == 0) return -1;
    auto range = _node_range(u);
    size_t first = get<0>(range);
    size_t last  = get<1>(range);
    // Try both with and without a flag
    for (symbol_type c = _with_edge_flag(x,false); c <= _with_edge_flag(x, true); c++) {
      size_t rnk = m_edges.rank(last+1, c);
      if (rnk == 0)
	continue;
      size_t most_recent = m_edges.select(rnk, c);
      // if within range, follow forward
      if (first <= most_recent && most_recent <= last) {
        // Don't have to check fwd for -1 since we checked for $ above
        return _forward(most_recent);
      }
    }
    return -1;
  }

  // Added for DCC. Will remove the other one later and rename this one.
  ssize_t interval_node_outgoing(const node_type & u, symbol_type x) const {
    //assert(u < num_nodes());
    assert(x < sigma + 1);
    if (x == 0) return -1;
    //auto range = _node_range(u);
    size_t first = get<0>(u);
    size_t last  = get<1>(u);
    // Try both with and without a flag
    for (symbol_type c = _with_edge_flag(x,false); c <= _with_edge_flag(x, true); c++) {
      size_t most_recent = m_edges.select(m_edges.rank(last+1, c), c);
      // if within range, follow forward
      if (first <= most_recent && most_recent <= last) {
        // Don't have to check fwd for -1 since we checked for $ above
        return _edge_to_node(_forward(most_recent));
      }
    }
    return -1;
  }



  // For DCC
  ssize_t _outgoing_edge_pair(size_t first, size_t last, symbol_type x) const {
    // Try both with and without a flag
    for (symbol_type c = _with_edge_flag(x,false); c <= _with_edge_flag(x, true); c++) {
      size_t most_recent = m_edges.select(m_edges.rank(last+1, c), c);
      // if within range, follow forward
      if (first <= most_recent && most_recent <= last) {
        // Don't have to check fwd for -1 since we checked for $ above
        return _forward(most_recent);
      }
    }
    return -1;
  }

inline  ssize_t incoming_edge(size_t v, symbol_type x) const {
    // This is very similar to indegree, so should maybe be refactored
    assert(v < num_nodes());
    assert(x < sigma + 1);
    // node u -> v : edge i -> j
    size_t j = _node_to_edge(v);
    symbol_type y = _symbol_access(j);
    if (y == 0) return -1;
    size_t i_first = _backward(j);
    size_t i_last  = _next_edge(i_first, y);
    size_t base_rank = m_edges.rank(i_first, _with_edge_flag(y, true));
    size_t last_rank = m_edges.rank(i_last, _with_edge_flag(y, true));
    size_t num_predecessors = last_rank - base_rank + 1;
    // binary search over first -> first + count;
    auto selector = [&](size_t i) -> size_t {
      return (i == 0)? i_first : m_edges.select(base_rank+i, _with_edge_flag(y,true));
    };
    auto accessor = [&](size_t i) -> symbol_type { return _first_symbol(selector(i)); };
    ssize_t sub_idx = function_binary_search(0, num_predecessors-1, x, accessor);
    if (sub_idx == -1) return -1;
    return selector(sub_idx);
  }

    
  // incoming
  ssize_t incoming(size_t v, symbol_type x) const {
    return _edge_to_node(incoming_edge(v, x));
  }

    
  // string -> node, edge
  // BGL style API

  label_type node_label(size_t v) const {
    size_t i = _node_to_edge(v);
    label_type label = label_type(k-1, _map_symbol(symbol_type{}));
    return _node_label_from_edge_given_buffer(i, label);
  }

  label_type node_label_from_edge(size_t i) const {
    label_type label = label_type(k-1, _map_symbol(symbol_type{}));
    return _node_label_from_edge_given_buffer(i, label);
  }

  label_type edge_label(size_t i) const {
    label_type label = label_type(k, _map_symbol(symbol_type{}));
    _node_label_from_edge_given_buffer(i, label);
    label[k-1] = _map_symbol(_strip_edge_flag(m_edges[i]));
    return label;
  }

  size_t num_edges() const { return m_symbol_ends[sigma]; /*_node_flags.size();*/ }
  size_t num_nodes() const { return m_num_nodes; /*m_node_rank(num_edges());*/ }

  private:
  size_t _symbol_start(symbol_type x) const {
    assert(x < sigma + 1);
    return (x==0)? 0 : m_symbol_ends[x - 1];
  }

  public:
  size_t _node_to_edge(size_t v) const {
    assert(v < num_nodes());
    return m_node_select(v+1);
  }

  size_t _edge_to_node(size_t i) const {
    assert(i < num_edges());
    return m_node_rank(i+1)-1;
  }

  // This should be moved to a helper file...
  symbol_type _strip_edge_flag(symbol_type x) const {
    return x >> 1;
  }

  // False -> normal edge (but still shifted to fit in the edge alphabet)
  // True -> minus flag
  symbol_type _with_edge_flag(symbol_type x, bool edge_flag) const {
    return (x << 1) | edge_flag;
  }

  symbol_type _symbol_access(size_t i) const {
    assert(i < num_edges());
    // I assume std binary search is optimised for small ranges (if sigma is small, e.g. DNA)
    return upper_bound(m_symbol_ends.begin(), m_symbol_ends.end(), i) - m_symbol_ends.begin();
  }

  node_type get_node(size_t v) const {
    auto r = _node_range(v);
    return node_type(get<0>(r), get<1>(r));
  }

  // provided for DCC paper
  inline symbol_type lastchar(const node_type & v) const {
    return _symbol_access(get<0>(v));
  }

  // more efficient than generating full label for incoming()
  symbol_type _first_symbol(size_t i) const {
    symbol_type x = 0;
    for (size_t pos = 1; pos <= k-1; pos++) {
      // Don't need to map it to the alphabet in this case
      x = _symbol_access(i);
      // All are $ before the last $
      if (x == 0) return x;
      i = _backward(i);
    }
    return x;
  }

  label_type _node_label_from_edge_given_buffer(size_t i, label_type & label) const {
    // Calculate backward k times and fill a buffer with symbol_access(edge)
    //label_type label = label_type(k-1, _map_symbol(symbol_type{}));
    for (size_t pos = 1; pos <= k-1; pos++) {
      symbol_type x = _symbol_access(i);
      label[k-pos-1] = _map_symbol(x);
      // All are $ before the last $
      if (x == 0) return label;
      i = _backward(i);
    }
    return label;
  }

  size_t _next_edge(size_t i, symbol_type x) const {
    if (i >= num_edges() - 1) return i;
    // Might not actually occur if out of rank bounds?
    size_t next_rank = 1 + m_edges.rank(1+i, _with_edge_flag(x, false));
    if (next_rank > m_edge_max_ranks[x]) return num_edges();
    return m_edges.select(next_rank, _with_edge_flag(x, false));
  }

  size_t _rank_distance(size_t a, size_t b) const {
    return m_node_rank(b) - m_node_rank(a);
  }

  // Return index of first possible edge obtained by following edge i
  public:
  ssize_t _forward(size_t i) const {
    symbol_type temp;
    return _forward(i, temp);
  }

  // This is so we can reuse the symbol lookup - save an access during traversal :)
  ssize_t _forward(size_t i, symbol_type & x) const {
    assert(i < num_edges());
    x = _strip_edge_flag(m_edges[i]);
    symbol_type fullx =_with_edge_flag(x, false);
    // if x == 0 ($) then we can't follow the edge
    // (should maybe make backward consistent with this, but using the edge 0 loop for node label generation).
    if (x == 0) return -1;

    // if this is flagged, then reset i to the corresponding unflagged symbol and use that as the starting point
    if (m_edges[i] & 1) {
      i = m_edges.select(m_edges.rank(i, fullx), fullx);
    }

    size_t start = _symbol_start(x);
    size_t nth   = m_edges.rank(i, fullx);
    size_t next  = m_node_select(m_node_rank(start+1) + nth);
    return next;
  }


  size_t _backward(size_t i) const {
    assert(i < num_edges());
    symbol_type x  = _symbol_access(i);
    // This handles x = $ so that we have the all-$ edge at position 0
    // but probably shouldn't be called this way in the normal case
    size_t x_start = _symbol_start(x);
    // rank is over [0,i) and select is 1-based
    size_t nth = _rank_distance(x_start+1, i+1);
    if (x == 0) return 0;
    // no minus flag because we want the FIRST
    return m_edges.select(nth+1, _with_edge_flag(x, false));
  }

  size_t backward(size_t v) const {
    return _edge_to_node(_backward(_node_to_edge(v)));
  }

  symbol_type _map_symbol(symbol_type x) const {
    return (m_alphabet.size() > 0)? m_alphabet[x] : x;
  }

  size_t _first_edge_of_node(size_t v) const {
    assert(v < num_nodes());
    // select is 1-based, but nodes are 0-based
    return m_node_select(v+1);
  }

  size_t _last_edge_of_node(size_t v) const {
    // find the *next* node's first edge and decrement
    // as long as a next node exists!
    assert(v + 1 <= num_nodes());
    if (v+1 == num_nodes()) return num_edges() - 1;
    else return _first_edge_of_node(v+1) - 1;
  }

  pair<size_t, size_t> _node_range(size_t v) const {
    return make_pair(_first_edge_of_node(v), _last_edge_of_node(v));
  }

  // TODO: add first_sibling and edge_range
  // TODO: update to use rank and select for larger alphabets
  size_t _last_sibling(size_t i) const {
    size_t last = i;
    while(last < num_edges() && m_node_flags[last] != 0) {
      last++;
    }
    // last should be one past end
    return last-1;
  }
    symbol_type _encode_symbol(uint8_t c) const {
        return lower_bound(m_alphabet.begin(), m_alphabet.end(), c) - m_alphabet.begin();
    }

    template <class InputIterator>
    boost::optional<node_type> index(InputIterator in) const {
        auto c = *in++;
        symbol_type first_symbol = _encode_symbol(c);
        // Range is from first edge of first, to last edge of last
        size_t start = _symbol_start(first_symbol);
        size_t end   = m_symbol_ends[first_symbol]-1;
        size_t first, last;

        // find c-labeled pred edge
        // if outside of range, find c- labeled pred edge
        for (size_t i = 0; i < k - 2; i++) {
            c = *in++;
            symbol_type x = _encode_symbol(c);
            // update range; Within current range, find first and last occurence of c or c-
            // first -> succ(x, first)
            for (uint8_t y=x<<1; y<(x<<1)+1; y++) {
                first = m_edges.select((m_edges.rank(start, y)) + 1, y);
                if (start <= first && first <= end) break;
            }
            if (!(start <= first && first <= end)) return boost::optional<node_type>();
            // last -> pred(x, last)
            if (start == end) {
                last = first;
            } else {
                for (uint8_t y=x<<1; y<(x<<1)+1; y++) {
                    last = m_edges.select((m_edges.rank(end + 1, y)), y);
                    if (start <= last && last <= end) break;
                }
            }
            assert(start <= last && last <= end);
            // Follow each edge forward
            start = _forward(first, x);
            end   = _forward(last, x);
            end   = _last_edge_of_node(_edge_to_node(end));
        }
        return boost::optional<node_type>(node_type(start, end));
    }

    
  size_type serialize(ostream& out, structure_tree_node* v=NULL, string name="") const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(k, out, child, "k");
    written_bytes += m_node_flags.serialize(out, child, "node_flags");
    written_bytes += m_node_rank.serialize(out, child, "node_rank");
    written_bytes += m_node_select.serialize(out, child, "node_select");
    written_bytes += m_edges.serialize(out, child, "edges");
    written_bytes += write_member(m_symbol_ends, out, child, "symbol_ends");
    written_bytes += write_member(m_edge_max_ranks, out, child, "edge_max_ranks");
    written_bytes += write_member(m_alphabet, out, child, "alphabet");
    written_bytes += write_member(m_num_nodes, out, child, "num_nodes");
    // helper bitvector m_doc_rmin_marked and m_doc_rmax_marked are not serialize
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  // THIS IS PROBABLY A BAD DESIGN CHOICE on behalf of sdsl
  // we shouldnt have a function which mutates the state? I dunno, I just have a bad feeling
  //! Loads the data structure from the given istream.
  void load(std::istream& in) {
    // DANGER WILL ROBINSON!
    // Yeah, don't normally const_cast!
    // but in this case the alternative is to leave them non-const,
    // thus sacrificing any gain...
    // plus I know what I'm doing here...
    read_member(deconst(k), in);
    deconst(m_node_flags).load(in);
    deconst(m_node_rank).load(in);
    deconst(m_node_rank).set_vector(&m_node_flags);
    deconst(m_node_select).load(in);
    deconst(m_node_select).set_vector(&m_node_flags);
    deconst(m_edges).load(in);
    read_member(deconst(m_symbol_ends), in);
    read_member(deconst(m_edge_max_ranks), in);
    read_member(deconst(m_alphabet), in);
    read_member(deconst(m_num_nodes), in);
  }

  size_type size() const { return num_edges(); }
};

#endif
