#pragma once
#ifndef _DEBRUIJN_GRAPH_H
#define _DEBRUIJN_GRAPH_H

#include <algorithm>
#include <fstream>
#include <iterator>
#include <vector>
#include <array>
#include <string>
#include <iostream>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "algorithm.hpp"
#include "utility.hpp"
#include "io.hpp"
#include "debug.hpp"

using namespace std;
using namespace sdsl;

// TODO: convert asserts into exceptions?
// TODO: named template paremeters (See Boost)
template <class  t_edge_vector_type  = wt_huff<rrr_vector<63>>,
          class  t_bit_vector_type   = sd_vector<>,
          class  t_dummy_vector_type = vector<kmer_t>> // used to store dummies
class debruijn_graph {
  public:
  typedef string label_type;
  const static size_t sigma = 4;
  typedef t_edge_vector_type edge_vector_type;
  typedef t_bit_vector_type  bit_vector_type;
  typedef t_dummy_vector_type  dummy_vector_type;
  typedef typename edge_vector_type::value_type symbol_type;
  typedef typename bit_vector_type::size_type size_type;
  typedef typename dummy_vector_type::value_type kmer_t;
  typedef size_t edge_type;
  typedef pair<edge_type, edge_type> node_type;

  const size_t           k{};

  const bit_vector_type                         m_node_flags{};
  const typename bit_vector_type::rank_0_type   m_node_rank{};
  const typename bit_vector_type::select_0_type m_node_select{};
  const edge_vector_type     m_edges{};
  // This is the "F table" in the blog/paper. It stores the starting positions of the sorted runs
  // of the k-1th symbol of the edge (or, last symbol of the node)
  // Could be implemented as the start positions, but here we store the cumulative sum (i.e. run ends)
  const array<size_t, 1+sigma> m_symbol_ends{};
  const array<size_t, 1+sigma> m_edge_max_ranks{};
  const label_type             m_alphabet{};
  const size_t                 m_num_nodes{};
  const t_bit_vector_type m_dummy_flags;
  const typename t_bit_vector_type::rank_1_type m_dummy_rank_1;
  // TODO: Add IO for these
  const typename t_bit_vector_type::rank_0_type m_dummy_rank_0;
  const typename t_bit_vector_type::select_0_type m_dummy_select_0;
  // TODO: support STXXL vectors as dummy labels may not need to be in memory
  const vector<kmer_t> m_dummies;

  public:
  debruijn_graph() {}

  // TODO: make each of these a range instead (so it doesnt depend on the type)
  debruijn_graph(size_t in_k, const t_bit_vector_type & node_flags, const t_edge_vector_type & edges, const array<size_t, 1+sigma>& symbol_ends, const label_type& alphabet,
      const t_bit_vector_type & dummy_flags, const dummy_vector_type & dummies)
    : k(in_k), m_node_flags(node_flags), m_node_rank(&m_node_flags), m_node_select(&m_node_flags),
      m_edges(edges),
      m_symbol_ends(symbol_ends),
      m_edge_max_ranks(_init_max_ranks(m_edges)),
      m_alphabet(alphabet),
      m_num_nodes(m_node_rank(m_node_flags.size())),
      m_dummy_flags(dummy_flags),
      m_dummy_rank_1(&m_dummy_flags),
      m_dummy_rank_0(&m_dummy_flags),
      m_dummy_select_0(&m_dummy_flags),
      m_dummies(dummies) {}

  private:
  array<size_t, 1+sigma> _init_max_ranks(const t_edge_vector_type & edges) {
    array<size_t, 1+sigma> max_ranks;
    size_t num_edges = edges.size();
    for (symbol_type x = 0; x<sigma+1;x++) {
      max_ranks[x] = edges.rank(num_edges, _with_edge_flag(x, false));
    }
    return max_ranks;
  }

  bool is_incoming_dummy(size_t edge) const {
    // TODO: get start edge of node
    // should be called at first edge of node
    return m_dummy_flags[edge];
  }

  public:
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
    // assert(v < num_nodes());
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
    if (is_incoming_dummy(j)) return 0;

    // edge label has to be the last node symbol of v
    symbol_type x = _symbol_access(j);
    if (x == 0) return 0;
    size_t i_first = _backward(j, false);
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
      if (rnk == 0) continue;
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

  // incoming
  ssize_t incoming(size_t v, symbol_type x) const {
    // This is very similar to indegree, so should maybe be refactored
    assert(v < num_nodes());
    assert(x < sigma + 1);
    // node u -> v : edge i -> j
    size_t j = _node_to_edge(v);
    //COSMO_LOG(info) << "j: " << j;
    if (is_incoming_dummy(j)) return -1;
    symbol_type y = _symbol_access(j);
    if (y == 0) return -1;
    size_t i_first = _backward(j, false);
    size_t i_last  = _next_edge(i_first, y);
    //COSMO_LOG(info) << "i_first: " << i_first;
    //COSMO_LOG(info) << "i_last: " << i_last;
    size_t base_rank = m_edges.rank(i_first, _with_edge_flag(y, true));
    size_t last_rank = m_edges.rank(i_last, _with_edge_flag(y, true));
    size_t num_predecessors = last_rank - base_rank + 1;
    // binary search over first -> first + count;
    auto selector = [&](size_t i) -> size_t {
      return (i == 0)? i_first : m_edges.select(base_rank+i, _with_edge_flag(y,true));
    };
    auto accessor = [&](size_t i) -> symbol_type { return _first_symbol(selector(i)); };
    // TODO: rewrite this with iterators and boost function adapter? many binary searches have issues
    ssize_t sub_idx = function_binary_search(0, num_predecessors-1, x, accessor);
    if (sub_idx == -1) return -1;
    return _edge_to_node(selector(sub_idx));
  }

  // string -> node, edge
  // BGL style API
  size_t _get_dummy_index(size_t i) const {
    size_t dummy_idx = m_dummy_rank_1(i) - (i - m_node_rank(i));
    return dummy_idx;
  }

  label_type node_label(size_t v) const {
    size_t i = _node_to_edge(v);
    if (is_incoming_dummy(i)) {
      return kmer_to_string(m_dummies[_get_dummy_index(i)], k-1);
    }
    label_type label = label_type(k-1, _map_symbol(symbol_type{}));
    return _node_label_from_edge_given_buffer(i, label);
  }

  label_type node_label_from_edge(size_t i) const {
    if (is_incoming_dummy(i)) {
      return kmer_to_string(m_dummies[_get_dummy_index(i)], k-1);
    }
    label_type label = label_type(k-1, _map_symbol(symbol_type{}));
    return _node_label_from_edge_given_buffer(i, label);
  }

  label_type edge_label(size_t i) const {
    if (is_incoming_dummy(i)) {
      return kmer_to_string(m_dummies[_get_dummy_index(i)], k-1);
    }
    label_type label = label_type(k, _map_symbol(symbol_type{}));
    _node_label_from_edge_given_buffer(i, label);
    label[k-1] = _map_symbol(_strip_edge_flag(m_edges[i]));
    return label;
  }

  // TODO: subtract number of $ edges?
  size_t num_edges() const { return m_node_flags.size(); /* - m_edge_max_ranks[0];*/ }
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
    assert(i < size());
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
    assert(i < size());
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
  
  // Just return the symbol relating to this edge (for walks..)
  typename label_type::value_type edge_symbol(size_t i) const
  {
    return _map_symbol(_strip_edge_flag(m_edges[i]));
  }

  // more efficient than generating full label for incoming()
  symbol_type _first_symbol(size_t i) const {
    symbol_type x = 0;
    for (size_t pos = 1; pos <= k-1; pos++) {
      x = _symbol_access(i);
      if (x == 0) return x;
      if (is_incoming_dummy(i)) {
        // Get character in question
        symbol_type x = ((m_dummies[_get_dummy_index(i)]<<((k - pos - 1)*2)&0x11)+1);
        return x;
      }
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

      if (is_incoming_dummy(i)) {
        // If last character, don't need to prepend anything
        if (pos + 1 == k) return label;
        string label_2 = kmer_to_string(m_dummies[_get_dummy_index(i)], k - 1);
        string prefix = kmer_to_string(m_dummies[_get_dummy_index(i)]<<((pos-1)*2), k - 1 - pos);
        for (size_t i = 0; i < k-pos-1; ++i) {
          label[i] = prefix[i];
        }
        return label;
      }

      // TODO: optimize these node/edge conversion calls (does backward convert to nodes in the process?)
      i = _node_to_edge(_edge_to_node(_backward(i,false)));
    }
    return label;
  }

  size_t _next_edge(size_t i, symbol_type x) const {
    if (i >= size() - 1) return i;
    // Might not actually occur if out of rank bounds?
    size_t next_rank = 1 + m_edges.rank(1+i, _with_edge_flag(x, false));
    if (next_rank > m_edge_max_ranks[x]) return size();
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
  // Return index of the FIRST edge of the node pointed to by edge i.
  ssize_t _forward(size_t i, symbol_type & x) const {
    assert(i < size());
    auto sym_with_flag = m_edges[i];
    x = _strip_edge_flag(sym_with_flag);
    size_t start = _symbol_start(x);
    symbol_type fullx = _with_edge_flag(x, false);
    // if x == 0 ($) then we can't follow the edge
    // (should maybe make backward consistent with this, but using the edge 0 loop for node label generation).
    if (x == 0) return -1;

    size_t nth = m_edges.rank(i+1, fullx);
    // find prev_rank of dummies from base, then select to the nth
    //nth += (m_dummy_rank(m_symbol_ends[x]) - m_dummy_rank(start));

    // if this is flagged, then reset i to the corresponding unflagged symbol and use that as the starting point
    //if (sym_with_flag & 1) {
      // i = m_edges.select(m_edges.rank(i, fullx), fullx);
      // The above is basically the same as subtracting 1...
      // ALTERNATIVE PATCH which doesn't need extra ranks/selects (observed and written by Mike Mueller)
      /* Since rank is on 0..i-1, nth always reflects the rank of the symbol below what i points to.
         We implicitly add one to the result though because we always start with the rank of the first node after start so
         nth only needs to be 0 in order to select that node.
         
         If, however, i points to an edge flagged with -, then this is an extra incoming edge that enters the
         same node as another similarly flagged edge, so we don't want to have that implicit +1  
         
         This assumes that the flagged edge matches an edge that occurs BEFORE it in the list, which is what I 
         observe in the paper */
     // nth--;
    //}

    // select to nth 0 from base in dummies
    size_t prev = m_dummy_rank_0(start);
    size_t next = m_dummy_select_0(prev + nth);
    //size_t next = m_node_select(m_node_rank(start+1) + nth);

    /*
    COSMO_LOG(info) << "";
    COSMO_LOG(info) << "edge  : " << i;
    COSMO_LOG(info) << "start : " << start;
    COSMO_LOG(info) << "nth   : " << nth;
    COSMO_LOG(info) << "#prev : " << prev;
    COSMO_LOG(info) << "next  : " << next;
    */

    return next;
  }

  ssize_t _backward(size_t i, bool find_first=true) const {
    //COSMO_LOG(info) << "";
    //COSMO_LOG(info) << "i       : " << i;

    assert(i < size());
    // Makes things a lot easier if we just start from the first edge
    i = (find_first)?_node_to_edge(_edge_to_node(i)):i;
    //COSMO_LOG(info) << "i'      : " << i;
    if (is_incoming_dummy(i)) return -1;
    symbol_type x  = _symbol_access(i);
    // This handles x = $ so that we have the all-$ edge at position 0
    // As we use the same symbol for outgoing dummy edges, which actually
    // DONT point back to the incoming dummy edges.
    // NOTE: this will only happen if we added all incoming dummy edge shifts
    if (x == 0) return 0;
    size_t x_start = _symbol_start(x);
    // rank is over [0,i) and select is 1-based
    size_t nth = m_dummy_rank_0(i+1) - m_dummy_rank_0(x_start);//_rank_distance(x_start, i+1);
    //size_t nth_adjust = 0;//(m_dummy_rank_1(i+1) - m_dummy_rank_1(x_start));
    size_t next = m_edges.select(nth, _with_edge_flag(x, false));

    /*
    COSMO_LOG(info) << "x_start : " << x_start;
    COSMO_LOG(info) << "nth     : " << nth;
    COSMO_LOG(info) << "nth_d   : " << nth_adjust;
    COSMO_LOG(info) << "next    : " << next;
    */

    return next;
  }

  size_t backward(size_t v) const {
    return _edge_to_node(_backward(_node_to_edge(v)));
  }

  symbol_type _map_symbol(symbol_type x) const {
    return (m_alphabet.size() > 0)? m_alphabet[x] : x;
  }

  size_t _first_edge_of_node(size_t v) const {
    assert(v < num_nodes());
    // Why +1? because select is 1-based, but nodes are 0-based
    return m_node_select(v+1);
  }

  size_t _last_edge_of_node(size_t v) const {
    // find the *next* node's first edge and decrement
    // as long as a next node exists!
    // TODO - test this
    assert(v + 1 <= num_nodes());
    if (v+1 == num_nodes()) return size() - 1;
    else return _first_edge_of_node(v+1) - 1;
  }

  pair<size_t, size_t> _node_range(size_t v) const {
    return make_pair(_first_edge_of_node(v), _last_edge_of_node(v));
  }

  // TODO: add first_sibling and edge_range
  // TODO: update to use rank and select for larger alphabets
  size_t _last_sibling(size_t i) const {
    size_t last = i;
    while(last < size() && m_node_flags[last] != 0) {
      last++;
    }
    // last should be one past end
    return last-1;
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
    written_bytes += m_dummy_flags.serialize(out, child, "dummy_flags");
    written_bytes += m_dummy_rank_1.serialize(out, child, "dummy_rank_1");
    written_bytes += m_dummy_rank_0.serialize(out, child, "dummy_rank_0");
    written_bytes += m_dummy_select_0.serialize(out, child, "dummy_select_0");
    written_bytes += sdsl::serialize(m_dummies, out, child, "dummies");
    // helper bitvector m_doc_rmin_marked and m_doc_rmax_marked are not serialize
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Loads the data structure from the given istream.
  void load(std::istream& in) {
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
    deconst(m_dummy_flags).load(in);
    deconst(m_dummy_rank_1).load(in);
    deconst(m_dummy_rank_1).set_vector(&m_dummy_flags);
    deconst(m_dummy_rank_0).load(in);
    deconst(m_dummy_rank_0).set_vector(&m_dummy_flags);
    deconst(m_dummy_select_0).load(in);
    deconst(m_dummy_select_0).set_vector(&m_dummy_flags);
    sdsl::load(deconst(m_dummies), in);
  }

  size_type size() const { return m_node_flags.size(); }

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
        size_t first = 0, last = 0;

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
                    auto rank_temp = m_edges.rank(end + 1, y);
                    last = m_edges.select((rank_temp), y);
                    if (start <= last && last <= end) break;
                }
            }
            if (!(start <= last && last <= end)) {
                assert(!"(start <= last && last <= end)");
            }
        
            // Follow each edge forward
            start = _forward(first, x);
            end   = _forward(last, x);
            end   = _last_edge_of_node(_edge_to_node(end));
        }
        return boost::optional<node_type>(node_type(start, end));
    }
    
};

#endif
