#ifndef SORT_HPP
#define SORT_HPP

#include "kmer.hpp"

/*
template <class KeyAccessor, typename record_t>
struct record_wrap {
  KeyAccessor k;
  typedef typename KeyAccessor::key_type key_type;
  key_type operator() (const record_t & obj) const {
    return k(get<0>(obj));
  }
  // TODO: hacky if we want to extend it to different record types.
  record_t min_value() const { return record_t(k.min_value(), 0); }
  record_t max_value() const { return record_t(k.max_value(), 0xFFFFFFFFFFFFFFFF); }
};
*/

template <typename kmer_t>
struct get_key_colex_node {
  typedef kmer_t key_type;
  key_type operator() (const kmer_t & obj) const {
    kmer_t   temp = get_start_node(obj);
    uint64_t edge = get_edge_label(obj);
    ((uint64_t*)&temp)[0] |= edge;
    return temp;
  }
  key_type min_value() const { return std::numeric_limits<key_type>::min(); }
  key_type max_value() const { return std::numeric_limits<key_type>::max(); }
};

template <typename kmer_t>
struct get_key_colex_edge {
  typedef kmer_t key_type;
  key_type operator() (const kmer_t & obj) const { return obj; }
  key_type min_value() const { return std::numeric_limits<key_type>::min(); }
  key_type max_value() const { return std::numeric_limits<key_type>::max(); }
};

template <typename kmer_t>
struct kmer_less {
  typedef kmer_t value_type;
  bool operator() (const value_type & a, const value_type & b) const {
    return (a < b);
  }
  value_type min_value() const { return std::numeric_limits<value_type>::min(); }
  value_type max_value() const { return std::numeric_limits<value_type>::max(); }
};

template <typename kmer_t>
struct node_less {
  typedef kmer_t value_type;
  bool operator() (const value_type & a, const value_type & b) const {
    return get_start_node(a) < get_start_node(b);
  }
  value_type min_value() const { return std::numeric_limits<value_type>::min(); }
  value_type max_value() const { return std::numeric_limits<value_type>::max(); }
};

template <typename dummy_t>
struct colex_dummy_less {
  typedef dummy_t value_type;
  typedef typename dummy_t::first_type kmer_t;
  get_key_colex_node<kmer_t> colex_node;

  bool operator() (const value_type & a, const value_type & b) const {
    kmer_t   a_node = get_start_node(a.first);
    kmer_t   b_node = get_start_node(b.first);
    uint64_t a_edge = get_edge_label(a.first);
    uint64_t b_edge = get_edge_label(b.first);

    // If string equal then sort by length, else sort by string
    if (a_node == b_node) {
      if (a.second == b.second) {
        return a_edge < b_edge;
      }
      else return (a.second < b.second);
    }
    else return (a_node < b_node);
  }
  value_type min_value() const { return std::make_pair(colex_node.min_value(), std::numeric_limits<size_t>::min()); }
  value_type max_value() const { return std::make_pair(colex_node.max_value(), std::numeric_limits<size_t>::max()); }
};

#endif
