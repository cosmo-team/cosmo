#ifndef SORT_HPP
#define SORT_HPP

#import "kmer.hpp"

template <typename kmer_t>
struct get_key_colex_node
{
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
struct get_key_colex_edge
{
  typedef kmer_t key_type;
  key_type operator() (const kmer_t & obj) const { return obj; }
  key_type min_value() const { return std::numeric_limits<key_type>::min(); }
  key_type max_value() const { return std::numeric_limits<key_type>::max(); }
};

#endif
