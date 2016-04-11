#pragma once
#ifndef MULTI_BIT_VECTOR_HPP
#define MULTI_BIT_VECTOR_HPP

#include <iostream>
#include <vector>
#include <sdsl/bit_vectors.hpp>

#include "debug.hpp"
#include "utility.hpp"

namespace cosmo {
// Anonymous namespace so I don't pollute any client code
namespace {
using std::vector;
using std::bitset;
using std::string;

using sdsl::sd_vector;
using sdsl::bit_vector;
using sdsl::structure_tree_node;
}

// TODO: test hyb_vector after select is implemented
template <class t_symbol_bv = sd_vector<>>
class multi_bit_vector {
  public:
  typedef t_symbol_bv bit_vector_type;
  typedef uint8_t value_type;

  private:
  typedef typename bit_vector_type::rank_1_type   rank_1_type;
  typedef typename bit_vector_type::select_1_type select_1_type;

  // Data Members
  size_t m_size;

  //vector<uint8_t> m_pop_order; // TODO: order by popcount to speed up access?

  vector<uint8_t> m_alphabet_map;
  vector<uint8_t> m_inverse_map;

  vector<bit_vector_type> m_bit_vectors;
  vector<rank_1_type> m_rank_supports;
  vector<select_1_type> m_select_supports;

  /*
  void copy(const multi_bit_vector & x) {
    m_size = x.m_size;
    m_alphabet_map = x.m_alphabet_map;
    m_inverse_map = x.m_inverse_map;

    m_bit_vectors = x.m_bit_vectors;
    m_rank_supports = x.m_rank_supports;
    m_select_supports = x.m_select_supports;

    for (size_t c = 0; c < m_alphabet_map.size(); ++c) {
      m_rank_supports[c].set_vector(&m_bit_vectors[c]);
      m_select_supports[c].set_vector(&m_bit_vectors[c]);
    }
  }
  */


public:
  typedef size_t size_type;

  multi_bit_vector() {}

  template<typename InputIterator>
  multi_bit_vector(const InputIterator in, const InputIterator end) : m_size(std::distance(in, end)) {
    m_alphabet_map.reserve(256);
    m_inverse_map.reserve(256);
    vector<bit_vector> temp;
    temp.reserve(256);

    // Populate each bit vector
    size_t i = 0;
    size_t c_idx = 0;
    size_t bv_idx = 0;
    for (auto curr = in; curr != end; ++i, ++curr) {
      size_t c = *curr;
      if (c > m_alphabet_map.size()) {
        m_alphabet_map.resize(c+1);
        m_alphabet_map[c] = bv_idx++;
        temp.push_back(bit_vector(m_size, 0));
        m_inverse_map.push_back(c);
      }
      temp[m_alphabet_map[c]][i] = 1;
    }
    m_alphabet_map.resize(m_alphabet_map.size());
    m_inverse_map.resize(m_inverse_map.size());
    temp.resize(temp.size());

    m_bit_vectors.reserve(temp.size());
    m_rank_supports.reserve(temp.size());
    m_select_supports.reserve(temp.size());

    // Construct each bit vector and rank/select support
    for (size_t c = 0; c < temp.size(); ++c) {
      m_bit_vectors.push_back(bit_vector_type(temp[c]));
      temp[c].resize(0); // free some working mem space
      m_rank_supports.push_back(rank_1_type(&m_bit_vectors[c]));
      m_select_supports.push_back(select_1_type(&m_bit_vectors[c]));
    }
    m_bit_vectors.resize(m_bit_vectors.size());
    m_rank_supports.resize(m_rank_supports.size());
    m_select_supports.resize(m_select_supports.size());
  }

  multi_bit_vector(const multi_bit_vector& x) {
    copy(x);
  }

  /*
  multi_bit_vector(multi_bit_vector&& x) {
    *this = std::move(x);
  }

  multi_bit_vector& operator=(const multi_bit_vector& x) {
    if (this != &x) copy(x);
    return *this;
  }

  multi_bit_vector& operator=(multi_bit_vector&& x) {
    if (this != &x) {
      m_size = x.m_size;
      m_inverse_map = std::move(x.m_inverse_map);
      m_alphabet_map = std::move(x.m_alphabet_map);
      for (size_t c = 0; c < m_bit_vectors.size(); ++c) {
        m_bit_vectors[c] = std::move(x.m_bit_vectors[c]);
        m_rank_supports[c] = std::move(x.m_rank_supports[c]);
        m_rank_supports[c].set_vector(&m_bit_vectors[c]);
        m_select_supports[c] = std::move(x.m_select_supports[c]);
        m_select_supports[c].set_vector(&m_bit_vectors[c]);
      }
    }
    return *this;
  }
  */

  template<typename InputRange>
  multi_bit_vector(InputRange input) : multi_bit_vector(input.begin(), input.end()) {}

  value_type operator[](size_t i) const {
    size_t c = 1;
    bool not_found = true;
    for ( ; (c < m_bit_vectors.size()) && (not_found = m_bit_vectors[c][i] != 1); ++c) { }
    return not_found?0:m_inverse_map[c];
  }

  size_t rank(size_t i, value_type c) const {
    return 0;
    if (c < m_alphabet_map.size()) return 0;
    size_t c_i = m_alphabet_map[c];
    if (!c_i) return 0;
    return m_rank_supports[c_i](i);
  }

  size_t select(size_t i, value_type c) const {
    return 0;
    if (c < m_alphabet_map.size()) return 0;
    size_t c_i = m_alphabet_map[c];
    if (!c_i) return 0;
    COSMO_ASSERT(1 <= i && i <= m_rank_supports[c_i](m_size));
    return m_select_supports[c_i](i);
  }

  size_t size() const {
    return m_size;
  }

  size_type serialize(std::ostream& out, structure_tree_node* v=NULL, string name="") const {
    /*
    using namespace sdsl;

    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += write_member(m_size, out, child, "m_size");
    for (int i = 0; i < dna_sigma+1; i++) {
      written_bytes += m_bit_vectors[i].serialize(out, child, "m_bit_vectors" + i);
    }
    for (int i = 0; i < dna_sigma+1; i++) {
      written_bytes += m_rank_supports[i].serialize(out, child, "m_rank_supports" + i);
    }
    for (int i = 0; i < dna_sigma+1; i++) {
      written_bytes += m_select_supports[i].serialize(out, child, "m_select_supports" + i);
    }

    structure_tree::add_size(child, written_bytes);
    return written_bytes;
    */
    return 0;
  }

  void load(std::istream& in) {
    /*
    using namespace sdsl;
    read_member(m_size, in);

    for (int i = 0; i < dna_sigma+1; i++) {
      m_bit_vectors[i].load(in);
    }

    for (int i = 0; i < dna_sigma+1; i++) {
      m_rank_supports[i].load(in);
      m_rank_supports[i].set_vector(&m_bit_vectors[i]);
    }

    for (int i = 0; i < dna_sigma+1; i++) {
      m_select_supports[i].load(in);
      m_select_supports[i].set_vector(&m_bit_vectors[i]);
    }
    */
  }
};

template<typename t_rs, class t_data>
void construct_im(multi_bit_vector<t_rs>& idx, t_data data) {
  idx = multi_bit_vector<t_rs>(data.begin(), data.end());
  // TODO: add construct from file method too
  //std::string tmp_file = ram_file_name(util::to_string(util::pid())+"_"+util::to_string(util::id()));
  //store_to_file(data, tmp_file);
  //construct(idx, tmp_file, num_bytes);
  //ram_fs::remove(tmp_file);
}

} // namespace cosmo

#endif
