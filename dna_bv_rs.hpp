#pragma once
#ifndef DNA_BV_RS_H
#define DNA_BV_RS_H

#include <iostream>
#include <sdsl/bit_vectors.hpp>

#include "utility.hpp"

namespace cosmo {
namespace index {
// Anonymous namespace so I don't pollute any client code
namespace {
using std::cout;
using std::endl;
using std::vector;
using std::array;
using std::bitset;
using std::string;

using sdsl::sd_vector;
using sdsl::rrr_vector;
using sdsl::hyb_vector;
using sdsl::bit_vector;
using sdsl::structure_tree_node;
}

// TODO: make t_symbol_bv, t_minus_bv, t_terminal_bv all parameters
// TODO: test hyb_vector after select is implemented
template <class t_symbol_bv = sd_vector<>>
class dna_bv_rs {
private:
  const static size_t dna_sigma = 8; // DNA + minus flags. $ considered separately

  typedef t_symbol_bv symbol_bv;
  typedef uint64_t value_type;

  // Data Members
  size_t m_size;
  array<symbol_bv, dna_sigma+1> m_bit_vectors;
  typedef typename symbol_bv::rank_1_type   rank_1_type;
  typedef typename symbol_bv::select_1_type select_1_type;
  array<rank_1_type, dna_sigma+1> m_rank_supports;
  array<select_1_type, dna_sigma+1> m_select_supports;

public:
  typedef size_t size_type;

  dna_bv_rs() {}

  template<typename InputIterator>
  dna_bv_rs(const InputIterator in, const InputIterator end) : m_size(std::distance(in, end)) {
    array<bit_vector, dna_sigma+1> temp;
    for (size_t c = 0; c < dna_sigma+1; ++c) {
      temp[c] = bit_vector(m_size, 0);
    }

    // Populate each bit vector
    size_t i = 0;
    for (auto curr = in; curr != end; ++i, ++curr) {
      size_t c = *curr;
      temp[c][i] = 1;
    }

    // Construct each bit vector
    for (size_t c = 0; c < dna_sigma + 1; ++c) {
      m_bit_vectors[c] = symbol_bv(temp[c]);
      temp[c].resize(0); // free some working mem space
      m_rank_supports[c] = rank_1_type(&m_bit_vectors[c]);
      m_select_supports[c] = select_1_type(&m_bit_vectors[c]);
    }
  }

  template<typename InputRange>
  dna_bv_rs(InputRange input) : dna_bv_rs(input.begin(), input.end()) {}

  // TODO: set a flag that maps ascii back...
  //dna_bv_rs(const std::string & input) : dna_bv_rs(input | transformed(nt_to_int())) {}
  //dna_bv_rs(const char * input) : dna_bv_rs(std::string(input)) {}

  value_type operator[](size_t i) const {
    size_t c = 1; // ignore the uncommon symbol $ until end
    bool not_found = false;
    // Shave the last iteration off... makes a tiny difference
    for ( ; (c < dna_sigma+1) && (not_found = m_bit_vectors[c][i] != 1); ++c) { }
    return not_found?0:c;//int_to_nt()(not_found?0:c);
  }

  size_t rank(size_t i, value_type c) const {
    ssize_t c_i = c;//nt_to_int()(c);
    if (c_i < 0) return 0;
    return m_rank_supports[c_i](i);
  }

  size_t select(size_t i, value_type c) const {
    ssize_t c_i = c;//nt_to_int()(c);
    if (c_i < 0) return 0;
    COSMO_ASSERT(1 <= i && i <= m_rank_supports[c_i](m_size));
    return m_select_supports[c_i](i);
  }

  size_t size() const {
    return m_size;
  }

  size_type serialize(std::ostream& out, structure_tree_node* v=NULL, string name="") const {
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
  }

  void load(std::istream& in) {
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
  }


};

} // namespace index
} // namespace cosmo

#endif
