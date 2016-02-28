#pragma once
#ifndef DNA_BV_RS_H
#define DNA_BV_RS_H

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

using sdsl::sd_vector;
using sdsl::hyb_vector;
using sdsl::bit_vector;
}

// TODO: make t_symbol_bv, t_minus_bv, t_terminal_bv all parameters
// TODO: test hyb_vector after select is implemented
template <class t_symbol_bv = sd_vector<>>
class dna_bv_rs {
private:
  const static size_t dna_sigma = 8; // DNA + minus flags. $ considered separately

  typedef t_symbol_bv symbol_bv;
  typedef uint8_t value_type;

  // Data Members
  size_t m_size;
  array<symbol_bv, dna_sigma+1> m_bit_vectors;
  typedef typename symbol_bv::rank_1_type   rank_1_type;
  typedef typename symbol_bv::select_1_type select_1_type;
  array<rank_1_type, dna_sigma+1> m_rank_supports;
  array<select_1_type, dna_sigma+1> m_select_supports;

public:
  template<typename InputIterator>
  dna_bv_rs(InputIterator in, InputIterator end) : m_size(std::distance(in, end)) {
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

  dna_bv_rs(const std::string & input) : dna_bv_rs(input | transformed(nt_to_int())) {}

  dna_bv_rs(const value_type * input) : dna_bv_rs(std::string(input)) {}

  value_type operator[](size_t i) const {
    size_t c = 0;
    for ( ; (c < dna_sigma + 1) && (m_bit_vectors[c][i] != 1); ++c) { }
    return int_to_nt()(c);
  }

  size_t rank(size_t i, value_type c) const {
    ssize_t c_i = nt_to_int()(c);
    if (c_i < 0) return 0;
    return m_rank_supports[nt_to_int()(c)](i);
  }

  size_t select(size_t i, value_type c) const {
    ssize_t c_i = nt_to_int()(c);
    if (c_i < 0) return 0;
    COSMO_ASSERT(1 <= i && i <= m_rank_supports[c_i](m_size));
    return m_select_supports[c_i](i);
  }

  size_t size() const {
    return m_size;
  }
};

} // namespace index
} // namespace cosmo

#endif
