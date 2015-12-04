#pragma once
#ifndef DNA_RS_H
#define DNA_RS_H

#include <bitset>
#include <vector>
#include <array>
#include <iterator>
#include <iostream>
#include <functional>
#include <numeric>
#include <algorithm>
#include <parallel/numeric>
#include <parallel/algorithm>
#include <boost/range/adaptors.hpp>
#include <boost/iterator/transform_iterator.hpp>
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
using sdsl::bit_vector;
using namespace boost::adaptors;
using namespace boost;
}

// TODO: access win/mac/linux programs that get this for us
// If this were in a GPU this should be 1024
const size_t cache_line_bits = 512;

template <size_t w, size_t i>
struct popcount_mask_r {
  const static uint64_t value = (1L << (i*w)) | popcount_mask_r<w, i-1>::value;
};

template <size_t w>
struct popcount_mask_r<w,0> {
  const static uint64_t value = 1;
};

template <size_t w>
struct popcount_mask {
  const static uint64_t value = popcount_mask_r<w,
    bitwidth<uint64_t>::width/w - 1>::value;
};

// Wider symbol popcount
template <size_t w>
struct popcounter {
  const size_t m_c;
  popcounter(size_t c=1) : m_c(c) {};

  template <typename block_t>
  inline size_t operator()(const block_t & x, size_t c) const {return popcount(x,c);}
  template <typename block_t>
  inline size_t operator()(const block_t & x) const {return popcount(x,m_c);}

  template <typename block_t>
  static inline block_t bitmap(const block_t & x, size_t c) {
    static_assert(w <= bitwidth<block_t>::width, "block width too small for symbol");
    assert(c < (1<<w));
    static_assert(popcount_mask<3>::value == 0x1249249249249249L, "");
    block_t mask = popcount_mask<w>::value;
    block_t temp = x ^ (c * mask);
    block_t result = temp;
    for (size_t i = 0; i < w; ++i) {
      result |= (temp >> i);
    }
    result &= mask;
    return result;
  }

  template <typename block_t>
  static inline size_t popcount(const block_t & x, size_t c) {
    block_t result = bitmap(x, c);
    const size_t block_bits  = bitwidth<block_t>::width; 
    const size_t block_width = block_bits/w; 
    return block_width - bitset<block_bits>(result).count();
  }
};

// Binary specialisation to skip some steps
// And just use (optimised) popcount
template <>
struct popcounter<1> {
  const size_t m_c;
  popcounter(size_t c=1) : m_c(c) {}

  template <typename block_t>
  inline size_t operator()(const block_t & x, size_t c) const {return popcount(x,c);}
  template <typename block_t>
  inline size_t operator()(const block_t & x) const {return popcount(x,m_c);}

  template <typename block_t>
  static inline block_t bitmap(const block_t & x, size_t c=1) {
    return (c)? x : ~x;
  }

  template <typename block_t>
  static inline size_t popcount(const block_t & x, size_t c=1) {
    const static size_t block_bits = bitwidth<block_t>::width;
    block_t temp = bitmap(x,c);
    return bitset<block_bits>(temp).count();
  }
};

template <size_t w=1, typename block_t>
inline size_t popcount(const block_t & x, size_t c) {
  return popcounter<w>::popcount(x, c);
}

template <size_t w=1, size_t c=1, typename block_t>
inline size_t popcount(const block_t & x) {
  return popcount(x, c);
}

// Over [0, i), like SDSL
template <size_t w=1, typename block_t>
inline size_t rank(const block_t & x, size_t i, size_t c=1) {
  if (i == 0) return 0;
  const static size_t block_bits  = bitwidth<block_t>::width;
  assert(i <= block_bits/w);
  block_t temp = popcounter<w>::bitmap(x, c);
  cout << bitset<block_bits>(temp) << "<" << endl;
  size_t shift = block_bits - i*w;
  temp <<= shift;
  return bitset<block_bits>(temp).count();
}

// TODO: write generalized vectorize instead
template <class popcounter>
struct vectorized_popcount {
  template <class InputIterator>
  inline size_t operator()(InputIterator start, InputIterator end, size_t c) const {
    typedef typename InputIterator::value_type block_t;
    auto f = [=](block_t x) -> size_t { return popcounter::popcount(x, c); };
    auto a = make_transform_iterator(start, f);
    auto b = make_transform_iterator(end, f);
    // TODO: detect gnu_parallel, thrust, etc
    // https://gcc.gnu.org/onlinedocs/libstdc++/manual/parallel_mode.html
    return __gnu_parallel::accumulate(a, b, 0);
  }

  template <class InputRange>
  inline size_t operator()(InputRange in, size_t c) const {
    return operator()(in.begin(), in.end(), c);
  }
};

template <size_t w, typename block_t>
block_t set_symbol(block_t & b, size_t i, size_t c) {
  block_t sym_mask = (1<<w) - 1;
  block_t mask = ~((sym_mask) << (i * w));
  b &= mask;
  b |= (c << (i*w));
  return b;
}

template <size_t w, typename block_t>
block_t set_symbol(const block_t & b, size_t i, size_t c) {
  block_t temp(b);
  return  set_symbol(temp, i, c);
}

template <size_t w, typename block>
int get_symbol(const block & x, size_t i) {
  block sym_mask = (1<<w) - 1;
  return (int)((x >> (i*w)) & sym_mask);
}

// Two layers of ranks. Large ranks are global, small ranks are local to large rank border.
// Individual blocks must be iteratively ranked in between small rank borders
template <size_t t_small_block_bits = cache_line_bits>
//          typename block_vector_type=vector>
//          typename rank_vector_type=vector>
class dna_rs {
private:
  // TODO: make these parameters to support
  // storing blocks on disk (using stxxl)
  // and ranks in ram (for example)
  const static size_t   block_bits = 64;
  const static size_t    dna_sigma = 8; // DNA + minus flags ($s stored separately)
  const static size_t     dna_bits = 3; // bits
  const static size_t  block_width = block_bits/dna_bits;

  typedef uint64_t             block_t;
  typedef uint16_t             small_rank_t;
  typedef uint64_t             large_rank_t;
  typedef vector<block_t>      block_vector;
  typedef sd_vector<>          terminal_bv;
  typedef vector<large_rank_t> large_rank_vector;
  typedef vector<small_rank_t> small_rank_vector;
  typedef typename terminal_bv::rank_1_type   terminal_rank_t;
  typedef typename terminal_bv::select_1_type terminal_select_t;

  const static size_t num_cache_lines = 1;
  const static size_t small_rank_bits       = bitwidth<small_rank_t>::width;
  // Make small blocks fit in cache line(s)
  // and align with data blocks (64 bits)
  const static size_t small_block_width     = (num_cache_lines * cache_line_bits)/block_width
                                            * block_width;
  // To make full use of the width of the small rank counts
  const static size_t large_block_max_width = 1<<small_rank_bits;
  // To align large blocks and small blocks
  const static size_t large_block_width     = large_block_max_width/small_block_width
                                            * small_block_width;

  // Data Members
  size_t            m_size; 
  block_vector      m_blocks;
  terminal_bv       m_terminals;
  terminal_rank_t   m_terminal_rank;
  terminal_select_t m_terminal_select;
  // Count vectors for each symbol.
  // Note: symbols dont change during a query, so
  // binary-search-based select should be faster
  // as related ranks are prefetched
  array<large_rank_vector, dna_sigma> m_large_ranks;
  array<small_rank_vector, dna_sigma> m_small_ranks;

public:
  template<typename InputIterator>
  dna_rs(InputIterator in, InputIterator end) {
    // Allocate space
    m_size = std::distance(in, end);
    m_blocks = block_vector(m_size, 0);
    // Large block boundaries at end of each large block
    size_t num_large_blocks = m_size/large_block_width;
    // Small block boundaries at end of each small block,
    // but not last (where large block boundary is)
    size_t num_small_blocks = m_size/small_block_width - num_large_blocks;

    bit_vector terminals(m_size,0);
    array<size_t, dna_sigma> global_counts;
    array<size_t, dna_sigma> local_counts;

    for (size_t i = 0 ; in != end; ++in, ++i) {
      int x = *in;
      assert(x >= 0 && x < dna_sigma + 1);
      // Mark $ signs in terminal bitvector
      if ( x == 0 ) {
        terminals[i] = 1;
      }
      // Offset the symbols to be one less
      // so we can index the counts properly
      else {
        --x;
        //global_counts[x]++;
        auto temp = m_blocks[i/block_width];
        m_blocks[i/block_width] = set_symbol<dna_bits>(temp, i%block_width, x);
      }
    }
    m_terminals = terminal_bv(terminals);
    sdsl::util::init_support(m_terminal_rank,   &m_terminals);
    sdsl::util::init_support(m_terminal_select, &m_terminals);
    // TODO: make members const?
  }

  template<typename InputRange>
  dna_rs(InputRange input) : dna_rs(input.begin(), input.end()) {}

  dna_rs(const std::string & input) : dna_rs(input | transformed(nt_to_int())) {}

  dna_rs(const char * input) : dna_rs(std::string(input)) {}

  char access(size_t i) const {
    int encoding = 0;
    if (1 == m_terminals[i]) encoding = 0;
    else {
      auto block = m_blocks[i/block_width];
      encoding = 1 + get_symbol<dna_bits>(block, i%block_width);
    }
    return int_to_nt()(encoding);
  }
 
  size_t rank(size_t i, char c) const {
    int x = nt_to_int()(c);
    if (x == 0) return m_terminal_rank(i);
    --x; // decrement symbol for vector representation
    auto f = vectorized_popcount<popcounter<dna_bits>>();
    // TODO: find start and end of small block instead
    size_t n_pre_blocks = i/block_width;
    cout << i << ", pre blocks:" << n_pre_blocks << endl;
    size_t pre_count = f(m_blocks.begin(), m_blocks.begin() + n_pre_blocks, x);
    cout << "precount: " << pre_count << endl;
    size_t r = cosmo::index::rank<dna_bits>(m_blocks[n_pre_blocks], i%block_width + 1, x);
    cout << "rank: " << r << endl;
    return pre_count + r;
  }

  size_t select(size_t i, char c) const {
    int x = nt_to_int()(c);
    if (x == 0) return m_terminal_select(i);
    // upper levels: gnu parallel search
    // GPU: (probably want wider small blocks?)
    // http://thrust.github.io/doc/group__vectorized__binary__search.html 
    // lower:
    // parallel partial sum/scan
    return 0;
  }

  // TODO: use table instead
  size_t count(char c) const {
    int x = nt_to_int()(c);
    if (x == 0) return m_terminal_rank(size()+1);
    --x;
    auto f = vectorized_popcount<popcounter<dna_bits>>();
    return f(m_blocks, x);
  }

  size_t size() const {
    return m_size;
  }
};

} // namespace index
} // namespace cosmo

#endif
