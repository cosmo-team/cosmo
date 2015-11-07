#ifndef IO_HPP
#define IO_HPP

#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <fstream>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>
#include "utility.hpp"
#include "dummies.hpp"
#include "kmer.hpp"
#include "debug.hpp"

static const size_t MAX_BITS_PER_KMER = 128;
static const size_t BUFFER_SIZE = 1024 * 1024;

using namespace std;

#define DSK_FILE_RECORD_SIZE(NBITS) (((NBITS)/8) + 4)

// Reads the DSK file input header
int dsk_read_header(int, uint32_t *, uint32_t *);
// Counts the number of records in the file - for allocation purposes
int dsk_num_records(int handle, uint32_t kmer_num_bits, size_t * num_records);
// Read kmers from file into the output array
size_t dsk_read_kmers(int handle, uint32_t kmer_num_bits, uint64_t * kmers_output);
//void merge_and_output(FILE * outfile, uint64_t * table_a, uint64_t * table_b, uint64_t * incoming_dummies, size_t num_records, size_t num_incoming_dummies, uint32_t k);

typedef uint8_t packed_edge;

#define PACKED_WIDTH (5)
#define PACKED_CAPACITY (bitwidth<uint64_t>::width/PACKED_WIDTH)

inline
packed_edge pack_edge(uint8_t symbol, bool start_flag, bool end_flag) {
  return (symbol << 2) | (start_flag << 1) | (end_flag);
}

inline
void append_packed_edge(uint64_t & block, packed_edge edge) {
  block = (block << PACKED_WIDTH) | edge;
}

class PackedEdgeOutputer {
  public:
  static const size_t capacity = PACKED_CAPACITY;

  ostream & _os;

  uint64_t _buf = 0;
  size_t _len   = 0;
  vector<size_t> _counts;
  bool closed = false;

  PackedEdgeOutputer(ostream & os) : _os(os) {
    _counts = vector<size_t>(DNA_RADIX + 1, 0);
  }

  ~PackedEdgeOutputer() {
    close();
  }

  void flush() {
    if (_len > 0) {
      // Make it so its easy to access the ith member in a block uniformly
      _buf <<= ((capacity - _len) * PACKED_WIDTH);
      _os.write((char*)&_buf, sizeof(uint64_t));
    }
    //_os.flush(); // let _os handle this
    _buf = 0;
    _len = 0;
  }

  private:
  void write_counts() {
    vector<size_t> accum(_counts);
    for (int i = 1; i < DNA_RADIX+1; i++) {
      accum[i] += accum[i-1]; // accumulate
    }
    // write out counts
    _os.write((char*)&accum[0], (DNA_RADIX+1) * sizeof(size_t));
  }

  public:
  void close() {
    if (closed) return;
    flush();
    write_counts();
    closed = true;
    //_os.close();
  }

  template <typename kmer_t>
  // Might be nicer if this was a Functor with operator() instead
  // but then need copy ctor etc
  void write(edge_tag tag, const kmer_t & x, const uint32_t k, bool first_start_node, bool first_end_node) {
    uint8_t f_sym = get_f(tag, x, k);
    uint8_t w_sym = get_w(tag, x);
    //printf("counts[%c]: %zu -> ", "$acgt"[f_sym], _counts[f_sym]);
    _counts[f_sym]++;
    //printf("%zu\n", _counts[f_sym]);

    packed_edge edge = pack_edge(w_sym, first_start_node, first_end_node);
    //cout << ((edge & 2) >> 1) << endl;
    append_packed_edge(_buf, edge);

    if (++_len == capacity) flush();
  }
};

inline packed_edge get_packed_edge_from_block(uint64_t block, size_t i) {
  return (block >> (PACKED_WIDTH * (PACKED_CAPACITY-i-1))) & 31;
}

template <typename BlockIterator>
inline packed_edge get_packed_edge(BlockIterator blocks, size_t i) {
  size_t block_idx = i/PACKED_CAPACITY;
  size_t local_idx = i%PACKED_CAPACITY;
  return get_packed_edge_from_block(blocks[block_idx], local_idx);
}

typedef std::tuple<uint8_t, bool, bool> edge_tuple;

static inline uint8_t unpack_symbol(packed_edge x) { return x >> 2; }
static inline bool    unpack_node_flag(packed_edge x) { return ((x & 2) >> 1); }
static inline bool    unpack_edge_flag(packed_edge x) { return (x & 1); }

inline edge_tuple unpack_to_tuple(packed_edge x) {
  return edge_tuple(unpack_symbol(x), unpack_node_flag(x), unpack_edge_flag(x));
}

template <typename BlockIterator>
inline edge_tuple get_edge(BlockIterator blocks, size_t i) {
  return unpack_to_tuple(get_packed_edge(blocks, i));
}

#endif
