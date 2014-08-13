#ifndef IO_H
#define IO_H

#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

#include "dummies.hpp"
#include "kmer.hpp"
#include "debug.h"

static const size_t MAX_BITS_PER_KMER = 128;
static const size_t BUFFER_SIZE = 0x8000; // 32Kb data buffer

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

#define PACKED_WIDTH 5

inline
packed_edge pack_edge(uint8_t symbol, bool start_flag, bool end_flag) {
  return (symbol << 2) | (start_flag << 1) | (end_flag);
}

inline
void append_packed_edge(uint64_t & block, packed_edge edge) {
  block = (block << PACKED_WIDTH) | edge;
}

class PackedEdgeOutputer {
  static const size_t capacity = 64/PACKED_WIDTH;

  ofstream _ofs;
  uint64_t _buf = 0;
  size_t _len   = 0;
  vector<size_t> _counts = vector<size_t>(DNA_RADIX + 1, 0);
  bool closed = true;

  public:
  PackedEdgeOutputer() {}

  void open(string filename) {
    _ofs.open(filename, ios::out | ios::binary);
    closed = false;
  }

  ~PackedEdgeOutputer() {
    if (!closed) close();
    closed = true;
  }

  void flush() {
    if (_len > 0) {
      // Make it so its easy to access the ith member in a block uniformly
      _buf <<= ((capacity - _len) * PACKED_WIDTH);
      _ofs.write((char*)&_buf, sizeof(uint64_t));
    }
    //_ofs.flush(); // let _ofs handle this
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
    _ofs.write((char*)&accum[0], (DNA_RADIX+1) * sizeof(size_t));
  }

  public:
  void close() {
    if (closed) return;
    closed = true;
    flush();
    write_counts();
    _ofs.close();
  }

  template <typename kmer_t>
  // Might be nicer if this was a Functor with operator() instead
  void write(edge_tag tag, const kmer_t & x, const uint32_t k, bool first_start_node, bool first_end_node) {
    uint8_t f_sym = get_f(tag, x, k);
    uint8_t w_sym = get_w(tag, x);
    //printf("counts[%c]: %zu -> ", "$acgt"[f_sym], _counts[f_sym]);
    _counts[f_sym]++;
    //printf("%zu\n", _counts[f_sym]);

    packed_edge edge = pack_edge(w_sym, first_start_node, first_end_node);
    append_packed_edge(_buf, edge);

    if (++_len == capacity) flush();
  }
};

template <typename kmer_t>
uint8_t get_w(edge_tag tag, const kmer_t & x) {
  if (tag == out_dummy) return 0;
  else return get_edge_label(x) + 1;
}

template <typename kmer_t>
uint8_t get_f(edge_tag tag, const kmer_t & x, const uint32_t k) {
  uint8_t sym;
  if (tag == in_dummy && k == 1) {
    return 0; // in_dummies might have $ if only an edge label
  }
  else if (tag == out_dummy) {
    sym = get_nt(x, 1); // since out_dummies are shifted
  }
  else {
    sym = get_nt(x, 1); // 2nd last symbol
  }
  return sym+1;
}

#endif
