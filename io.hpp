#ifndef IO_HPP
#define IO_HPP

#include <algorithm>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <fstream>
#include <bitset>
#include <queue>
#include <sstream>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>

#include "utility.hpp"
#include "dummies.hpp"
#include "config.hpp"
#include "kmer.hpp"
#include "debug.hpp"
#include "kmc_api/kmc_file.h"

static const size_t MAX_BITS_PER_KMER = 128;
static const size_t BUFFER_SIZE = 1024 * 1024;

using namespace std;


namespace {

static inline uint64_t nibblet_reverse(const uint64_t &word)
{

    uint64_t ret_val = word;
    for (int i=0; i < 16; ++i) {
        //std::cout <<"************* "<< i << " ************" << std::endl << std::endl;
        //print_kmers(std::cout, &word , 1, 32);

        uint64_t right_mask = 3ull << (2*i);
        //print_kmers(std::cout, &right_mask , 1, 32);
        
        uint64_t left_mask = (3ull << 62ull) >> (2*i);
        //print_kmers(std::cout, &left_mask , 1, 32);
        
        uint64_t left_val = (word & left_mask) >> (62ull - (2*i));
        //print_kmers(std::cout, &left_val , 1, 32);
        
        uint64_t right_val = (word & right_mask) >> (2*i);
        //print_kmers(std::cout, &right_val , 1, 32);
        
        uint64_t new_left_bits = (right_val << 62ull) >> (2*i);
        //print_kmers(std::cout, &new_left_bits , 1, 32);
        
        uint64_t new_right_bits = left_val << (2*i);
        //print_kmers(std::cout, &new_right_bits , 1, 32);
        
        ret_val = (ret_val & ~left_mask) | new_left_bits;
        ret_val = (ret_val & ~right_mask) | new_right_bits;
        //std::cout << std::endl;
    }
    return ret_val;
}

inline void clear_bv(color_bv &bv) {
  bv.reset();
}

inline void set_bit(color_bv &bv, uint32_t j) {
  bv[j] = 1; // |= 1LL << j % 64;
}

// code from http://www.cplusplus.com/reference/queue/priority_queue/priority_queue/
// with the polarity reversed
// FIXME: don't compare strings!!!
typedef std::tuple<unsigned, CKmerAPI, uint64 > queue_entry;    

class mylessthan
{
    bool reverse;
public:
    mylessthan(const bool& revparam=false)
        {reverse=revparam;}
    inline bool operator() (const queue_entry& lhs, const queue_entry&rhs) const
        {
            // if (reverse) return (const_cast<CKmerAPI*>(&(lhs.second))->to_string() > const_cast<CKmerAPI*>(&(rhs.second))->to_string());
            // else return (const_cast<CKmerAPI*>(&(lhs.second))->to_string() < const_cast<CKmerAPI*>(&(rhs.second))->to_string());
            if (reverse) return ( *const_cast<CKmerAPI*>(&(std::get<1>(rhs))) < *const_cast<CKmerAPI*>(&(std::get<1>(lhs))));
            else return (*const_cast<CKmerAPI*>(&(std::get<1>(lhs))) < *const_cast<CKmerAPI*>(&(std::get<1>(rhs))));
        }
};


typedef std::priority_queue<queue_entry, std::vector<queue_entry>, mylessthan> mypq_type;

static inline int push(mypq_type& queue, const std::vector<CKMCFile *>& kmer_data_bases, const unsigned i, const unsigned k)
{
    int num_pushed = 0;
    CKmerAPI kmer_object(k);
    uint64 counter = 0;// for coverage
    if (kmer_data_bases[i]->ReadNextKmer(kmer_object, counter)) {
        queue_entry entry = std::make_tuple(i, kmer_object, counter);
        queue.push(entry);
        //std::cout << "queue.push" << print_entry(entry) << std::endl;
        num_pushed += 1;
    }
    return num_pushed;

}

static inline queue_entry pop_replace(mypq_type& queue, const std::vector<CKMCFile *>& kmer_data_bases, const unsigned k)
{
    queue_entry popped_value = queue.top();
    queue.pop();

    push(queue, kmer_data_bases, std::get<0>(popped_value), k);
    return popped_value;
}

}

inline bool kmc_read_header(std::string db_fname, uint32_t & k, size_t &min_union, size_t &max_union, size_t &num_colors, std::vector<CKMCFile *> &kmer_data_bases) {
    std::ifstream db_list(db_fname.c_str());
    std::string fname;
    unsigned colornum = 0;
    min_union = 0;
    max_union = 0;

    while ( db_list >>  fname ) {
        COSMO_LOG(info) << "Color " << colornum << ": " << fname;
        colornum++;
        
        CKMCFile * kmer_data_base = new CKMCFile;
        kmer_data_bases.push_back(kmer_data_base);

        if (!kmer_data_base->OpenForListing(fname))
        {
            COSMO_LOG(error) << "ERROR: Could not open KMC2 database '" << fname << "'";
            return false;
        }
        else
        {
            //TODO : check if the file is sorted FIXME
            
            //uint32 _kmer_length;
            uint32 _mode;
            uint32 _counter_size;
            uint32 _lut_prefix_length;
            uint32 _signature_len;
            uint32 _min_count;
            uint64 _max_count;
            uint64 _total_kmers;

            kmer_data_base->Info(k, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
            //kmer_num_bits = 64 * ((k * 2 - 1) / 64 + 1); // FIXME: double check ceil(quotients) is what we're getting and not more, this may be too conservative
            if (_total_kmers > min_union) {
              min_union = _total_kmers;
            }
            max_union += _total_kmers;
            //std::string str;
            COSMO_LOG(info) << "k: " << k << " max kmer count: " << _max_count << " total kmers: " << _total_kmers << std::endl;
        }
    }
    assert(kmer_data_bases.size() > 0);
    num_colors = kmer_data_bases.size();
    return true;
        
}

void update_multiplicity(const queue_entry& current, std::vector<uint64>& max_multiplicity, std::vector<uint64>& multiplicity_histogram)
{
    unsigned count = std::get<2>(current);
    unsigned queuenum = std::get<0>(current);
    if (max_multiplicity[queuenum] < count) max_multiplicity[queuenum] = count;
    unsigned bin = count ;
    if (bin < 256) {
        multiplicity_histogram[bin] += 1;
    }
}



template <class Visitor>
size_t kmc_read_kmers(std::vector<CKMCFile *> &kmer_data_bases, uint32_t k, Visitor visit) {

    // keep some stats per kmer database. FIXME: This is somewhat redundant with header info we report above
    std::vector<uint64> max_multiplicity(kmer_data_bases.size());
    std::vector<uint64> kmers_read(kmer_data_bases.size());
    std::vector<uint64> total_kmers(kmer_data_bases.size());
    std::vector<uint64> multiplicity_histogram(256);
    unsigned long long kmer_token_count = 0;

    for (unsigned int i = 0; i < multiplicity_histogram.size(); ++i) {
        multiplicity_histogram[i] = 0;
    }
    
    for (unsigned int i = 0; i < max_multiplicity.size(); ++i) {
        max_multiplicity[i] = 0;
        kmers_read[i] = 0;
        uint32 _mode;
        uint32 _counter_size;
        uint32 _lut_prefix_length;
        uint32 _signature_len;
        uint32 _min_count;
        uint64 _max_count;
        uint64 _total_kmers;
        kmer_data_bases[i]->Info(k, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
        total_kmers[i] = _total_kmers;
    }
                                           
    
    const mylessthan gt_comparitor(true);
    const mylessthan lt_comparitor(false);

    mypq_type queue(gt_comparitor);

    color_bv color = 0;
    
    // initialize the queue with a file identifier (as a proxy for the input sequence itself) and the value at the head of the file for peeking
    for (unsigned i = 0; i < kmer_data_bases.size(); ++i) {
        if (push(queue, kmer_data_bases, i, k) == 0) {
            std::cerr << "WARNING: File number " << i << " contains no k-mers." << std::endl;
        }
    }

    // pop the first element into 'current' to initialize our state (and init any other state here such as this one's color)
    queue_entry current = pop_replace(queue, kmer_data_bases, k);
    update_multiplicity(current, max_multiplicity, multiplicity_histogram);
    kmer_token_count += std::get<2>(current);
    kmers_read[std::get<0>(current)]++;
    
    color.set(std::get<0>(current)); // FIXME: make sure not using << operator elsewhere!
    //std::cout << "current = " << print_entry(current) << " = queue.pop()" << std::endl;    

    //
    unsigned long long num_merged_kmers = 0;
    while (!queue.empty()) {
        // std::string s1 = const_cast<CKmerAPI*>(&(queue.top().second))->to_string();
        // std::string s2 = const_cast<CKmerAPI*>(&(current.second))->to_string();
        // if (s1 == s2) { // if this is the same kmer we've seen before
        if (*const_cast<CKmerAPI*>(&(std::get<1>(queue.top()))) == *const_cast<CKmerAPI*>(&(std::get<1>(current)))) { // if this is the same kmer we've seen before
            queue_entry additional_instance = pop_replace(queue, kmer_data_bases, k);
            update_multiplicity(additional_instance, max_multiplicity, multiplicity_histogram);
            kmer_token_count += std::get<2>(additional_instance);
            kmers_read[std::get<0>(additional_instance)]++;
            color.set(std::get<0>(additional_instance));
            //std::cout << "additional_instance = " << print_entry(additional_instance) << " = queue.pop()" << std::endl;
        } else { // if the top of the queue contains a new instance

            // emit our current state
            std::vector<unsigned long long /*uint64*/> kmer;            
            std::get<1>(current).to_long(kmer);
            num_merged_kmers++;
            if (num_merged_kmers % 1000000 == 0) {
                uint64 global_max_mult = 0;
                float progress_sum = 0.0;
                for (unsigned int i = 0; i < max_multiplicity.size(); ++i) {
                    if (max_multiplicity[i] > global_max_mult) global_max_mult = max_multiplicity[i];
                    if (total_kmers[i] > 0)
                        progress_sum += (float)(kmers_read[i]) / (float)(total_kmers[i]);
                }

                float progress = progress_sum / max_multiplicity.size();
                std::cout << "Number of merged k-mers: " << num_merged_kmers
                          << " Max multiplicity seen: " << global_max_mult
                          << " fraction complete: " << progress  
                          << " Estimated total: " << (unsigned long long) (num_merged_kmers / progress)
                          << std::endl;
                std::cout << "Global Mult hist: [";
                for (unsigned int i = 0; i < multiplicity_histogram.size(); ++i) {
                    std::cout << multiplicity_histogram[i]  << ", ";
                }
                std::cout << "]" << std::endl;


            }

            //std::cout << const_cast<CKmerAPI*>(&(current.second))->to_string() << " : " << color << std::endl;
            //FIXME: the following line is to fix a bug in the kmc2 API where it shifts word[0] << 64 when k=63 which results in "word[1] =  word[0] + word[1]" the next line compensates; not sure how pervasive this is in the k>32 space.
            // if (kmer.size() >= 2) {
            //   kmer[1] -= kmer[0];
            // }

            std::reverse(kmer.begin(), kmer.end());
            visit(*((kmer_t*)&kmer[0]), color);
            color.reset();


            // now initialize our current state with the top
            current = pop_replace(queue, kmer_data_bases, k);
            kmer_token_count += std::get<2>(current);
            kmers_read[std::get<0>(current)]++;
            update_multiplicity(current, max_multiplicity, multiplicity_histogram);
            //std::cout << "current = " << print_entry(current) << " = queue.pop()" << std::endl;

            color.set(std::get<0>(current)); // FIXME: make sure not using << operator elsewhere!
        }
    }

    // and finally emit our current state    
    std::vector<unsigned long long /*uint64*/> kmer;            
    std::get<1>(current).to_long(kmer);
        

    //FIXME: the following line is to fix a bug in the kmc2 API where it shifts word[0] << 64 when k=63 which results in "word[1] =  word[0] + word[1]" the next line compensates; not sure how pervasive this is in the k>32 space.
    // if (kmer.size() >= 2) {
    //   kmer[1] -= kmer[0];
    // }
    std::reverse(kmer.begin(), kmer.end());
    visit(*((kmer_t*)&kmer[0]), color);
    color.reset();

    num_merged_kmers++;
            
    

		
    // CKmerAPI kmer_object(k);



        
    // while (kmer_data_bases[0]->ReadNextKmer(kmer_object, counter)) {


    //     std::vector<unsigned long long /*uint64*/> kmer;            
    //     kmer_object.to_long(kmer);


    //     color_bv color = 1;
    //     kmer_colors[numkmers] = color;

    //     for (int block=0; block < kmer.size(); ++block) {
    //         kmers_output[numkmers*kmer.size() + block] = kmer[block]; // FIXME: check if kmer_output is big endian or little endian

    //         if (block == 1) {
    //             assert(k == 63); //FIXME: the following line is to fix a bug in the kmc2 API where it shifts word[0] << 64 when k=63 which results in "word[1] =  word[0] + word[1]" the next line compensates; not sure how pervasive this is in the k>32 space.
    //             kmers_output[numkmers*kmer.size() + block] -=kmer[0];
    //         }

    //     }
    //     numkmers++;
    // }

    for (auto kmer_data_base: kmer_data_bases) {
        kmer_data_base->Close();
        delete kmer_data_base;
    }
    std::cerr << "Total k-mer token count: " << kmer_token_count << std::endl;
    return num_merged_kmers;
    
}
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
