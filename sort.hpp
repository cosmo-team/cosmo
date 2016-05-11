#pragma once
#ifndef SORT_HPP
#define SORT_HPP
#include <type_traits>
#include <thread>
#include <set>
#include <stxxl.h>
#include <stxxl/bits/containers/sorter.h>
#include <stxxl/deque>

#include <boost/tuple/tuple.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>

#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "debruijn_graph_shifted.hpp"
#include "dummies.hpp"
#include "kmer.hpp"
#include "debug.hpp"
#include "config.hpp"
#include "utility.hpp"
#include "multi_bit_vector.hpp"

namespace cosmo {

namespace {

template <class dbg_t>
struct make_dbg {
  template <class t_bit_vector_type, class t_edge_vector_type, class label_type, class dummy_vector_type>
  dbg_t operator()(size_t k, const t_bit_vector_type & node_flags, const t_edge_vector_type & edges, const array<size_t, 5>& symbol_ends, const label_type& alphabet,
    const t_bit_vector_type & dummy_flags, const dummy_vector_type & dummies) {
    return dbg_t(k, node_flags, edges, symbol_ends, alphabet, dummy_flags, dummies);
  }
};

template <>
struct make_dbg<debruijn_graph_shifted<>> {
  template <class t_bit_vector_type, class t_edge_vector_type, class label_type, class dummy_vector_type>
  debruijn_graph_shifted<> operator()(size_t k, const t_bit_vector_type & node_flags, const t_edge_vector_type & edges, const array<size_t, 5>& symbol_ends, const label_type& alphabet,
    const t_bit_vector_type &, const dummy_vector_type &) {
    return debruijn_graph_shifted<>(k, node_flags, edges, symbol_ends, alphabet);
  }
};

// Should be unrolled and inlined, making good use of registers
template <typename kmer_t>
inline void kmer_shifts(const kmer_t & x, kmer_t * a, const size_t num_shifts) {
  //const size_t num_shifts = bitwidth<kmer_t>::width / 2;

  // NOTE: misses off final shift (which makes 0)
  for (size_t i = 1; i < num_shifts; ++i) {
    a[i] = (x << (2*i));
  }
}

template <typename kmer_t>
inline void kmer_shifts(const kmer_t & x, kmer_t * a) {
  const size_t num_shifts = bitwidth<kmer_t>::width / 2;
  return kmer_shifts(x, a, num_shifts);
}

template <typename T, typename iterator>
struct typed_iterator : iterator {
  typedef input_iterator_tag iterator_category;
  typedef T value_type;
  typedef size_t difference_type;
  typedef T& reference;
  typedef T* pointer;
};

// TODO: patch STXXL to have type traits and be copy/movable instead
template <typename T, typename iterator>
typed_iterator<T, iterator> & make_typed_iterator(const iterator & it) {
  return (typed_iterator<T, iterator>&) it;
}

} // namespace

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
    //return get_start_node(a) < get_start_node(b);
    auto a_node = get_start_node(a);
    auto b_node = get_start_node(b);
    if (a_node < b_node) return true;
    if (a_node == b_node) return (a < b);
    return false;
  }
  value_type min_value() const { return std::numeric_limits<value_type>::min(); }
  value_type max_value() const { return std::numeric_limits<value_type>::max(); }
};

template <typename dummy_t>
struct dummy_less {
  typedef dummy_t value_type;
  typedef typename std::remove_reference<decltype(get<0>(dummy_t()))>::type kmer_t;
  typedef typename std::remove_reference<decltype(get<1>(dummy_t()))>::type length_t;

  //node_less<kmer_t> n_less;
  bool operator() (const value_type & a, const value_type & b) const {
    //return a < b;//n_less(get<0>(a), get<0>(b)) && (get<1>(a) < get<1>(b));
    kmer_t   a_node = get_start_node(get<0>(a));
    kmer_t   b_node = get_start_node(get<0>(b));
    uint64_t a_edge = get_edge_label(get<0>(a));
    uint64_t b_edge = get_edge_label(get<0>(b));

    // If string equal then sort by length, else sort by string
    if (a_node == b_node) {
      if (get<1>(a) == get<1>(b)) {
        return a_edge < b_edge;
      }
      else return (get<1>(a) < get<1>(b));
    }
    else return (a_node < b_node);
  }
  value_type min_value() const { return dummy_t(); }
  value_type max_value() const { return dummy_t(std::numeric_limits<kmer_t>::max(),
                                                std::numeric_limits<length_t>::max()); }
};

template <typename record_t>
struct record_less {
  typedef record_t value_type;
  typedef typename std::remove_reference<decltype(get<0>(record_t()))>::type kmer_t;

  node_less<kmer_t> less;

  bool operator() (const value_type & a, const value_type & b) const {
    return less(get<0>(a), get<0>(b));
  }
  // payload is ignored
  value_type min_value() const { return record_t(less.min_value()); }
  value_type max_value() const { return record_t(less.max_value()); }
};

template <class A, class B>
auto cons(A a, B b) -> decltype(boost::tuples::cons<A,B>(a,b)) {
  return boost::tuples::cons<A, B>(a,b);
}

/*
template<class InputIt1, class InputIt2, class Visitor>
OutputIt find_dummies(InputIt1 first1, InputIt1 last1,
                      InputIt2 first2, InputIt2 last2,
                      Visitor visit) {
  // While we still have records in Table A
  while (first1 != last1) {
    // If Table B is depleted then
    if (first2 == last2) {
      // Rest of Table A not in Table B
      return std::copy(first1, last1, d_first);
    // If A_x not in B
    } else if (*first1 < *first2) {
      // Output incoming dummy
      *d_first++ = *first1++;
    // If B_x not in A
    } else if (*first2 < *first1) {
      // Output outgoing dummy
      *d_first++ = *first2;
    // A_x == B_x (hence is in both A and B - no dummy needed)
    } else {
      // Output record
      ++first1;
    }
    // Progress inner table (outer table is progressed in every above case)
    ++first2;
  }
  // Table A depleted, so rest of Table B not in Table A
  return std::copy(first2, last2, d_first);
}
*/

template <class t_dbg_t, typename t_kmer_t, typename ... t_payload_t>
struct dbg_builder {
  typedef t_dbg_t dbg_t;
  typedef t_kmer_t kmer_t;
  typedef boost::tuple<t_payload_t...> payload_t;
  typedef boost::tuple<kmer_t, t_payload_t...> record_t;

  struct output_t {
    edge_tag tag;
    kmer_t edge;
    decltype(record_t().get_tail()) payload;
    bool is_first_prefix;
    bool is_first_suffix;
    size_t k; // for variable order (when we later have dummy shifts)
    size_t lcs; // for variable order
  };

  typedef stxxl::vector<record_t, 1, stxxl::lru_pager<8>, block_size> record_vector_t;
  typedef stxxl::vector<kmer_t,   1, stxxl::lru_pager<8>, block_size>   kmer_vector_t;
  typedef stxxl::vector<dummy_t,  1, stxxl::lru_pager<8>, block_size>  dummy_vector_t;
  typedef node_less<kmer_t>       node_comparator_t;
  typedef record_less<record_t> record_comparator_t; // TODO: do i need both of these?
  typedef kmer_less<kmer_t>       edge_comparator_t;
  typedef dummy_less<dummy_t>    dummy_comparator_t;
  typedef stxxl::sorter<record_t, record_comparator_t, block_size> record_sorter_t;
  //typedef stxxl::sorter<kmer_t, node_comparator_t, block_size>     node_sorter_t;
  typedef stxxl::sorter<kmer_t,     edge_comparator_t, block_size>   edge_sorter_t;
  typedef stxxl::sorter<dummy_t,   dummy_comparator_t, block_size>  dummy_sorter_t;

  private:
  size_t k;
  size_t M;
  const bool   do_swap;
  const bool   shift_dummies;
  const bool   variable_order; // TODO: make this based off input dbg type
  string out_file_base;
  reverse_complement<kmer_t> rc;
  reverse_nt<kmer_t>         revnt;
  swap_gt<kmer_t>            swap;

  // Split memory in half so these can sort concurrently
  record_sorter_t record_sorter;
  edge_sorter_t     edge_sorter;

  public:
  template <typename parameters_t>
  dbg_builder(const parameters_t & parameters) : k(parameters.k), M(parameters.m), do_swap(parameters.swap),
                                                 shift_dummies(parameters.shift_dummies),
                                                 variable_order(parameters.variable_order),
                                                 rc(k),
                                                 record_sorter(record_comparator_t(), M/2),
                                                 edge_sorter(edge_comparator_t(),   M/2) {
    out_file_base = parameters.output_prefix + parameters.output_base; //sdsl::util::basename(parameters.input_filename);
  }

  void push(const record_t & in_rec) {
    auto edge    = in_rec.get_head();
    auto payload = in_rec.get_tail();

    kmer_t x = revnt(edge);
    if (do_swap) x = swap(x);
    kmer_t y = rc(x);

    record_sorter.push(cons(x, payload));
    record_sorter.push(cons(y, payload));
    edge_sorter.push(x);
    edge_sorter.push(y);
  }

  template <typename ... inputs_t>
  void push(inputs_t ... inputs) {
    push(record_t(inputs...));
  }

  template <class Visitor1, class Visitor2>
  dbg_t build(Visitor1 pre_merge, Visitor2 visit) {
    //static_assert(sizeof(typename InputIterator::value_type) == sizeof(record_t), "Record size is different to iterator size.");
    COSMO_LOG(info) << "Building from " << record_sorter.size() << " edges...";

    COSMO_LOG(trace) << "Merging runs...";
    record_vector_t kmers_a; // colex node
    kmer_vector_t   kmers_b; // colex edge

    // TODO: Work out way to not materialize table (use sorters directly? wrap stream in iterator)
    std::thread t1([&](){
      record_sorter.sort();
      kmers_a.resize(record_sorter.size());
      COSMO_LOG(trace) << "Materializing table A...";
      stxxl::stream::materialize(record_sorter, kmers_a.begin(), kmers_a.end());
      record_sorter.finish_clear();
    });
    std::thread t2([&](){
      edge_sorter.sort();
      kmers_b.resize(edge_sorter.size());
      COSMO_LOG(trace) << "Materializing table B...";
      stxxl::stream::materialize(edge_sorter, kmers_b.begin(), kmers_b.end());
      edge_sorter.finish_clear();
    });
    t1.join();
    t2.join();

    auto get_key([](auto x) -> kmer_t {
      return get<0>(x);
    });

    // Can write dummy labels and positions instead (for non variable order)
    // stxxl::syscall_file dum_file(basename(parameters.input_filename) + ".dummies", stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
    stxxl::deque<size_t> incoming_dummy_positions;
    stxxl::deque<kmer_t> incoming_dummies;
    dummy_vector_t incoming_dummies_v;

    // NOTE: need to sort incoming dummies if we want all shifts (currently needed for variable order)
    dummy_sorter_t incoming_dummy_sorter(dummy_comparator_t(), M);

    stxxl::deque<kmer_t> outgoing_dummies_q;
    edge_sorter_t outgoing_dummy_sorter(edge_comparator_t(), M);
    kmer_vector_t outgoing_dummies_v;

    // Detect incoming dummies
    // TODO: add flags to turn off adding revcomps and instead find unary dBG paths (ala bcalm)
    // TODO: when single stranded support is added (i.e. no need for revcomps)
    // then we should run both B-A and A-B to find outgoing and incoming dummy edges respectively
    // (as we wont be able to use the revcomp trick). Can be run in two threads.
    // TODO: OR, we can use set_symmetric_difference with each wrapped with a tag that is ignored during the comparison
    {
      COSMO_LOG(trace) << "Searching for nodes requiring dummy edges...";
      typename record_vector_t::bufreader_type a_reader(kmers_a);
      typename kmer_vector_t::bufreader_type   b_reader(kmers_b);

      auto a = boost::make_iterator_range(make_typed_iterator<record_t>(a_reader.begin()),
                                          make_typed_iterator<record_t>(a_reader.end()))
                                        | transformed(get_key);
      auto b = boost::make_iterator_range(make_typed_iterator<kmer_t>(b_reader.begin()),
                                          make_typed_iterator<kmer_t>(b_reader.end()));

      //kmer_t shifts[bitwidth<kmer_t>::width/2];
      //const size_t in_mem_shifts = 16;
      //const size_t max_shift = k - in_mem;
      size_t num_incoming_dummies = 0;
      if (shift_dummies) {
        // B - A
        // We can add reverse complements and sort them to get outgoing dummies
        // So the set difference direction doesn't matter
        /*
        typedef boost::dynamic_bitset<> bits;
        array<bits, in_mem> low_shifts;
        for (size_t i = 0; i < in_mem; ++i) {
          size_t num_shifts = 1<<(2*(i+1));
          low_shifts[i] = bits(num_shifts);
        }
        */

        kmer_t prev_kmer{};
        size_t idx = 0;
        find_outgoing_dummy_nodes<kmer_t>(a, b, k, [&](kmer_t x) {
          num_incoming_dummies++;
          outgoing_dummies_q.push_back(get_end_node(x,k)>>2); // add edge symbol (A) for (possibly) simpler merging
          kmer_t y = ((rc(x))<<2)>>2; // dummy node with A as outgoing edge (requires shifting)
          // kmer_shifts(y, shifts);
          size_t lcp_len = (idx++>0)*lcp(y, prev_kmer, k);
          prev_kmer = y;
          size_t this_max = k-lcp_len; //std::min(k - lcp_len, max_shift);
          //cerr << kmer_to_string(y, k) << " " << lcp_len << endl;
          for (size_t i = 1; i < this_max; ++i) {
            auto dum = y<<(2*i);
            //cerr << kmer_to_string(dum, k, k-i) << " (" << k-i << "/" << k-lcp_len << ")" << endl;
            incoming_dummy_sorter.push(dummy_t(dum, k-i));
          }
          // generate final shifts in memory if we need to
          /*
          for (size_t i = this_max; i < k - lcp_len; ++i) {
            size_t this_k = k - i;
            auto idx = y<<(2*i);
            auto edge = idx >> (bitwidth<kmer_t>::width - 2);
            idx = (idx << 2) >> (bitwidth<kmer_t>::width - this_k * 2);
            idx |= edge;
            //if (this_k > in_mem) {
              //COSMO_LOG(info) << "this_k, i, lcp_len, this_max : " << this_k << " " << i << " " << lcp_len << " " << this_max;
            //}
            //low_shifts[this_k-1][idx] = 1;
            auto original = (idx << (bitwidth<kmer_t>::width - 2)) | ((idx >> 2) << (bitwidth<kmer_t>::width - this_k * 2));
          }
          */
        });
      // TODO: producer consumer queues or streams (so we can start merging straight away in another thread)
      } else {
        // A - B
        // We do it this way to get indices of incoming dummies instead
        // Although it means we have to sort the outgoing dummies
        find_incoming_dummy_nodes<kmer_t>(a, b, k, [&](size_t idx, kmer_t x) {
          incoming_dummy_positions.push_back(idx);
          auto y = (rc(x<<2)); // dummy node with A as outgoing edge
          outgoing_dummy_sorter.push(y);
        });
      }

      kmers_b.clear();

      if (shift_dummies) {
        // Create shifts
        // Can sorting be skipped in favour of counting symbols in each position + radix sort-like approach?
        COSMO_LOG(info) << "# out dummies      : " << outgoing_dummies_q.size();
        COSMO_LOG(info) << "# inc dummies      : " << num_incoming_dummies;
        /*
        kmer_t prev_kmer{};
        size_t idx = 0;
        const size_t in_mem = 12; // all shifts but the last in_mem shifts
        const size_t max_shift = k - in_mem;
        for (auto x:incoming_dummies) {
          size_t lcp_len = (idx++>0)*lcp(x, prev_kmer, k);
          prev_kmer = x;
          // TODO: if I end up sorting just the dummies again
          // (no shifts i.e. if I work out how to generate them in order),
          // move this back up to the finding dummies loop, and add LCP field
          size_t this_max = std::min(k - lcp_len, max_shift);
          for (size_t i = 1; i < this_max; ++i) {
            incoming_dummy_sorter.push(dummy_t(x<<(2*i), k-i));
            //cerr << kmer_to_string(x, k) << " " << lcp_len << endl;
          }
          cerr << kmer_to_string(x, k) << " " << lcp_len << " " << this_max << endl;
          */
        //}
        //incoming_dummies.clear();
        COSMO_LOG(info) << "# inc dummy shifts : " << incoming_dummy_sorter.size();
        COSMO_LOG(trace) << "Sorting dummies...";
        incoming_dummy_sorter.sort();
        COSMO_LOG(trace) << "Materializing dummies...";
        incoming_dummies_v.resize(incoming_dummy_sorter.size());
        stxxl::stream::materialize(incoming_dummy_sorter, incoming_dummies_v.begin(), incoming_dummies_v.end());
        incoming_dummy_sorter.finish_clear();
        /*
        for (auto x : incoming_dummies_v) {
          cerr << kmer_to_string(get<0>(x), k, get<1>(x)) << endl;
        }
        */
      } else {
        COSMO_LOG(info) << "# out dummies: " << outgoing_dummy_sorter.size();
        COSMO_LOG(info) << "# inc dummies: " << incoming_dummy_positions.size();
        COSMO_LOG(trace) << "Sorting dummies...";
        outgoing_dummy_sorter.sort();
        COSMO_LOG(trace) << "Materializing dummies...";
        outgoing_dummies_v.resize(outgoing_dummy_sorter.size());
        stxxl::stream::materialize(outgoing_dummy_sorter, outgoing_dummies_v.begin(), outgoing_dummies_v.end());
        outgoing_dummy_sorter.finish_clear();
      }
    }

    /*
    cout << "outgoing dummies:" << endl;
    for (auto x : outgoing_dummies) {
      cout << kmer_to_string(x, k) << endl;
    }
    cout << endl;
    */
    // Call pre-merge visitor if it exists
    pre_merge(kmers_a.size() + outgoing_dummies_q.size() + outgoing_dummies_v.size() + incoming_dummies_v.size());

    typename record_vector_t::bufreader_type a_reader(kmers_a);
    typename kmer_vector_t::bufreader_type   o_reader(outgoing_dummies_v);
    typename dummy_vector_t::bufreader_type  i_reader(incoming_dummies_v);

    auto a = boost::make_iterator_range(make_typed_iterator<record_t>(a_reader.begin()),
                                        make_typed_iterator<record_t>(a_reader.end()));
    auto o = boost::make_iterator_range(make_typed_iterator<kmer_t>(o_reader.begin()),
                                        make_typed_iterator<kmer_t>(o_reader.end()));
    auto i = boost::make_iterator_range(make_typed_iterator<dummy_t>(i_reader.begin()),
                                        make_typed_iterator<dummy_t>(i_reader.end()));

    COSMO_LOG(trace) << "Merging in dummies...";

    using namespace sdsl;
    string temp_edge_file = out_file_base + ".w.temp";
    stxxl::syscall_file edge_file(temp_edge_file, stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT | stxxl::file::TRUNC);
    edge_file.set_size(kmers_a.size() + outgoing_dummies_q.size() + outgoing_dummies_v.size() + incoming_dummies_v.size());
    // using pointers so can call destructor later, which seems to stop the file from growing randomly and causing errors. A little hacky.
    stxxl::vector<uint8_t> * output = new stxxl::vector<uint8_t>(&edge_file);
    typename stxxl::vector<uint8_t>::bufwriter_type * edge_writer = new stxxl::vector<uint8_t>::bufwriter_type(*output);
    //int_vector<8> output(kmers_a.size() + outgoing_dummies_q.size() + outgoing_dummies_v.size() + incoming_dummies_v.size(), 0);
    // TODO: make dBG have stxxl vector for dummies instead
    //typename stxxl::vector<kmer_t>::bufwriter_type dummy_writer(dummies_ext);
    vector<kmer_t> dummies;
    dummies.reserve(incoming_dummy_positions.size());
    // TODO: look into swapping bit_vector's underlying vector with stxxl vector
    bit_vector * node_starts = new bit_vector(kmers_a.size() + outgoing_dummies_v.size() + outgoing_dummies_q.size() + incoming_dummies_v.size());
    bit_vector * dummy_flags = new bit_vector((!shift_dummies) * (kmers_a.size() + outgoing_dummies_v.size()));
    array<size_t, 5> counts{0,0,0,0,0};

    size_t idx = 0;
    kmer_t prev_edge{};
    uint8_t prev_k = 0;
    FirstStartNodeFlagger is_first_start_node(k);
    FirstEndNodeFlagger   is_first_end_node(k);


    auto merge_visitor = [&] (edge_tag tag, record_t rec, uint8_t this_k) {
      auto x = rec.get_head();
      auto payload = rec.get_tail();
      bool is_first_prefix = is_first_start_node(x, this_k);
      bool is_first_suffix = is_first_end_node(tag, x, this_k);
      size_t common_suffix_length = 0;
      if (variable_order && idx > 0) {
        common_suffix_length = node_lcs(x, prev_edge, std::min(this_k, prev_k));
        // If dummy edge length is equal, and the LCS length is the length of the non-dummy node suffix,
        // include the $ signs as well (not sure if this is needed)
        // TODO : test turning this off - functionality the same?
        if (prev_k == this_k && this_k == common_suffix_length + 1) {
          common_suffix_length = k - 1;
        }
      }
      output_t result{tag, x, payload, is_first_prefix, is_first_suffix, this_k, common_suffix_length};

      visit(result);

      char label = get_edge_label(x);
      int  node_last_sym = (tag == in_dummy && this_k == 1)? 0:1+get_nt(x, 1);
      bool is_out_dummy = (tag == out_dummy);
      // TODO: change to use ascii (consistently) coz these symbols are annoyiinnng
      char w_idx = is_out_dummy?0:((1+label)<<1) | (!result.is_first_suffix);
      //cerr << (int) w_idx << endl;
      //char w = is_out_dummy?'$':(DNA_ALPHA "ACGT")[w_idx];
      //output[idx]=w_idx;
      *edge_writer << w_idx;
      counts[node_last_sym]++;

      (*node_starts)[idx]   = !is_first_prefix; // inverted so sparse bv is sparser

      if (!shift_dummies) {
        if (result.tag == in_dummy) {
          dummies.push_back(x<<2); // remove edge symbol and add to vector
          //dummy_writer << get<0>(x);
        }
        (*dummy_flags)[idx] = (result.tag == in_dummy) || !is_first_prefix;
      }

      prev_edge = x;
      prev_k = this_k;
      idx++;
    };

    if (shift_dummies) {
      merge_dummies_with_shifts(a, outgoing_dummies_q, i, k, merge_visitor);
    } else {
      merge_dummies(a, o, incoming_dummy_positions, k, merge_visitor);
    }

    // Cleanup
    kmers_a.clear();
    outgoing_dummies_v.clear();
    incoming_dummy_positions.clear();

    // Prefix sum counts
    for (size_t i = 1; i < 5; ++i) {
      counts[i] += counts[i - 1];
    }

    delete edge_writer;
    delete output;

    // TODO: get these types from the input dBG instead
    COSMO_LOG(trace) << "Building dBG...";
    typedef wt_huff<rrr_vector<63>> wt_t;
    wt_t edges;
    construct(edges, temp_edge_file, 1);
    edge_file.close_remove();

    //construct_im(edges_mbv, output);
    // TODO: add parameter to keep temp files
    //boost::filesystem::remove(out_file_base+".edges");
    sd_vector<> node_bv(*node_starts);
    delete node_starts;
    sd_vector<> dum_pos_bv;
    if (!shift_dummies) {
      dum_pos_bv = sd_vector<>(*dummy_flags);
    }
    delete dummy_flags;

    //COSMO_LOG(info) << "size of WT: " << size_in_mega_bytes(edges) << " MB";
    //COSMO_LOG(info) << "size of MBV: " << size_in_mega_bytes(edges_mbv) << " MB";
    //COSMO_LOG(info) << "size of node BV: " << size_in_mega_bytes(node_bv) << " MB";

    /*
    if (!shift_dummies) {
      COSMO_LOG(info) << "size of dummy BV: " << size_in_mega_bytes(dum_pos_bv) << " MB";
      COSMO_LOG(info) << "size of dummy vec: " << size_in_mega_bytes(dummies) << " MB";
    }
    */

    return make_dbg<dbg_t>()(k, node_bv, edges, counts, "$ACGT", dum_pos_bv, dummies);
  }

  template <class Visitor>
  dbg_t build(Visitor visit) { return build([](auto){}, visit); }

  dbg_t build() { return build([](auto){}); }
};


} // namespace cosmo

#endif
