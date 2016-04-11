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
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>

#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "dummies.hpp"
#include "kmer.hpp"
#include "debug.hpp"
#include "config.hpp"
#include "utility.hpp"

namespace cosmo {

namespace {

// Should be unrolled and inlined, making good use of registers
template <typename kmer_t>
void kmer_shifts(const kmer_t & x, kmer_t * a) {
  size_t width = bitwidth<kmer_t>::width / 2;
  for (size_t i = 0; i < width; ++i) {
    a[i] = (x << (i*2));
  }
}

template <typename T, typename iterator>
struct typed_iterator : iterator {
  typedef input_iterator_tag iterator_category;
  typedef T value_type;
  typedef size_t difference_type;
  typedef T& reference;
  typedef T* pointer;
};

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
    return get_start_node(a) < get_start_node(b);
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
    return a < b;//n_less(get<0>(a), get<0>(b)) && (get<1>(a) < get<1>(b));
    /*
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
    */
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
  typedef stxxl::sorter<kmer_t,     edge_comparator_t, block_size>   edge_sorter_t;
  typedef stxxl::sorter<dummy_t,   dummy_comparator_t, block_size>  dummy_sorter_t;

  private:
  size_t k;
  size_t M;
  bool   do_swap;
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

    std::function<kmer_t(record_t)> record_key([](record_t x) -> kmer_t {
      return get<0>(x);
    });

    std::function<kmer_t(dummy_t)> dummy_key([](dummy_t x) -> kmer_t {
      return get<0>(x);
    });

    // Can write dummy labels and positions instead (for non variable order)
    // stxxl::syscall_file dum_file(basename(parameters.input_filename) + ".dummies", stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
    stxxl::deque<size_t> incoming_dummy_positions;

    // NOTE: need to sort incoming dummies if we want all shifts (currently needed for variable order)
    //dummy_vector_t incoming_dummies;
    //dummy_sorter_t incoming_dummy_sorter(dummy_comparator_t(), M);

    // Try with just queue
    stxxl::deque<kmer_t> outgoing_dummies_q;
    edge_sorter_t outgoing_dummy_sorter(edge_comparator_t(), M);
    kmer_vector_t outgoing_dummies;

    // Detect incoming dummies
    //std::thread t3([&]()
    {
      COSMO_LOG(trace) << "Searching for nodes requiring incoming dummy edges...";
      typename record_vector_t::bufreader_type a_reader(kmers_a);
      typename kmer_vector_t::bufreader_type   b_reader(kmers_b);

      auto a = boost::make_iterator_range(make_typed_iterator<record_t>(a_reader.begin()),
                                          make_typed_iterator<record_t>(a_reader.end()))
                                        | transformed(record_key);
      auto b = boost::make_iterator_range(make_typed_iterator<kmer_t>(b_reader.begin()),
                                          make_typed_iterator<kmer_t>(b_reader.end()));

      // TODO: use set_symmetric_difference with each wrapped with a tag that is ignored during the comparison
      // TODO: producer consumer queues or streams
      find_incoming_dummy_nodes<kmer_t>(a, b, k, [&](size_t idx, kmer_t x) {
        incoming_dummy_positions.push_back(idx);
        // We can add reverse complements and sort them to get outgoing dummies
        auto y = (rc(x<<2));
        outgoing_dummy_sorter.push(y);
      });

      kmers_b.clear();

      COSMO_LOG(trace) << "Sorting outgoing dummies...";
      outgoing_dummy_sorter.sort();
      outgoing_dummies.resize(outgoing_dummy_sorter.size());
      COSMO_LOG(trace) << "Materializing outgoing dummies...";
      stxxl::stream::materialize(outgoing_dummy_sorter, outgoing_dummies.begin(), outgoing_dummies.end());
      outgoing_dummy_sorter.finish_clear();
    }
    //);
    //t3.join();
    // TODO: add flags to turn off adding revcomps and instead find unary dBG paths (ala bcalm)
    // TODO: add flag to toggle outputting bitvector marking incoming dummies, or full dummy shifts
    COSMO_LOG(info) << "# inc dummies: " << incoming_dummy_positions.size();
    COSMO_LOG(info) << "# out dummies: " << outgoing_dummies.size();

    /*
    cout << "incoming dummies:" << endl;
    for (auto x : incoming_dummies) {
      cout << kmer_to_string(get<0>(x), k, get<1>(x)) << endl;
    }
    cout << endl;
    */

    /*
    cout << "outgoing dummies:" << endl;
    for (auto x : outgoing_dummies) {
      cout << kmer_to_string(x, k) << endl;
    }
    cout << endl;
    */
    pre_merge(kmers_a.size() + outgoing_dummies.size());

    typename record_vector_t::bufreader_type a_reader(kmers_a);
    typename kmer_vector_t::bufreader_type   o_reader(outgoing_dummies);
    //typename dummy_vector_t::bufreader_type  i_reader(incoming_dummies);

    auto a = boost::make_iterator_range(make_typed_iterator<record_t>(a_reader.begin()),
                                        make_typed_iterator<record_t>(a_reader.end()));
                                       //| transformed(record_key);
    auto o = boost::make_iterator_range(make_typed_iterator<kmer_t>(o_reader.begin()),
                                        make_typed_iterator<kmer_t>(o_reader.end()));
    //auto i = boost::make_iterator_range(make_typed_iterator<dummy_t>(i_reader.begin()),
    //                                   | transfonode dummy_key);

    COSMO_LOG(trace) << "Merging in dummies...";

    using namespace sdsl;
    int_vector<8> output(kmers_a.size() + outgoing_dummies.size(), 0);
    //stxxl::syscall_file dummy_file(out_file_base + ".dummies", stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
    // TODO: make dBG have stxxl vector for dummies instead
    //stxxl::vector<kmer_t> dummies_ext(&dummy_file);
    //dummies_ext.resize(incoming_dummy_positions.size());
    //typename stxxl::vector<kmer_t>::bufwriter_type dummy_writer(dummies_ext);
    vector<kmer_t> dummies;
    dummies.reserve(incoming_dummy_positions.size());
    // TODO: look into swapping bit_vector's underlying vector with stxxl vector
    bit_vector node_starts(kmers_a.size() + outgoing_dummies.size());
    bit_vector dummy_flags(kmers_a.size() + outgoing_dummies.size());
    array<size_t, 5> counts{0,0,0,0,0};

    size_t idx = 0;
    FirstStartNodeFlagger is_first_start_node(k);
    FirstEndNodeFlagger   is_first_end_node(k);
    merge_dummies(a, o, incoming_dummy_positions, [&](edge_tag tag, record_t rec){
      auto x = rec.get_head();
      auto payload = rec.get_tail();
      bool is_first_prefix = is_first_start_node(x, k);
      bool is_first_suffix = is_first_end_node(tag, x, k);
      output_t result{tag, x, payload, is_first_prefix, is_first_suffix};

      visit(result);

      char label = get_edge_label(x);
      int  node_last_sym = get_nt(x, 1);
      bool is_out_dummy = (tag == out_dummy);
      char w_idx = is_out_dummy?0:label | ((!result.is_first_suffix)<<2);
      //char w = is_out_dummy?'$':(DNA_ALPHA "ACGT")[w_idx];
      output[idx]=w_idx;
      counts[1 + node_last_sym]++;
      if (result.tag == in_dummy) {
        dummies.push_back(x);
        //dummy_writer << get<0>(x);
      }
      node_starts[idx]   = !is_first_prefix; // inverted so sparse bv is sparser
      dummy_flags[idx] = (result.tag == in_dummy);
      idx++;
    });

    // Cleanup
    kmers_a.clear();
    outgoing_dummies.clear();
    incoming_dummy_positions.clear();

    // TODO: get these types from the input dBG instead
    COSMO_LOG(trace) << "Building dBG...";
    typedef wt_huff<rrr_vector<63>> wt_t;
    wt_t edges;
    //construct(edges, out_file_base+".edges", 1);
    construct_im(edges, output);
    // TODO: add parameter to keep temp files
    //boost::filesystem::remove(out_file_base+".edges");
    sd_vector<> node_bv(node_starts);
    sd_vector<> dum_pos_bv(dummy_flags);

    // Prefix sum counts
    counts[0] = 0;
    for (size_t i = 1; i < 5; ++i) {
      counts[i] += counts[i - 1];
    }

    COSMO_LOG(info) << "size of WT: " << size_in_mega_bytes(edges) << " MB";
    COSMO_LOG(info) << "size of node BV: " << size_in_mega_bytes(node_bv) << " MB";
    COSMO_LOG(info) << "size of dummy BV: " << size_in_mega_bytes(dum_pos_bv) << " MB";
    COSMO_LOG(info) << "size of dummy vec: " << size_in_mega_bytes(dummies) << " MB";

    dbg_t graph(k, node_bv, edges, counts, "$acgtACGT", dum_pos_bv, dummies);
    COSMO_LOG(info) << "size of DBG: " << size_in_mega_bytes(graph) << " MB";
    return graph;
  }

  template <class Visitor>
  dbg_t build(Visitor visit) { return build([](auto x){}, visit); }

  dbg_t build() { return build([](auto x){}); }
};


} // namespace cosmo

#endif
