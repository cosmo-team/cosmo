#pragma once
#ifndef SORT_HPP
#define SORT_HPP

#include <type_traits>
#include <thread>
#include <set>
#include <stxxl.h>
#include <stxxl/bits/containers/sorter.h>

#include <boost/tuple/tuple.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>

#include "dummies.hpp"
#include "kmer.hpp"
#include "debug.hpp"
#include "config.hpp"

namespace cosmo {

namespace {

// Should be unrolled and inlined, making good use of registers
template <typename kmer_t>
void kmer_shifts(kmer_t & x, kmer_t * a) {
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

template <typename t_kmer_t, typename ... t_payload_t>
struct kmer_sorter {
  typedef t_kmer_t kmer_t;
  typedef boost::tuple<t_payload_t...> payload_t;
  typedef boost::tuple<kmer_t, t_payload_t...> record_t;

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

  // take InputIterator for file whose value_type is record_type
  // the first element of the record tuple should be the kmer
  // Also accepts pure kmer_t input
  template <class InputRange, class Visitor, typename parameters_t>
  void sort(InputRange & s, const parameters_t & parameters, Visitor visit) const {
    //static_assert(sizeof(typename InputIterator::value_type) == sizeof(record_t), "Record size is different to iterator size.");
    size_t k = parameters.k;
    size_t M = parameters.m;

    auto rc    = reverse_complement<kmer_t>(k);
    auto revnt = reverse_nt<kmer_t>();
    auto swap  = swap_gt<kmer_t>();

    // Split memory in half so these can sort concurrently
    record_sorter_t record_sorter(record_comparator_t(), M/2);
    edge_sorter_t     edge_sorter(edge_comparator_t(),   M/2);

    //output record_types in sorted order with dummies merged
    COSMO_LOG(trace) << "Creating runs...";
    for ( auto in_rec : s ) {
      auto edge    = in_rec.get_head();
      auto payload = in_rec.get_tail();

      kmer_t x = revnt(edge);
      if (parameters.swap) x = swap(x);
      kmer_t y = rc(x);

      record_sorter.push(cons(x, payload));
      record_sorter.push(cons(y, payload));
      edge_sorter.push(x);
      edge_sorter.push(y);
    }
    COSMO_LOG(info) << "Added " << record_sorter.size()/2 << " edges, not including revcomps.";

    COSMO_LOG(trace) << "Merging runs...";
    record_vector_t kmers_a; // colex node
    kmer_vector_t   kmers_b; // colex edge

    // TODO: Test using sorters directly (wrap in stream then iterator) instead of materializing
    // to avoid IO
    std::thread t1([&](){
      record_sorter.sort();
      kmers_a.resize(record_sorter.size());
      COSMO_LOG(trace) << "Writing table A to temporary storage...";
      stxxl::stream::materialize(record_sorter, kmers_a.begin(), kmers_a.end());
      record_sorter.finish_clear();
    });
    std::thread t2([&](){
      edge_sorter.sort();
      kmers_b.resize(edge_sorter.size());
      COSMO_LOG(trace) << "Writing table B to temporary storage...";
      stxxl::stream::materialize(edge_sorter, kmers_b.begin(), kmers_b.end());
      edge_sorter.finish_clear();
    });
    t1.join();
    t2.join();

    cout << "A:" << endl;
    for (auto x : kmers_a) {
      cout << kmer_to_string(x, k) << endl;
    }
    cout << endl;

    cout << "b:" << endl;
    for (auto x : kmers_b) {
      cout << kmer_to_string(x, k) << endl;
    }
    cout << endl;

    exit(1);
    std::function<kmer_t(record_t)> record_key([](record_t x) -> kmer_t {
      return get<0>(x);
    });

    std::function<kmer_t(dummy_t)> dummy_key([](dummy_t x) -> kmer_t {
      return get<0>(x);
    });

    // Can write dummy labels and positions instead (for non variable order)
    //stxxl::syscall_file dum_file(basename(parameters.input_filename) + ".dummies", stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
    //stxxl::vector<size_t> dummy_positions;
    //stxxl::vector<size_t>::bufwriter_type position_writer(dummy_positions);
    kmer_vector_t incoming_dummies;
    kmer_vector_t outgoing_dummies;
    size_t num_incoming_dummies = 0;
    size_t num_outgoing_dummies = 0;

    /*
    std::thread t3([&](){
      COSMO_LOG(trace) << "Searching for nodes requiring incoming dummy edges...";
      typename kmer_vector_t::bufwriter_type   writer(incoming_dummies);
      typename record_vector_t::bufreader_type a_reader(kmers_a);
      typename kmer_vector_t::bufreader_type   b_reader(kmers_b);

      auto a = boost::make_iterator_range(make_typed_iterator<record_t>(a_reader.begin()),
                                          make_typed_iterator<record_t>(a_reader.end()))
                                        | transformed(record_key);
      auto b = boost::make_iterator_range(make_typed_iterator<kmer_t>(b_reader.begin()),
                                          make_typed_iterator<kmer_t>(b_reader.end()));

      find_incoming_dummy_nodes<kmer_t>(a, b, k, [&](size_t idx, kmer_t x) {
        //position_writer << idx;
        num_incoming_dummies++;
        writer << x;
      });
      //dummy_positions.resize(num_incoming_dummies);
      incoming_dummies.resize(num_incoming_dummies);
      COSMO_LOG(info) << "num_incoming_dummies: " << num_incoming_dummies;
    });
    */

    std::thread t4([&](){
      COSMO_LOG(trace) << "Searching for nodes requiring outgoing dummy edges...";
      typename kmer_vector_t::bufwriter_type   writer(outgoing_dummies);
      typename record_vector_t::bufreader_type a_reader(kmers_a);
      typename kmer_vector_t::bufreader_type   b_reader(kmers_b);

      auto a = boost::make_iterator_range(make_typed_iterator<record_t>(a_reader.begin()),
                                          make_typed_iterator<record_t>(a_reader.end()))
                                         | transformed(record_key);
      auto b = boost::make_iterator_range(make_typed_iterator<kmer_t>(b_reader.begin()),
                                          make_typed_iterator<kmer_t>(b_reader.end()));

      //kmer_t shifts[bitwidth<kmer_t>::width/2];
      //const size_t low_k = 16;
      //#define sum_pow_4(i) (((1<<(2*((i)+2)))-4)/3)
      find_outgoing_dummy_nodes<kmer_t>(a, b, k, [&](kmer_t x) {
        // TODO: generate all shifts
        num_outgoing_dummies++;
        if (x == 0) {
          COSMO_LOG(error) << "All A KMER WHYYYY";
        }
        writer << x;
      });
      outgoing_dummies.resize(num_outgoing_dummies);
      COSMO_LOG(info) << "num_outgoing_dummies: " << num_outgoing_dummies;
    });
    //t3.join();
    t4.join();

    /*
    cout << "incoming dummies:" << endl;
    for (auto x : incoming_dummies) {
      cout << kmer_to_string(rc(x), k) << endl;
    }
    cout << endl;
    */

    cout << "outgoing dummies:" << endl;
    for (auto x : outgoing_dummies) {
      cout << kmer_to_string(x, k) << endl;
    }
    cout << endl;

    // second table no longer needed
    kmers_b.clear();

    //COSMO_LOG(trace) << "Sorting incoming dummies...";
    //dummy_sorter.sort();

    //COSMO_LOG(trace) << "Writing incoming dummies...";
    //incoming_dummies.resize(dummy_sorter.size());
    //stxxl::stream::materialize(dummy_sorter, incoming_dummies.begin(), incoming_dummies.end());
    //dummy_sorter.finish_clear();

    //auto payload = payload_t();
    //std::function<kmer_t(record_t)> capture_payload([&](record_t x){
      //payload(x.get_tail()); // might be null_type if only kmers provided, but thats ok
      //return get<0>(x);
    //});

    /*
    typename record_vector_t::bufreader_type a_reader(kmers_a);
    typename kmer_vector_t::bufreader_type   o_reader(outgoing_dummies);
    typename dummy_vector_t::bufreader_type  i_reader(incoming_dummies);

    auto a = boost::make_iterator_range(make_typed_iterator<record_t>(a_reader.begin()),
                                        make_typed_iterator<record_t>(a_reader.end()))
                                       | transformed(record_key);
    auto o = boost::make_iterator_range(make_typed_iterator<kmer_t>(o_reader.begin()),
                                        make_typed_iterator<kmer_t>(o_reader.end()));
    auto i = boost::make_iterator_range(make_typed_iterator<dummy_t>(i_reader.begin()),
                                        make_typed_iterator<dummy_t>(i_reader.end()))
                                       | transformed(dummy_key);

    COSMO_LOG(trace) << "Merging and writing files...";
    //stxxl::vector<kmer_t> chars;
    stxxl::vector<char> output;
    output.resize(kmers_a.size() + outgoing_dummies.size() + incoming_dummies.size());
    typename stxxl::vector<char>::bufwriter_type w(output);

    merge_dummies(a, o, i, [&](kmer_t x){
      w << get_nt(x,0);
    });
    */
    /*
    merge_dummies(a, o, i, k,
    [&](edge_tag tag, const kmer_t & x, size_t this_k, size_t lcs_len, bool first_end_node) {
      char_writer << get_nt(x, 0);
      #ifdef VAR_ORDER
      out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
      char l(lcs_len);
      lcs.write((char*)&l, 1);
      #else
      out.write(tag, x, this_k, lcs_len, first_end_node);
      #endif
      cols.write((char*)&color, sizeof(color_t));
      prev_k = this_k;
    });
    */
  }
};


} // namespace cosmo

#endif
