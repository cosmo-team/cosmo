#pragma once
#ifndef DUMMIES_HPP
#define DUMMIES_HPP

#include <boost/heap/priority_queue.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/variant.hpp>
#include <parallel/algorithm>
#include <boost/range/adaptor/transformed.hpp>     // Map function to inputs
#include <boost/range/adaptors.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/uniqued.hpp>         // Uniquify
#include <boost/range/algorithm/set_algorithm.hpp> // set_difference
#include <boost/function_output_iterator.hpp>      // for capturing output of set_algorithms
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <algorithm>
#include <utility>                                 // make_pair (for boost ranges)
#include <functional>                              // function (to avoid errors with the lambdas)
#include <cstring>                                 // memset

#include "debug.hpp"
#include "config.hpp"
#include "kmer.hpp"

using namespace boost::adaptors;
using boost::get;

enum edge_tag { out_dummy, in_dummy, standard};

template <typename kmer_t>
struct uniq {
  kmer_t prev;
  bool   first = true;

  bool operator()( kmer_t x ) {
    if (x == prev && !first) {
      return false;
    }
    else {
      prev = x;
      first = false;
      return true;
    }
  }
};

template <typename kmer_t, typename InputRange1, typename InputRange2, typename Func>
void find_incoming_dummy_nodes(const InputRange1 a_range, const InputRange2 b_range, uint32_t k, Func out_f) {
  //typedef decltype(*a_range.begin()) kmer_t;
  //typedef typename OutputIterator::value_type pair_t;
  size_t idx = 0;
  kmer_t temp;
  auto a_lam = std::function<kmer_t(kmer_t)>([&](kmer_t x) -> kmer_t {
    temp = x;
    return get_start_node(x);
  });
  auto b_lam   = std::function<kmer_t(kmer_t)>([k](kmer_t x) -> kmer_t {return get_end_node(x,k);});
  auto a = a_range | transformed(a_lam) | filtered(uniq<kmer_t>());
  auto b = b_range | transformed(b_lam) | filtered(uniq<kmer_t>());

  auto pairer  = [&](kmer_t) { out_f(idx++, temp); };
  auto paired_out = boost::make_function_output_iterator(pairer);
  //std::set_difference(a_start, a.end(), b.begin(), b.end(), paired_out);
  boost::set_difference(a, b, paired_out);
  // GPU: http://thrust.github.io/doc/group__set__operations.html
  //__gnu_parallel::set_difference(a.begin(), a.end(), b.begin(), b.end(), paired_out);
}

// TODO: use set_symmetric_difference with each wrapped with a tag that is ignored during the comparison
// that way we can find both dummy types at same time
template <typename kmer_t, typename InputRange1, typename InputRange2, typename Func>
void find_outgoing_dummy_nodes(const InputRange1 a_range, const InputRange2 b_range, uint32_t k, Func out_f) {
  kmer_t temp;
  auto a_lam = std::function<kmer_t(kmer_t)>([](kmer_t x) -> kmer_t {
    return get_start_node(x);
  });
  auto b_lam = std::function<kmer_t(kmer_t)>([&](kmer_t x) -> kmer_t {return get_end_node((temp = x),k);});
  auto a = a_range | transformed(a_lam) | filtered(uniq<kmer_t>());
  auto b = b_range | transformed(b_lam) | filtered(uniq<kmer_t>());

  auto pairer = [&](kmer_t) { out_f(temp); };
  auto paired_out = boost::make_function_output_iterator(pairer);
  boost::set_difference(b, a, paired_out);
  // GPU: http://thrust.github.io/doc/group__set__operations.html
  //__gnu_parallel::set_difference(a.begin(), a.end(), b.begin(), b.end(), paired_out);
}

template <class Visitor>
class Unique {
  Visitor _v;
  bool first_iter = true;
  edge_tag last_tag;
  uint32_t last_k;

  public:
    Unique(Visitor v) : _v(v) {}

    template <typename kmer_t>
    void operator()(edge_tag tag, const kmer_t & x, uint32_t k) {
      static kmer_t last_kmer;

      if (first_iter || tag != last_tag || x != last_kmer || k != last_k) {
        _v(tag, x, k);
      }
      first_iter = false;
      last_tag   = tag;
      last_kmer  = x;
      last_k     = k;
    }
};

class FirstStartNodeFlagger {
  uint32_t _graph_k;

  public:
    FirstStartNodeFlagger(uint32_t k) : _graph_k(k) {}
    template <typename kmer_t>
    bool operator()(const kmer_t & x, const uint32_t k) {
      // Note: for multithreading, this might be dangerous? could make class members instead
      static uint32_t last_k = 0;

#ifdef VAR_ORDER
      static kmer_t last_edge = 0;
      size_t result = node_lcs(x, last_edge, std::min(k, last_k));
      // If dummy edge length is equal, and the LCS length is the length of the non-dummy suffix,
      // include the $ signs as well
      if (last_k == k && k - 1 == result)
        result = _graph_k - 1;
      last_edge = x;
      // This should bump to K when it is == this_k - 1
      //std::cerr << "kmer: " << kmer_to_string(x,k) << std::endl;
      //std::cerr << "lcs: " << result << std::endl;
#else
      static kmer_t last_start_node = 0;
      static bool first_iter = true;
      kmer_t this_start_node = get_start_node(x);
      size_t result = (first_iter || this_start_node != last_start_node || k != last_k);
      first_iter = false;
      last_start_node  = this_start_node;
#endif

      last_k     = k;

      return result;
    }
};

class FirstEndNodeFlagger {
  uint32_t _graph_k;
  bool first_iter = true;

  public:
    FirstEndNodeFlagger(uint32_t k) : _graph_k(k) {}
    template <typename kmer_t>
    bool operator()(edge_tag tag, const kmer_t & x, const uint32_t k) {
      static kmer_t last_suffix;
      static uint32_t last_k;
      static bool edge_seen[DNA_RADIX];
      #define reset_flags() memset(edge_seen, 0, DNA_RADIX)

      bool edge_flag = true;
      kmer_t this_suffix = get_start_node_suffix(x, _graph_k);

      // reset "edge seen" flags
      if (this_suffix != last_suffix || k != last_k || first_iter) {
        reset_flags();
        edge_flag = true;
        first_iter = false;
      }
      // Only get the edge label if not a dummy out edge
      if (tag != out_dummy){
        uint8_t edge = get_edge_label(x);
        edge_flag = !edge_seen[edge];
        edge_seen[edge] = true;
      }

      last_suffix = this_suffix;
      last_k      = k;

      return edge_flag;
    }
};

template <class Visitor>
auto uniquify(Visitor v) -> Unique<decltype(v)> {
  return Unique<decltype(v)>(v);
}

struct incrementer : boost::static_visitor<> {
  template <typename T>
  void operator()(T & t) const {
    t->advance_begin(1);
  }
};

struct is_empty2 : boost::static_visitor<bool> {
  template <typename T>
  bool operator()(T & t) const {
    return t->empty();
  }
};

template <typename R>
struct fronter : boost::static_visitor<R> {
  template <typename T>
  R operator()(T t) const {
    return (R)*(t->begin());
  }
};

template <typename value_type, typename first_range, typename ... ranges>
struct heap_item {
  typedef heap_item<value_type, first_range, ranges...> this_type;
  typedef boost::variant<first_range, ranges...> range_variant;
  //typedef first_range:: value_type;
  range_variant m_range;

  typename boost::heap::fibonacci_heap<this_type>::handle_type handle;

  //heap_item(range_variant r) : m_range(&r) {}
  template <typename range>
  heap_item(range r) : m_range(r) {}
  //template <typename range>

  // new incremented heap item
  this_type & operator++() {
    boost::apply_visitor(incrementer(), m_range);
    return *this;
  }

  bool empty() const {
    return boost::apply_visitor(is_empty2(), m_range);
  }

  value_type front() const {
    const auto f = fronter<value_type>{};
    return boost::apply_visitor(f, m_range);
  }

  friend bool operator<(const this_type & l, const this_type & r) {
    // TODO : fix the comparison here (front() will return an element_t
    // which is a tuple with record_t, this_k, dummy_tag
    auto l_element = l.front();
    auto l_record  = get<0>(l_element);
    auto l_edge    = get<0>(l_record);
    auto l_node    = get_start_node(l_edge);
    auto l_k       = get<1>(l_element);
    //auto l_tag     = get<2>(l_element);

    auto r_element = r.front();
    auto r_record  = get<0>(r_element);
    auto r_edge    = get<0>(r_record);
    auto r_node    = get_start_node(r_edge);
    auto r_k       = get<1>(r_element);
    //auto r_tag     = get<2>(r_element);

    // compare nodes. If nodes are equal, compare length.
    // if length is equal, return edge < edge
    // else return length < length
    // else return node < node
    // swapped because we want min heap, but a heap defaults to max
    if (r_node == l_node) {
      if (r_k == l_k) {
        return r_edge < l_edge;
      }
      else return (r_k < l_k);
    }
    else return (r_node < l_node);
  }
};

template <typename InputRange1, typename InputRange2, typename InputRange3, class Visitor>
void merge_dummies(const InputRange1 & a_range, const InputRange2 & o_range, const InputRange3 & i_range, uint8_t k, Visitor visit) {
  typedef typename InputRange1::value_type record_t;
  typedef typename InputRange2::value_type kmer_t;

  size_t idx  = 0;
  auto i_next = i_range.cbegin();
  auto a = a_range | transformed([&](record_t x){
    bool is_dummy = (i_next != i_range.cend()) && (idx == *(i_next));
    return boost::make_tuple(x, is_dummy?in_dummy:standard);
  });

  auto o = o_range | transformed([](kmer_t x) {
    return boost::make_tuple(record_t(x), out_dummy);
  });

  using boost::get;
  auto v   = [&](auto x) {
    if (get<1>(x) != out_dummy) {
      ++idx;
      if (get<1>(x) == in_dummy) ++i_next;
    }
    visit(boost::get<1>(x), boost::get<0>(x), k);
  };
  auto out = boost::make_function_output_iterator(v);
  boost::set_union(a, o, out, [](boost::tuple<record_t, edge_tag> x,
                                 boost::tuple<record_t, edge_tag> y){
    return (get<0>(get<0>(x))<<2) < (get<0>(get<0>(y))<<2);
  });
}

// TODO: Try a parallell merge (and set difference) on the internal STXXL blocks
template <typename InputRange1, typename InputRange2, typename InputRange3, class Visitor>
void merge_dummies_with_shifts(const InputRange1 & a_range, const InputRange2 & o_range, const InputRange3 & i_range, const uint8_t k, Visitor visit) {
  using namespace boost::heap;
  typedef typename InputRange1::value_type record_t;
  typedef boost::tuple<record_t, uint8_t, edge_tag> element_t;

  auto a = a_range | transformed([k](auto x) {
    return boost::make_tuple(x, k, standard);
  });

  auto o = o_range | transformed([k](auto x) {
    return boost::make_tuple(record_t(x), k, out_dummy);
  });

  auto i = i_range | transformed([](auto x) {
    return boost::make_tuple(record_t(get<0>(x)), get<1>(x), in_dummy);
  });

  // Make sure queue can hold references to each variant of range
  typedef heap_item<element_t, decltype(a)*, decltype(o)*, decltype(i)*> h_t;

  // Make min-heap priority queue with elements being the lists
  // uses boost::heap::fibonacci_heap for O(1) update of iterators
  // (although it uses nodes, so the indirection might be slower for a small number of input files)
  //fibonacci_heap<h_t> q;
  //TODO: work out if using a fib heap is possible?
  //typedef typename fibonacci_heap<h_t>::handle_type handle_t;
  boost::heap::priority_queue<h_t> q;

  // Add each range to queue
  if(!a.empty()) {
    q.push(h_t(&a));
    //handle_t h = q.push(h_t(&a));
    //(*h).handle = h;
  }
  if(!o.empty()) {
    q.push(h_t(&o));
    //handle_t h = q.push(h_t(&o));
    //(*h).handle = h;
  }
  if(!i.empty()) {
    q.push(h_t(&i));
    //handle_t h = q.push(h_t(&i));
    //(*h).handle = h;
  }

  while(!q.empty()) {
    h_t input = q.top();
    q.pop();
    auto element = input.front();
    auto record = get<0>(element);
    auto this_k = get<1>(element);
    auto tag    = get<2>(element);
    visit(tag, record, this_k);
    ++input;
    // Fibonacci heaps let us update with O(1), so instead of popping
    // and repushing, we only pop if its empty
    // and signal that the next element in that node will have increased (in sorted order)
    //if (input.empty()) q.pop();
    //else q.increase(input.handle);
    if (!input.empty()) q.push(input);
  }
}

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
