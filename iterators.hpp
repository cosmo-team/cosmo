#pragma once
#ifndef ITERATORS_HPP
#define ITERATORS_HPP

#include <algorithm>
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>
#include "debug.h"

template <class InputIterator1, class InputIterator2>
class set_difference_iterator 
  : public boost::iterator_facade<set_difference_iterator<InputIterator1, InputIterator2>,
                                  // InputIterator1 and InputIterator2 have the same value_type
                                  typename std::iterator_traits<InputIterator1>::value_type,
                                  std::input_iterator_tag> {
  public:
  typedef typename std::iterator_traits<InputIterator1>::value_type value_type;

  set_difference_iterator() : _first1(), _last1(), _first2(), _last2(), _next() {}
  set_difference_iterator(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2)
    : _first1(first1), _last1(last1), _first2(first2), _last2(last2), _next() { increment(); }
  // Note that we have to increment the iterator ONCE to point to the correct next elemenet (for comparisons)

  private:
  // This is to separate incrementing and dereferencing code. The Iterators have the same value type.

  InputIterator1 _first1;
  InputIterator1 _last1;
  InputIterator2 _first2;
  InputIterator2 _last2;

  InputIterator1 _next;

  friend class boost::iterator_core_access;

  void increment() {
    // Assumes that we call the increment() func once in the ctor
    // so first1 and first2 are really saying where to point to NEXT
    while (_first1 != _last1 && _first2 != _last2)
    {
      if (*_first1 < *_first2) {
        _next = _first1;
        ++_first1;
        return;
      }
      else if (*_first2 < *_first1) {
        ++_first2;
      }
      else {
        ++_first1;
        ++_first2;
      }
    }
    _next = _first1;
    ++_first1;
  }

  value_type& dereference() const {
    return *_next;
  }

  bool equal(set_difference_iterator<InputIterator1, InputIterator2> const& other) const {
    return _next == other._next;
  }
};

template <class InputIterator1, class InputIterator2,
          class Compare = std::less<typename std::iterator_traits<InputIterator1>::value_type> >
class merge_iterator
  : public boost::iterator_facade<merge_iterator<InputIterator1, InputIterator2>,
                                  // InputIterator1 and InputIterator2 have the same value_type
                                  typename std::iterator_traits<InputIterator1>::value_type,
                                  std::input_iterator_tag> {
  public:
  typedef typename std::iterator_traits<InputIterator1>::value_type value_type;

  // constructors
  merge_iterator() : _first1(), _last1(), _first2(), _last2(), _next(), _less() {}
  merge_iterator(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2, Compare less = Compare())
    : _first1(first1), _last1(last1), _first2(first2), _last2(last2), _less(less), _next() { increment(); }
  // Note that we have to increment the iterator ONCE to point to the correct next elemenet (for comparisons)

  private:
  InputIterator1 _first1;
  InputIterator1 _last1;
  InputIterator2 _first2;
  InputIterator2 _last2;

  value_type * _next;

  Compare _less;

  friend class boost::iterator_core_access;

  void increment() {
    if (_first1 == _last1 && _first2 == _last2) _next = 0;
    else if (_first1 == _last1) {_next = &(*_first2); ++_first2;}
    else if (_first2 == _last2) {_next = &(*_first1); ++_first1;}
    else if (_less(*_first1, *_first2)) {_next = &(*_first1); ++_first1;}
    else {_next = &(*_first2); ++_first2;}
  }

  value_type& dereference() const {
    return *_next;
  }

  bool equal(merge_iterator<InputIterator1, InputIterator2, Compare> const& other) const {
    return _next == other._next;
  }
};

#endif
