#pragma once
#ifndef _WT_ALGORITHM_H
#define _WT_ALGORITHM_H

#include <vector>
#include <string>
#include <sdsl/bits.hpp>
#include <boost/optional.hpp>
#include <sdsl/wt_algorithm.hpp>

template <class t_wt>
size_t _prev_lte_rec(const t_wt & wt, const typename t_wt::node_type & node, size_t i, typename t_wt::value_type c) {
  if (i == 0) return 0;

  // at a leaf, we return the current position
  if (wt.is_leaf(node)) {
    return i;
  }

  auto children = wt.expand(node);
  auto left_child = std::get<0>(children);
  auto right_child = std::get<1>(children);

  uint64_t mask = (1ULL) << (wt.max_level - 1);
  if (c & (mask >> node.level)) {
    // Left tree elements are guaranteed to be < c, so no need to traverse further down that side
    size_t left_rank = wt.node_rank0(node, i);
    // select to count to the position of most recent element
    size_t left_pos  = (left_rank)? wt.node_select0(node, left_rank) + 1: 0;
    // unless left_pos == i, we have larger elements in our range that *could* be < c
    // we have to traverse the right subtree to find their values
    size_t right_rank     = wt.node_rank1(node, i);
    size_t right_sub_rank = (right_rank)? _prev_lte_rec(wt, right_child, right_rank, c) : 0;
    size_t right_pos  = (right_sub_rank)? wt.node_select1(node, right_sub_rank) + 1 : 0;
    return std::max(left_pos, right_pos);
  } else {
    size_t rank = wt.node_rank0(node, i);
    size_t sub  = (rank)? _prev_lte_rec(wt, left_child, rank, c) : 0;
    size_t pos  = (sub)? wt.node_select0(node, sub) + 1 : 0;
    return pos;
  }
}

template <class t_wt>
size_t _next_lte_rec(const t_wt & wt, const typename t_wt::node_type & node, size_t i, typename t_wt::value_type c) {
  string pad(node.level * 2, ' ');
  pad += "  ";

  if (i > node.size) {
    return node.size+1;
  }

  // at a leaf, we return the current position
  if (wt.is_leaf(node)) {
    return i;
  }
  // Basic idea is the same, but need to select to next element.
  // and check if it goes past the global rank of that symbol

  auto children = wt.expand(node);
  auto left_child = std::get<0>(children);
  auto right_child = std::get<1>(children);

  uint64_t mask = (1ULL) << (wt.max_level - 1);
  if (c & (mask >> node.level)) {
    size_t node_bit = (wt.node_access(node, i-1)==1);
    // right branch
    // find the first 0 that occurs at i' >= i
    size_t r0 = wt.node_rank0(node, i); // # 0s in [0,i)
    size_t extend0 = (node_bit == 1);
    size_t p = wt.node_select0(node, r0+extend0) + 1;
    if (p == i) return p;
    size_t r1 = wt.node_rank1(node, i); // # 1s in [0, i)
    size_t extend1 = (node_bit==0);// || node_bit==0);// || (r1 < i && node_bit==0));
    size_t rr = _next_lte_rec(wt, right_child, r1 + extend1, c);
    size_t pos = wt.node_select1(node, rr) + 1;
    return std::min(std::min(p, pos), (size_t)node.size+1);
  } else {
    // bit = 0 - must recurse on left branch
    // number of 0s in [0, i)
    size_t r0  = wt.node_rank0(node, i);
    //size_t extend = (lr == 0 || (r0 < i && wt.node_access(node, i-1)==1));
    size_t extend = (wt.node_access(node, i-1)==1);
    size_t lr  = _next_lte_rec(wt, left_child, r0+extend, c);
    size_t pos = wt.node_select0(node, lr) + 1;
    if (pos > node.size) return node.size+1;
    return pos;
  }
}

// Answers on [0..i) range (so 0 returns 0, much like sdsl rank)
template <class t_wt>
size_t prev_lte(const t_wt & wt, size_t i, typename t_wt::value_type c) {
  static_assert(t_wt::lex_ordered, "prev_lte requires a lex_ordered WT");

  if (i == 0) return 0;

  // Probs should do this better (w.r.t. cache)
  auto x = sdsl::symbol_lte(wt, c);
  if (!get<0>(x)) return 0;
  c = get<1>(x);

  return _prev_lte_rec(wt, wt.root(), i, c);
}

// Answers on [0..i) range (so 0 returns 0, much like sdsl rank)
template <class t_wt>
size_t next_lte(const t_wt & wt, size_t i, typename t_wt::value_type c) {
  static_assert(t_wt::lex_ordered, "next_lte requires a lex_ordered WT");

  if (i == 0) return 0;

  auto x = sdsl::symbol_lte(wt, c);
  if (!get<0>(x)) return 0;
  c = get<1>(x);

  return _next_lte_rec(wt, wt.root(), i, c);
}

// make iterators instead
// random access next_lte (select_lte) is probably quite possible as well
template <class t_wt>
vector<size_t> range_lte(const t_wt & wt, size_t i, size_t j, typename t_wt::value_type c) {
  vector<size_t> range;
  for (size_t next = i; next <= j; next++) {
    next = next_lte(wt,next,c);
    if (next > j) break;
    range.push_back(next);
  }
  return range;
}

#endif
