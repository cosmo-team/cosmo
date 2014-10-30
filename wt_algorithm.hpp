#pragma once
#ifndef _WT_ALGORITHM_H
#define _WT_ALGORITHM_H

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
  if (i == 0) return 0;

  // at a leaf, we return the current position
  if (wt.is_leaf(node)) {
    return i;
  }

  auto children = wt.expand(node);
  auto left_child = std::get<0>(children);
  auto right_child = std::get<1>(children);

  uint64_t mask = (1ULL) << (wt.max_level - 1);
  if (!(c & (mask >> node.level))) { // left
    size_t right_rank = wt.node_rank1(node, i);
    size_t right_pos  = (right_rank)? wt.node_select1(node, right_rank + 1) + 1: 0;
    size_t left_rank      = wt.node_rank0(node, i);
    size_t left_sub_rank  = (left_rank)? _prev_lte_rec(wt, left_child, left_rank, c) : 0;
    size_t left_pos       = (left_sub_rank)? wt.node_select0(node, left_sub_rank+1) + 1 : 0;
    return std::min(left_pos, right_pos);
  } else { // right
    size_t rank = wt.node_rank1(node, i);
    size_t sub  = (rank)? _prev_lte_rec(wt, right_child, rank, c) : 0;
    size_t pos  = (sub)? wt.node_select1(node, sub) + 1 : 0;
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

/*
template <class t_wt>
size_t next_lte(const t_wt & wt, size_t i, typename t_wt::value_type c) {
  static_assert(t_wt::lex_ordered, "next_lte requires a lex_ordered WT");
  return 0;
}
*/

/*
template <class t_wt>
size_t range_lte(const t_wt & wt, size_t i, size_t j, typename t_wt::value_type c) {
  return 0;
}
*/

#endif
