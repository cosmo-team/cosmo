#pragma once
#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <climits> // For CHAR_BIT

template <typename T>
size_t bitwidth() {
  return sizeof(T) * CHAR_BIT;
}

template <typename T>
size_t bitwidth(const T&) {
  return sizeof(T) * CHAR_BIT;
}

#endif
