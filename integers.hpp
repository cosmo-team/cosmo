#pragma once
#ifndef INTEGERS_HPP
#define INTEGERS_HPP

#include <ttmath/ttmath.h> // http://www.ttmath.org/

namespace cosmo {

// Have to define class so I can overload the istream >> operator
// TODO: Add check for 32 bit OS
/*
template <size_t num_blocks>
struct big_uint : ttmath::UInt<num_blocks> {
};

// Input operator used for reading from the file
template <size_t num_blocks>
std::istream& operator >> (std::istream& in, big_uint<num_blocks>& x) {
  in.read((char*)&x, sizeof(big_uint<num_blocks>));
  return in;
}
*/

typedef ttmath::UInt<2> uint128_t;

}

#endif
