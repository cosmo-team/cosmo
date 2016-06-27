#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <cassert>

//#include <boost/log/trivial.hpp>
#include <boost/assert.hpp>

#define BOOST_LOG_DYN_LINK 1
// Aliases (in case we want to change the behaviour later)
#define COSMO_ASSERT BOOST_ASSERT
#define COSMO_LOG(x) std::cerr
//BOOST_LOG_TRIVIAL
#define TRACE COSMO_LOG(trace)
// log levels: trace, debug, info, warning, error, fatal

#endif
