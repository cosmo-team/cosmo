#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <cassert>

#include <boost/log/trivial.hpp>
#include <boost/assert.hpp>

#define BOOST_LOG_DYN_LINK 1
// Aliases (in case we want to change the behaviour later)
#define COSMO_LOG BOOST_LOG_TRIVIAL
#define COSMO_ASSERT BOOST_ASSERT

// log levels: trace, debug, info, warning, error, fatal

#endif
