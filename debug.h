#ifndef DEBUG_H
#define DEBUG_H

#include <assert.h>

// Default is to be opposite NDEBUG, but debug levels may be added sometime...
// Other ideas: Specifying multiple logging sources (with C++)
#ifndef NDEBUG
#define DEBUG
#endif

#ifdef DEBUG
    #include <stdio.h>
    #define TRACE(...) {fprintf(stderr, __VA_ARGS__);}
#else
    #define TRACE(...)
#endif

#endif /* end of include guard */
