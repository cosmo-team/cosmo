#ifndef DEBUG_H
#define DEBUG_H

#include <assert.h>

// Default is to be opposite NDEBUG
#define DEBUG !NDEBUG

#if DEBUG
    #include <stdio.h>
    #define TRACE(...) {fprintf(stderr, __VA_ARGS__);}
#else
    #define TRACE(...)
#endif

#endif /* end of include guard */
