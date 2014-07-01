#ifndef SORT_H
#define SORT_H

#include "common.h"
#include "uint128_t.h"

void colex_partial_radix_sort_64(uint64_t * a, uint64_t * b, size_t num_records, uint32_t k, uint32_t j, uint64_t ** new_a, uint64_t** new_b);
// comparators for quicksort
//int compare_64(const void *, const void *);
//int compare_128(const void *, const void *);
void qsort_varlen(uint64_t * array, unsigned char * lengths, int lower, int upper);
//void radix_lsd_64();
//void radix_lsd_64();

#endif
