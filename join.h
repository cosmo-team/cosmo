#ifndef JOIN_H
#define JOIN_H

#include "common.h"

size_t count_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k);
//size_t count_outgoing_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k);

#endif
