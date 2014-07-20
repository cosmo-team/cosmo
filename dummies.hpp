#ifndef DUMMIES_HPP
#define DUMMIES_HPP

// reimplement this with lambdas
size_t count_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k);
void get_incoming_dummy_edges_64(uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k, uint64_t * incoming_dummies, size_t num_incoming_dummies);
void prepare_incoming_dummy_edges_64(uint64_t * dummy_nodes, unsigned char * k_values, size_t num_dummies, uint32_t k);
void merge_dummies(FILE * outfile, uint64_t * table_a, uint64_t * table_b, size_t num_records, uint32_t k, uint64_t * incoming_dummies, size_t num_incoming_dummies, unsigned char * dummy_lengths);

#endif
