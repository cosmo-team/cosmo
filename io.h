#ifndef IO_H
#define IO_H

#include "common.h"

static const size_t MAX_BITS_PER_KMER = 128;
static const size_t BUFFER_SIZE = 0x8000; // 32Kb data buffer

#define DSK_FILE_RECORD_SIZE(NBITS) (((NBITS)/8) + 4)

// Reads the DSK file input header
int dsk_read_header(int, uint32_t *, uint32_t *);
// Counts the number of records in the file - for allocation purposes
int dsk_num_records(int handle, uint32_t kmer_num_bits, size_t * num_records);
// Read kmers from file into the output array
size_t dsk_read_kmers(int handle, uint32_t kmer_num_bits, uint64_t * kmers_output);
void merge_and_output(FILE * outfile, uint64_t * table_a, uint64_t * table_b, uint64_t * incoming_dummies, size_t num_records, size_t num_incoming_dummies, uint32_t k);

#endif
