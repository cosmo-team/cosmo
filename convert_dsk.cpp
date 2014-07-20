// STL Headers
//#include <fstream>
#include <iostream>
#include <algorithm>

// C STDLIB Headers
#include <cstdio>
#include <cstdlib>

// Custom Headers
#include "uint128_t.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "debug.h"
#include "nanotime.h"

using namespace std;

static const char * USAGE = "<DSK output file>";

template <typename kmer_t>
void convert_representation(kmer_t * kmers_in, kmer_t * kmers_out, size_t num_kmers) {
  // Swap G and T and reverses the nucleotides so that they can
  // be compared and sorted as integers to give colexicographical ordering
  transform((uint64_t*)kmers_in, (uint64_t*)(kmers_in + num_kmers), (uint64_t*)kmers_out, swap_gt);
  transform(kmers_in, kmers_out + num_kmers, kmers_in, reverse_nt<kmer_t>());
}

template <typename kmer_t>
void print_kmers(ostream & out, kmer_t * kmers, size_t num_kmers, uint32_t k, uint8_t * dummy_lengths = 0) {
  for (size_t i = 0; i<num_kmers; i++) {
    out << kmer_to_string(kmers[i], k, ((dummy_lengths)? dummy_lengths[i]:k)) << endl;
  }
}

int main(int argc, char * argv[]) {
  // Parse argv
  if (argc != 2) {
    fprintf(stderr, "Usage: %s\n", USAGE);
    exit(EXIT_FAILURE);
  }
  const char * file_name = argv[1];

  // Open File
  int handle = -1;
  if ( (handle = open(file_name, O_RDONLY)) == -1 ) {
    fprintf(stderr, "ERROR: Can't open file: %s\n", file_name);
    exit(EXIT_FAILURE);
  }

  // Read Header
  uint32_t kmer_num_bits = 0;
  uint32_t k = 0;
  if ( !dsk_read_header(handle, &kmer_num_bits, &k) ) {
    fprintf(stderr, "ERROR: Error reading file %s\n", file_name);
    exit(EXIT_FAILURE);
  }
  uint32_t kmer_num_blocks = (kmer_num_bits / 8) / sizeof(uint64_t);
  TRACE(">> READING DSK FILE\n");
  TRACE("kmer_num_bits, k = %d, %d\n", kmer_num_bits, k);
  TRACE("kmer_num_blocks = %d\n", kmer_num_blocks);

  if (kmer_num_bits > MAX_BITS_PER_KMER) {
    fprintf(stderr, "ERROR: Kmers larger than %zu bits are not currently supported."
                    " %s uses %d bits per kmer (possibly corrupt?).\n",
        MAX_BITS_PER_KMER, file_name, kmer_num_bits);
    exit(EXIT_FAILURE);
  }

  // Read how many items there are (for allocation purposes)
  size_t num_kmers = 0;
  if ( dsk_num_records(handle, kmer_num_bits, &num_kmers) == -1) {
    fprintf(stderr, "Error seeking file %s\n", file_name);
    exit(EXIT_FAILURE);
  }
  if (num_kmers == 0) {
    fprintf(stderr, "ERROR: File %s has no kmers (possibly corrupt?).\n", file_name);
    exit(EXIT_FAILURE);
  }
  TRACE("num_kmers = %zd\n", num_kmers);

  // ALLOCATE SPACE FOR KMERS (done in one malloc call)
  // x 4 because we need to add reverse complements, and then we have two copies of the table
  uint64_t * kmer_blocks = (uint64_t*)calloc(num_kmers * 4, sizeof(uint64_t) * kmer_num_blocks);
  //kmer_t * table_a = kmers;
  //kmer_t * table_b = kmers + num_records * 2; // x2 because of reverse complements

  // READ KMERS FROM DISK INTO ARRAY
  size_t num_records_read = dsk_read_kmers(handle, kmer_num_bits, kmer_blocks);
  close(handle);
  switch (num_records_read) {
    // Case 0 would have been captured above;
    case 0:
      fprintf(stderr, "Error reading file %s\n", argv[1]);
      exit(EXIT_FAILURE);
      break;
    default:
      TRACE("num_records_read = %zu\n", num_records_read);
      assert (num_records_read == num_kmers);
      break;
  }

  // TODO: move this to a template or somethin...
  if (kmer_num_bits == 64) {
    typedef uint64_t kmer_t;
    kmer_t * kmers = (kmer_t*)kmer_blocks;
    convert_representation(kmers, kmers, num_kmers);
    transform(kmers, kmers + num_kmers, kmers + num_kmers, reverse_complement<kmer_t>(k));
    print_kmers(cout, kmers, num_kmers * 2, k);
  }
  else if (kmer_num_bits == 128) {
    typedef uint128_t kmer_t;
    kmer_t * kmers = (kmer_t*)kmer_blocks;
    convert_representation(kmers, kmers, num_kmers);
    transform(kmers, kmers + num_kmers, kmers + num_kmers, reverse_complement<kmer_t>(k));
    print_kmers(cout, kmers, num_kmers * 2, k);
  }

  // add reverse_complements
  // sort
  // allocate space for dummies
  // extract dummies
  // sort dummies
  // merge
  // output (iterator?) -> write to file
  // free dummies

  free(kmer_blocks);
  return 0;
}
