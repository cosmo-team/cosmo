#include "common.h"
#include "transform.h"
#include "sort.h"
#include "join.h"
#include "io.h"
#include "uint128_t.h"
#include "nanotime.h"

const char * USAGE = "<DSK output file>";

int main(int argc, char * argv[]) {
  // parse argv file
  if (argc != 2) {
    fprintf(stderr, "Usage: %s %s\n", argv[0], USAGE);
    exit(1);
  }
  char * file_name = argv[1];

  // open file
  int handle = -1;
  if ( (handle = open(file_name, O_RDONLY)) == -1 ) {
    fprintf(stderr, "Can't open file: %s\n", file_name);
    exit(1);
  }

  // read header
  uint32_t kmer_num_bits = 0;
  uint32_t k = 0;
  if ( !dsk_read_header(handle, &kmer_num_bits, &k) ) {
    fprintf(stderr, "Error reading file %s\n", file_name);
    exit(1);
  }
  uint32_t kmer_num_blocks = (kmer_num_bits / 8) / sizeof(uint64_t);
  TRACE("kmer_num_bits, k = %d, %d\n", kmer_num_bits, k);
  TRACE("kmer_num_blocks = %d\n", kmer_num_blocks);

  if (kmer_num_bits > MAX_BITS_PER_KMER) {
    fprintf(stderr, "Kmers larger than %zu bits are not currently supported. %s uses %d bits per kmer.\n",
        MAX_BITS_PER_KMER, file_name, kmer_num_bits);
    exit(1);
  }

  // read how many items there are
  size_t num_records = 0;
  if ( dsk_num_records(handle, kmer_num_bits, &num_records) == -1) {
    fprintf(stderr, "Error seeking file %s\n", file_name);
    exit(1);
  }
  if (num_records == 0) {
    fprintf(stderr, "File %s has no kmers.\n", file_name);
    exit(1);
  }
  TRACE("num_records = %zd\n", num_records);

  // ALLOCATE SPACE FOR KMERS (done in one malloc call)
  typedef uint64_t kmer_t[kmer_num_blocks];
  // x 4 because we need to add reverse complements, and then we have two copies of the table
  kmer_t * kmers = calloc(num_records * 4, sizeof(kmer_t));
  kmer_t * table_a = kmers;
  kmer_t * table_b = kmers + num_records * 2; // x2 because of reverse complements

  // READ KMERS FROM DISK INTO ARRAY
  size_t num_records_read = dsk_read_kmers(handle, kmer_num_bits, (uint64_t*) kmers);
  switch (num_records_read) {
    case -1:
      fprintf(stderr, "Error reading file %s\n", argv[1]);
      exit(1);
      break;
    default:
      TRACE("num_records_read = %zu\n", num_records_read);
      assert (num_records_read == num_records);
      break;
  }

  close(handle);

  // ADD REVERSE COMPLEMENTS
  add_reverse_complements((uint64_t*)kmers, (uint64_t*)(kmers + num_records), num_records, k);

  nanotime_t start, end;
  start = get_nanotime();
  if (kmer_num_bits == 64) {
    //qsort(kmers, num_records*2, sizeof(uint64_t), compare_64);
    // First step so that the edges are sorted after the node is sorted later
    colex_partial_radix_sort_64((uint64_t*)table_a, (uint64_t*)table_b, num_records*2, 1, 0, (uint64_t**)&table_a, (uint64_t**)&table_b);
    // from (k, 0] instead of (k, 1] because we want the colex(row) AND the second last result, which will be colex(node), edge
    // Note that new_a and new_b map to table_b and table_a respectively.
    colex_partial_radix_sort_64((uint64_t*)table_a, (uint64_t*)table_b, num_records*2, k, 0, (uint64_t**)&table_b, (uint64_t**)&table_a);
    // At this point, a will have colex(node), edge sorting (as required for output)
    // and b will have colex(row)
  }
  else if (kmer_num_bits == 128) {
    fprintf(stderr, "NOT YET IMPLEMENTED FOR 128 BIT KMERS\n");
    exit(1);
    //qsort(kmers, num_records*2, sizeof(uint128_t), compare_128);
  }
  end = get_nanotime();
  double ms = (double)(end - start)/(1000000);
  fprintf(stderr, "Sort Time: %.2f ms\n", ms);

  /*
  #ifndef NDEBUG
  printf("TABLE A:\n");
  //print_kmers_hex(stdout, (uint64_t*)table_a, num_records * 2, kmer_num_bits);
  print_kmers_acgt(stdout, (uint64_t*)table_a, num_records * 2, k);
  printf("TABLE B:\n");
  //print_kmers_hex(stdout, (uint64_t*)table_b, num_records * 2, kmer_num_bits);
  print_kmers_acgt(stdout, (uint64_t*)table_b, num_records * 2, k);
  #endif
  */

  // outgoing dummy edges are output in correct order while merging, whereas incoming dummy edges are not in the correct
  // position, but are sorted relatively, hence can be merged if collected in previous passes
  size_t num_incoming_dummies = count_incoming_dummy_edges_64((uint64_t*)table_a, (uint64_t*)table_b, num_records*2, k);
  //TRACE("num_incoming_dummy_edges = %zu\n", num_incoming_dummies);
  uint64_t * incoming_dummies = calloc(num_incoming_dummies*(k-1)*2, sizeof(kmer_t));
  uint64_t * dummies_a = incoming_dummies;
  uint64_t * dummies_b = incoming_dummies + num_incoming_dummies * (k-1);
  unsigned char * incoming_dummy_lengths = calloc(num_incoming_dummies*(k-1)*2, sizeof(kmer_t));
  unsigned char * lengths_a = incoming_dummy_lengths;
  unsigned char * lengths_b = incoming_dummy_lengths + num_incoming_dummies * (k-1);
  get_incoming_dummy_edges_64((uint64_t*)table_a, (uint64_t*)table_b, num_records*2, k, incoming_dummies, num_incoming_dummies);
  prepare_incoming_dummy_edges_64(incoming_dummies, incoming_dummy_lengths, num_incoming_dummies, k-1);
  colex_varlen_partial_radix_sort_64(dummies_a, dummies_b, lengths_a, lengths_b, num_incoming_dummies*(k-1), 1, 0, &dummies_a, &dummies_b, &lengths_a, &lengths_b);
  colex_varlen_partial_radix_sort_64(dummies_a, dummies_b, lengths_a, lengths_b, num_incoming_dummies*(k-1), k-1, 1, &dummies_a, &dummies_b, &lengths_a, &lengths_b);
/*
#ifndef NDEBUG
  printf("Incoming Dummies:\n");
  print_dummies_acgt(stdout, dummies_a, lengths_a, num_incoming_dummies*(k-1), k);
#endif
*/
  // output in ascii first
  //merge_and_output(stderr, table_a, table_b, incoming_dummies, num_records*2, k);
  fprintf(stderr, "MERGING DUMMIES\n");
  merge_dummies(stderr, (uint64_t*)table_a, (uint64_t*)table_b, num_records*2, k, dummies_a, num_incoming_dummies*(k-1), lengths_a);
  free(incoming_dummies);
  free(incoming_dummy_lengths);
  free(kmers);

  return 0;
}
