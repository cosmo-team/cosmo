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
    // At this point, a will have colex(node), edge sorting (as required)
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

  #ifndef NDEBUG
  printf("TABLE A:\n");
  //print_kmers_hex(stdout, (uint64_t*)table_a, num_records * 2, kmer_num_bits);
  print_kmers_acgt(stdout, (uint64_t*)table_a, num_records * 2, k);
  printf("TABLE B:\n");
  //print_kmers_hex(stdout, (uint64_t*)table_b, num_records * 2, kmer_num_bits);
  print_kmers_acgt(stdout, (uint64_t*)table_b, num_records * 2, k);
  #endif

  size_t num_incoming_dummy_edges = count_incoming_dummy_edges((uint64_t*)table_a, (uint64_t*)table_b, k);
  //size_t num_outgoing_dummy_edges = count_outgoing_dummy_edges((uint64_t*)table_a, (uint64_t*)table_b, k);

  // TODO: implement joining algorithm for 64 bit kmers (counting first)
  // join.h, join.c
  // count_incoming_dummy_edges
  // count_outgoing_dummy_edges
  // allocate space for in and out dummies
  // find_incoming_dummy_edges
  // find_outgoing_dummy_edges
  // dummy: position integer and kmers
  // TODO: Support 128-bit kmers (for sorting and printing)
  // TODO: implement joining algorithm for 128 bit kmers
  // TODO: count how many dummy edges there are and allocate that much space
  // find_dummy_edges(table_a, table_b, k, dummy_out_ptr, dummy_in_ptr)
  // TODO: count how many of each symbol there are in the 2nd last column
  // TODO: filter leftmost k-1 symbols to add LAST flag
  // TODO: filter for the rightmost k-1 symbols to add minus (have a seen boolean for each symbol)
  // TODO: implement buffered output for our fmt:
  // Header: k, num records (incl dummy, etc), ACGT table
  // uint64_t vector of 4x5 bits = 4 edges.
  // make struct for 5 bits: last_flag(1), symbol(2), minus_flag(1), dummy_flag(1)
  // make struct for record (4 bits waste + array for 4 records) - check size

  free(kmers);
  return 0;
}
