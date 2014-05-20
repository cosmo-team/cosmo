#include "common.h"
#include "transform.h"
#include "sort.h"
#include "io.h"
#include "uint128_t.h"
#include "nanotime.h"

const char * USAGE = "<DSK output file>";

int compare_64(const void *, const void *);
int compare_128(const void *, const void *);

int compare_64(const void * lhs, const void * rhs)
{
  uint64_t a = *(uint64_t*)lhs;
  uint64_t b = *(uint64_t*)rhs;
  if (a < b) {
    return -1;
  }
  else if (a > b){
    return 1;
  }
  return 0;
}

int compare_128(const void * lhs, const void * rhs) {
  uint128_t a = *(uint128_t*)lhs;
  uint128_t b = *(uint128_t*)rhs;
  if (a.upper < b.upper)      return -1;
  else if (a.upper > b.upper) return  1;
  else if (a.lower < b.lower) return -1;
  else if (a.lower > b.lower) return  1;
  return 0;
}

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
  kmer_t * kmers = calloc(num_records * 2, sizeof(kmer_t));
  kmer_t * table_a = kmers;
  kmer_t * table_b = kmers + num_records * 2; // x2 because of reverse complements

  // READ KMERS FROM DISK INTO ARRAY
  size_t num_records_read = num_records_read = dsk_read_kmers(handle, kmer_num_bits, (uint64_t*) kmers);
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

  // TODO: Manipulate bits so we can sort in colex(node) order
  // Which is to say, just rotate RIGHT one place
  // Then the second sort can be a stable sort on the LSD
  // x & 0x3 << (k - 1) * 2 | x >> 2

  #ifndef NDEBUG
  print_kmers_hex(stdout, (uint64_t*)kmers, num_records * 2, kmer_num_bits);
  #endif

  nanotime_t start, end;
  start = get_nanotime();
  // TODO: Change sorting function to be LSD radix
  // Can use the second table as a temp one
  if (kmer_num_bits == 64) {
    //qsort(kmers, num_records*2, sizeof(uint64_t), compare_64);
    lsd_radix_sort_64((uint64_t*)table_a, (uint64_t*)table_b, num_records*2, k);
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
  printf("SORTING\n");
  print_kmers_hex(stdout, (uint64_t*)kmers, num_records * 2, kmer_num_bits);
  #endif

  // SECOND SORT PHASE
  // Store a copy of the kmers to sort in colex(row) order, for joining
  // in order to calculate dummy edges
  //memcpy(table_b, table_a, num_records * 2 * sizeof(kmer_t));

  free(kmers);
  return 0;
}
