#include "common.h"
#include "io.h"
#include "lut.h"
#include "uint128_t.h"
#include "nanotime.h"

#define CPU // CUDA?

const char * USAGE = "<DSK output file>";

uint64_t swap_gt_64(uint64_t);
int compare_64(const void *, const void *);
int compare_128(const void *, const void *);

// Swaps G (11 -> 10) and T (10 -> 11) representation so radix ordering is lexical
inline uint64_t swap_gt_64(uint64_t x) {
  return x ^ ((x & 0xAAAAAAAAAAAAAAAA) >> 1);
}

void swap_64(uint64_t *x, size_t i, size_t j);
inline void swap_64(uint64_t *x, size_t i, size_t j) {
  uint64_t temp = x[i];
  x[i] = x[j];
  x[j] = temp;
}

void qsort_64(uint64_t * x, size_t l, size_t u);
void qsort_64(uint64_t * x, size_t l, size_t u) {
  size_t i, j;
  uint64_t t;
  if (l >= u)
    return;
  t = x[l];
  i = l;
  j = u+1;
  for (;;) {
    do i++; while (i <= u && x[i] < t);
    do j--; while (x[j] > t);
    if (i > j)
      break;
    swap_64(x, i, j);
  }
  swap_64(x, l, j);
  qsort_64(x, l, j-1);
  qsort_64(x, j+1, u);
}

// Doesn't reverse on bit level, reverses at the two-bit level
uint64_t block_reverse_64(uint64_t x);
uint64_t block_reverse_64(uint64_t x) {
  uint64_t output;

  unsigned char * p = (unsigned char *) &x;
  unsigned char * q = (unsigned char *) &output;
  q[7] = reverse_8(p[0]);
  q[6] = reverse_8(p[1]);
  q[5] = reverse_8(p[2]);
  q[4] = reverse_8(p[3]);
  q[3] = reverse_8(p[4]);
  q[2] = reverse_8(p[5]);
  q[1] = reverse_8(p[6]);
  q[0] = reverse_8(p[7]);
  return output;
}

uint64_t block_revcomp_64(uint64_t x);
inline uint64_t block_revcomp_64(uint64_t x) {
  uint64_t output;

  unsigned char * p = (unsigned char *) &x;
  unsigned char * q = (unsigned char *) &output;
  q[7] = revcomp_8(p[0]);
  q[6] = revcomp_8(p[1]);
  q[5] = revcomp_8(p[2]);
  q[4] = revcomp_8(p[3]);
  q[3] = revcomp_8(p[4]);
  q[2] = revcomp_8(p[5]);
  q[1] = revcomp_8(p[6]);
  q[0] = revcomp_8(p[7]);
  return output;
}

// Different to block_revcomp_64 because it shifts the correct amount after
uint64_t reverse_complement_64(uint64_t x, uint32_t k);
inline uint64_t reverse_complement_64(uint64_t x, uint32_t k) {
  return block_revcomp_64(x) >> (64 - k * 2);
}

uint128_t reverse_complement_128(uint128_t x, uint32_t k);
inline uint128_t reverse_complement_128(uint128_t x, uint32_t k) {
  uint64_t temp = block_revcomp_64(x.upper);
  x.upper = block_revcomp_64(x.lower);
  x.lower = temp;
  return right_shift_128(x, (128 - k*2));
}

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
  uint64_t a_upper = ((uint64_t*)lhs)[0];
  uint64_t a_lower = ((uint64_t*)lhs)[1];
  uint64_t b_upper = ((uint64_t*)rhs)[0];
  uint64_t b_lower = ((uint64_t*)rhs)[1];
  return (a_upper - b_upper) + (a_lower - b_lower);
}

void print_kmers_hex(FILE * outfile, uint64_t * kmers, size_t num_kmers, uint32_t kmer_num_bits);
void print_kmers_hex(FILE * outfile, uint64_t * kmers, size_t num_kmers, uint32_t kmer_num_bits) {
  assert(kmer_num_bits <= 128);
  for (size_t i = 0; i < num_kmers; i++) {
    if (kmer_num_bits == 64) {
      fprintf(outfile, "%016llx\n", kmers[i]);
    }
    else if (kmer_num_bits == 128) {
      uint64_t upper = kmers[i * 2];
      uint64_t lower = kmers[i * 2 + 1];
      fprintf(outfile, "%016llx %016llx\n", upper, lower);
    }
  }
}

void print_kmers_dec(FILE * outfile, uint64_t * kmers, size_t num_kmers, uint32_t kmer_num_bits);
void print_kmers_dec(FILE * outfile, uint64_t * kmers, size_t num_kmers, uint32_t kmer_num_bits) {
  assert(kmer_num_bits <= 128);
  for (size_t i = 0; i < num_kmers; i++) {
    if (kmer_num_bits == 64) {
      fprintf(outfile, "%020llu\n", kmers[i]);
    }
    else if (kmer_num_bits == 128) {
      uint64_t upper = kmers[i * 2];
      uint64_t lower = kmers[i * 2 + 1];
      fprintf(outfile, "%020llu %020llu\n", upper, lower);
    }
  }
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

  // allocate space for kmers
  typedef uint64_t kmer_t[kmer_num_blocks];
  kmer_t * kmers = calloc(num_records * 2, sizeof(kmer_t));

  // read items into array
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

  print_kmers_hex(stdout, (uint64_t*)kmers, num_records, kmer_num_bits);

  #ifndef NDEBUG
  TRACE("TESTING REVERSE COMPLEMENTS\n");
  if (kmer_num_bits == 64) {
    uint64_t x = ((uint64_t*)kmers)[0];
    uint64_t y = reverse_complement_64(x, k);
    TRACE("     x  = %016llx\n", x);
    TRACE("  rc(x) = %016llx\n", y);
  }
  else if (kmer_num_bits == 128) {
    uint128_t x = ((uint128_t*)kmers)[0];
    uint128_t y = reverse_complement_128(x, k);
    TRACE("     x  = %016llx %016llx\n", x.upper, x.lower);
    TRACE("  rc(x) = %016llx %016llx\n", y.upper, y.lower);
  }
  #endif

  //kmer_t * reverse_complements = kmers + num_records;
  //add_reverse_complements(kmers, reverse_complements, num_records);

  // According to the paper linked below, merge sort is better for keys of 8 bytes
  // or more
  // http://203.144.248.23/ACM.FT/1810000/1807207/p351-satish.pdf
  // and according to this, insertion sort is good for almost sorted lists:
  // http://stackoverflow.com/questions/1513566/which-sorting-algorithm-is-best-suited-to-re-sort-an-almost-fully-sorted-list
  // After the first sorting phase, we could do radix sort or insertion sort on the first char
  // TODO: Add radix sort implementation (or do my own) and time it
  // is it at least as fast as quicksort?
  // TODO: Radix sort upper digits, then insertion sort the lower
  // allocate same space
  // iterate over and count based on bytes for first positions
  // size_t counts[256];
  // for this to work, we need to internally reverse the NTs in the bytes,
  // but not inside the whole key type
  // TODO: Use GPU implementation
  // http://nvlabs.github.io/cub/classcub_1_1_block_radix_sort.html
  // TODO: time on both GPU and CPU
  // TODO: Reverse bytes first. Could do it in comparator, but probably need to
  // Examine keys multiple times
  if (kmer_num_bits == 64) {
    nanotime_t start, end;
    start = get_nanotime();
    #ifdef LIB_QSORT
    qsort(kmers, num_records, sizeof(uint64_t), compare_64);
    #else
    qsort_64((uint64_t*)kmers, 0, num_records - 1);
    #endif
    end = get_nanotime();
    double ms = (double)(end - start)/(1000000);
    #ifdef LIB_QSORT
    fprintf(stderr, "LIB QSORT timing: %.2f ms\n", ms);
    #else
    fprintf(stderr, "MY QSORT timing: %.2f ms\n", ms);
    #endif
  }
  else if (kmer_num_bits == 128) {
    // TODO: Fix compare_128
    qsort(kmers, num_records, sizeof(uint128_t), compare_128);
  }

  //print_kmers_hex(stdout, (uint64_t*)kmers, num_records, kmer_num_bits);

  // Store a copy of the kmers to sort in colex(row) order, for joining
  // in order to calculate dummy edges
  kmer_t * table_a = kmers;
  kmer_t * table_b = kmers + num_records * 2; // x2 because of reverse complements
  memcpy(table_b, table_a, num_records * 2 * sizeof(kmer_t));

  // TODO: change buffer size to 2x (then to 4x)
  // convert, then reverse all 2-bit fields
  // add reverse complements to second half of array
  // Use quicksort to sort them all
  // Memcpy, transform and resort - insertion sort based on last field?
  // filter dummy edges in two passes, keeping a list of positions
  // update positions so we can write them out
  // write out and keep edge counter, compare to dummy edge positions
  // - fill output buffer, each entry takes 5 bits of a 64 bit int
  // Find a server to run it on
  // compare time vs minia for SAME K VALUE
  // 2 data sets
  free(kmers);

  return 0;
}
