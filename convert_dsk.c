#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>

#include "debug.h"
//#include "nanotime.h"

const char * USAGE = "<DSK output file>";
const size_t MAX_BITS_PER_KMER = 128;
const size_t BUFFER_SIZE = 0x8000; // 32Kb data buffer

int dsk_read_header(int, uint32_t *, uint32_t *);
uint64_t swap_gt_64(uint64_t);
int compare_64(const void *, const void *);
int compare_128(const void *, const void *);

inline uint64_t swap_gt_64(uint64_t x) {
  return x ^ ((x & 0xAAAAAAAAAAAAAAAA) >> 1);
}

int compare_64(const void * lhs, const void * rhs)
{
  uint64_t a = *(uint64_t*)lhs;
  uint64_t b = *(uint64_t*)rhs;
  return a - b;
}

int compare_128(const void * lhs, const void * rhs) {
  uint64_t a_upper = ((uint64_t*)lhs)[0];
  uint64_t a_lower = ((uint64_t*)lhs)[1];
  uint64_t b_upper = ((uint64_t*)rhs)[0];
  uint64_t b_lower = ((uint64_t*)rhs)[1];
  return (a_upper - b_upper) + (a_lower - b_lower);
}

int dsk_read_header(int handle, uint32_t * kmer_num_bits, uint32_t * k) {
  int success = read(handle, (char*)kmer_num_bits, sizeof(uint32_t)) != -1 &&
                read(handle, (char*)k, sizeof(uint32_t)) != -1;
  return success;
}

size_t dsk_record_size(uint32_t);
inline size_t dsk_record_size(uint32_t kmer_num_bits) {
  // record fmt: kmer, count
  return (kmer_num_bits/8) + 4;
}

int dsk_num_records(int handle, uint32_t kmer_num_bits, size_t * num_records);
int dsk_num_records(int handle, uint32_t kmer_num_bits, size_t * num_records) {
  size_t record_size = dsk_record_size(kmer_num_bits);
  off_t original_pos = -1;
  off_t end_pos = -1;
  if ( (original_pos = lseek(handle, 0 ,SEEK_CUR)) == -1 ||
       (end_pos = lseek(handle, 0, SEEK_END)) == -1  ) {
    return -1;
  }
  if ( lseek(handle, original_pos, SEEK_SET) == -1 ) {
    return -1;
  }
  *num_records = (end_pos - original_pos)/record_size;
  return 0;
}

size_t dsk_read_kmers(int handle, uint32_t kmer_num_bits, uint64_t * kmers_output);
size_t dsk_read_kmers(int handle, uint32_t kmer_num_bits, uint64_t * kmers_output) {
  // TODO: Add a parameter to specify a limit to how many records we read (eventually multipass merge-sort?)

  // read the items items into the array via a buffer
  char input_buffer[BUFFER_SIZE];

  // Only read a multiple of records... this makes it easier to iterate over
  size_t record_size = dsk_record_size(kmer_num_bits);
  size_t read_size = (BUFFER_SIZE / record_size) * record_size;

  // Technically we *could* read in bytes at a time if need-be... but nah.
  assert(read_size > 0);

  ssize_t num_bytes_read = 0;
  size_t next_slot = 0;

  // This if statement would be more readable inside the loop, but it's moved out here for performance.
  // TODO: This might run better if we filled in another buffer on the stack with kmers only, then memcpy to the output?
  if (kmer_num_bits <= 64) {
    do {
      // Try read a batch of records.
      if ( (num_bytes_read = read(handle, input_buffer, read_size)) == -1 ) {
        return -1;
      }

      // Did we read anything?
      if (num_bytes_read ) {
        TRACE("num_bytes_read = %zd\n", num_bytes_read);

        // Iterate over kmers, skipping counts
        for (ssize_t offset = 0; offset < num_bytes_read; offset += sizeof(uint64_t) + sizeof(uint32_t), next_slot += 1) {
          kmers_output[next_slot] = *((uint64_t*)(input_buffer + offset));
        }
      }
    } while ( num_bytes_read );
  }
  else if (64 < kmer_num_bits <= 128) {
    do {
      // Try read a batch of records.
      if ( (num_bytes_read = read(handle, input_buffer, read_size)) == -1 ) {
        return -1;
      }

      // Did we read anything?
      if (num_bytes_read ) { 
        TRACE("num_bytes_read = %zd\n", num_bytes_read);

        // Iterate over kmers, skipping counts
        for (ssize_t offset = 0; offset < num_bytes_read; offset += 2 * sizeof(uint64_t) + sizeof(uint32_t), next_slot += 2) {
            // Swapping lower and upper block (to simplify sorting later)
            kmers_output[next_slot + 1]     = *((uint64_t*)(input_buffer + offset));
            kmers_output[next_slot] = *((uint64_t*)(input_buffer + offset + sizeof(uint64_t)));
        }
      }
    } while ( num_bytes_read );
  }
  else assert (kmer_num_bits <= 128);
  // Return the number of kmers read (whether 64 bit or 128 bit)
  return next_slot / ((kmer_num_bits/8)/sizeof(uint64_t));
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
  kmer_t * kmers = malloc(num_records * sizeof(kmer_t));

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

  // print em
  for (size_t i = 0; i < num_records; i++)
  {
    if (kmer_num_blocks == 1) {
      printf("%016llx\n", *kmers[i]);
    }
    else if (kmer_num_blocks == 2) {
      uint64_t upper = kmers[i][0];
      uint64_t lower = kmers[i][1];
      printf("%016llx %016llx\n", upper, lower);
    }
  }
  /*
  printf("SORTING...\n");
  if (kmer_num_blocks == 1) {
    qsort(kmers, num_records, sizeof(uint64_t), compare_64);
  }
  else if (kmer_num_blocks == 2) {
    qsort(kmers, num_records, 2 * sizeof(uint64_t), compare_128);
  }

  for (size_t i = 0; i < num_records; i++)
  {
    if (kmer_num_blocks == 1) {
      printf("%016llx\n", *kmers[i]);
    }
    else if (kmer_num_blocks == 2) {
      uint64_t upper = kmers[i][0];
      uint64_t lower = kmers[i][1];
      printf("%016llx %016llx\n", lower, upper);
    }
  }
  */
  // TODO: refactor
  //
  // Convert to normal format
  // change buffer size to 2x
  // iterate over and read them into buffer
  // add their reverse complements
  // Convert
  // Sort
  // Filter dummy edges
  // - count how many
  // - create array
  // - create dummy edges
  // Sort
  // Filter for dummy edges
  // Merge
  // Write
  free(kmers);

  return 0;
}
