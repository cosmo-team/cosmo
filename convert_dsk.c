#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

#include "debug.h"
//#include "nanotime.h"

const char * USAGE = "<DSK output file>";

const size_t BUFFER_SIZE = 0x8000; // 32Kb data buffer

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
  if ( (read(handle, (char*)&kmer_num_bits, sizeof(uint32_t))) == -1 ||
       (read(handle, (char*)&k, sizeof(uint32_t))) == -1 ) {
    fprintf(stderr, "Error reading file %s\n", file_name);
    exit(1);
  }
  TRACE("kmer_num_bits, k = %d, %d\n", kmer_num_bits, k);
  uint32_t kmer_num_blocks = (kmer_num_bits / 8) / sizeof(uint64_t);
  TRACE("kmer_num_blocks = %d\n", kmer_num_blocks);

  // read how many items there are
  size_t record_size = (kmer_num_bits / 8) + 4; // record fmt: kmer, count
  off_t original_pos = -1;
  off_t end_pos = -1;
  if ( (original_pos = lseek(handle, 0 ,SEEK_CUR)) == -1 ||
       (end_pos = lseek(handle, 0, SEEK_END)) == -1  ) {
    fprintf(stderr, "Error seeking file %s\n", file_name);
    exit(1);
  }
  if ( lseek(handle, original_pos, SEEK_SET) == -1 ) {
    fprintf(stderr, "Error seeking file %s\n", file_name);
    exit(1);
  }
  size_t num_records = (end_pos - original_pos)/record_size;
  TRACE("num_records = %zu\n", num_records);

  // allocate an array for input
  if (kmer_num_bits > 64) {
    uint64_t b[2];
    read(handle, (char*)&b, sizeof(uint64_t) * 2);
    TRACE("%016llx%016llx\n", b[0], b[1]);
    fprintf(stderr, "Not yet implemented for k > 32\n");
    exit(1);
  }
  typedef uint64_t kmer_t[kmer_num_blocks];
  kmer_t * kmers = calloc(num_records, sizeof(kmer_t));

  // read the items items into the array via a buffer
  char input_buffer[BUFFER_SIZE];
  memset(input_buffer, 0, BUFFER_SIZE);
  // Only read a multiple of records... we are going to remove the counts anyway
  size_t read_size = BUFFER_SIZE / record_size;
  ssize_t num_bytes_read = 0;
  do {
    if ( (num_bytes_read = read(handle, input_buffer, read_size)) == -1 ) {
      fprintf(stderr, "Error reading file %s\n", argv[1]);
      exit(1);
    }
    // TODO: iterate over, adding the kmers, as curr_byte < num_bytes_read...
    // TODO: make function that fills the buffer 
    // curr_byte += sizeof(uint64_t) * 2
    // curr_byte += sizeof(uint32_t)
    if (num_bytes_read ) { 
      TRACE("%016llx\n", (uint64_t)((uint64_t*)input_buffer)[0]);
    }
  } while ( num_bytes_read );

  close(handle);

  // TODO: version control
  // refactor
  // change buffer size to 2x
  // iterate over and read them into buffer
  // + their reverse complements
  // print them out to make sure its being read
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
