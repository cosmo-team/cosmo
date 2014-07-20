#include "io.hpp"

int dsk_read_header(int handle, uint32_t * kmer_num_bits, uint32_t * k) {
  return read(handle, (char*)kmer_num_bits, sizeof(uint32_t)) != -1 &&
         read(handle, (char*)k, sizeof(uint32_t)) != -1;
}

int dsk_num_records(int handle, uint32_t kmer_num_bits, size_t * num_records) {
  size_t record_size = DSK_FILE_RECORD_SIZE(kmer_num_bits);
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

size_t dsk_read_kmers(int handle, uint32_t kmer_num_bits, uint64_t * kmers_output) {
  // TODO: Add a parameter to specify a limit to how many records we read (eventually multipass merge-sort?)
  // THIS IS ALSO A SECURITY CONCERN if we don't trust the DSK input (i.e. e.g. accept DSK files in a web service)

  // read the items items into the array via a buffer
  char input_buffer[BUFFER_SIZE];

  // Only read a multiple of records... this makes it easier to iterate over
  size_t record_size = DSK_FILE_RECORD_SIZE(kmer_num_bits);
  size_t read_size = (BUFFER_SIZE / record_size) * record_size;
  assert(read_size > 0);

  ssize_t num_bytes_read = 0;
  size_t next_slot = 0;

  // This if statement would be more readable inside the loop, but it's moved out here for performance.
  if (kmer_num_bits <= 64) {
    do {
      // Try read a batch of records.
      if ( (num_bytes_read = read(handle, input_buffer, read_size)) == -1 ) {
        return 0;
      }

      // Did we read anything?
      if (num_bytes_read ) {
        // Iterate over kmers, skipping counts
        for (ssize_t offset = 0; offset < num_bytes_read; offset += sizeof(uint64_t) + sizeof(uint32_t), next_slot += 1) {
          kmers_output[next_slot] = *((uint64_t*)(input_buffer + offset));
        }
      }
    } while ( num_bytes_read );
  }
  else if (64 < kmer_num_bits && kmer_num_bits <= 128) {
    do {
      // Try read a batch of records.
      if ( (num_bytes_read = read(handle, input_buffer, read_size)) == -1 ) {
        return 0;
      }

      // Did we read anything?
      if (num_bytes_read ) {
        // Iterate over kmers, skipping counts
        for (ssize_t offset = 0; offset < num_bytes_read; offset += 2 * sizeof(uint64_t) + sizeof(uint32_t), next_slot += 2) {
            // Swapping lower and upper block (to simplify sorting later)
            kmers_output[next_slot + 1] = *((uint64_t*)(input_buffer + offset));
            kmers_output[next_slot]     = *((uint64_t*)(input_buffer + offset + sizeof(uint64_t)));
        }
      }
    } while ( num_bytes_read );
  }
  else assert (kmer_num_bits <= 128);
  // Return the number of kmers read (whether 64 bit or 128 bit)
  return next_slot / ((kmer_num_bits/8)/sizeof(uint64_t));
}
