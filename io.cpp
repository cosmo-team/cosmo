#include "io.hpp"

int dsk_read_header(int handle, uint32_t * kmer_num_bits, uint32_t * k) {
  return read(handle, (char*)kmer_num_bits, sizeof(uint32_t)) != -1 &&
         read(handle, (char*)k, sizeof(uint32_t)) != -1;
}

int cortex_read_header(int handle, uint32_t * kmer_num_bits, uint32_t * k) {
  // check header
  char magic_number[6];
  int version;
  int rc;
  rc = read(handle, &magic_number, sizeof(char) * 6);

  rc = read(handle, &version,sizeof(int));
  if (read < 0 || strncmp(magic_number, "CORTEX", 6))
    return 0;

  read(handle, (char*)k ,sizeof(uint32_t));
  read(handle, (char*)kmer_num_bits, sizeof(uint32_t));
  *kmer_num_bits *= 64;
  return 1;

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
  *num_records = (end_pos - original_pos)/record_size * 4; // worst case there are 4 edges for each kmer
  return 0;
}

int cortex_num_records(int handle, uint32_t kmer_num_bits, size_t * num_records, uint32_t * num_colors) {
  int number_of_colours;
  off_t original_pos = -1;
  off_t end_pos = -1;
  int rc;
  rc  = read(handle, (char*)&number_of_colours, sizeof(uint32_t));
  size_t record_size = kmer_num_bits / 8 + number_of_colours * 5;

  // skip over per color header info
  int i;
  int  max_tot = 0;
  for (i=0; i<number_of_colours;i++) {
    int mean_read_len=0;
    rc = read(handle,&mean_read_len, sizeof(int));
    long long tot=0;
    rc = read(handle, &tot, sizeof(long long));
    if (tot > max_tot)
      max_tot = (int) tot;
  }
  for (i=0; i<number_of_colours;i++) {
    int sample_id_lens; 
    rc = read(handle,&sample_id_lens,sizeof(int));
    char tmp_name[100]; // hack // should null this out to avoid extra stack junk being printed
		rc = read(handle,tmp_name,sizeof(char) * sample_id_lens);
  }

  for (i=0; i<number_of_colours;i++) {
    long double seq_err;
    rc = read(handle,&seq_err,sizeof(long double));
  }

  for (i=0; i<number_of_colours;i++) {
    char dummy;
    int d;

    rc = read(handle,&dummy,sizeof(char));
    rc = read(handle,&dummy,sizeof(char));
    rc = read(handle,&dummy,sizeof(char));
    rc = read(handle,&dummy,sizeof(char));
    rc = read(handle,&d,sizeof(int));
    rc = read(handle,&d,sizeof(int));
    int len_name_of_graph;
    rc = read(handle,&len_name_of_graph,sizeof(int));
    char name_of_graph_against_which_was_cleaned[100]; // hack
    rc = read(handle,name_of_graph_against_which_was_cleaned,sizeof(char)*len_name_of_graph);
 }

  char magic_number[6];
  rc = read(handle, magic_number,sizeof(char)*6);

  if ( (original_pos = lseek(handle, 0 ,SEEK_CUR)) == -1 ||
       (end_pos = lseek(handle, 0, SEEK_END)) == -1  ) {
    return -1;
  }
  if ( lseek(handle, original_pos, SEEK_SET) == -1 ) {
    return -1;
  }

  *num_records = (end_pos - original_pos)/record_size;
  *num_colors = number_of_colours;

  return 0;
}

void clear_bv(color_bv &bv)
{
    bv = 0;
}

void set_bit(color_bv &bv, uint32_t j)
{
    bv |= 1LL << j % 64;
}

// Only doing this complicated stuff to hopefully get rid of the counts in an efficient way
// (that is, read a large chunk of the file including the counts, then discard them)
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


size_t cortex_read_kmers(int handle, uint32_t kmer_num_bits, uint32_t num_colors, uint32_t k, uint64_t * kmers_output, std::vector<color_bv>  &kmer_colors) {
    // TODO: Add a parameter to specify a limit to how many records we read (eventually multipass merge-sort?)
    (void) k;
    size_t next_slot = 0;
    int coverage[100]; // hack
    while (1) {
        unsigned int i;
        uint64_t kmer;
        char individual_edges_reading_from_binary[100]; // bit field for which edges enter and leave
        int rc;

        rc = read(handle, &kmer,sizeof(uint64_t));
        if (rc <= 0)
            break;
        rc = read(handle, &coverage, sizeof(int) * num_colors);
        rc = read(handle, individual_edges_reading_from_binary, sizeof(char) * num_colors);
        char edge = 0;
        for (i=0; i<num_colors;i++) {
            edge |= individual_edges_reading_from_binary[i] & 0xF;
        }
        kmers_output[next_slot] = kmer;

        color_bv color_acc;
        uint32_t j;
        // A (0) -> 0001, C (1) -> 0010, G (2) -> 0100, T (3) -> 1000
        for (i=0; i< 4; i++) {
            char mask = 1 << i;
            if (edge & mask) {
                /*
                  if (kmer_num_bits <= 64) {
                  uint64_t k_plus_1 = kmer | i << (k * 2);
                  kmers_output[next_slot] = k_plus_1;
                  }
                */
                // now write out whether each color has this kmer edge
                clear_bv(color_acc);
                for (j=0; j<num_colors;j++) {
                    if (individual_edges_reading_from_binary[j] & mask)
                        set_bit(color_acc, j);
                        //color_acc |= 1LL << j % 64;
                    //if ((j > 0 && j % 64 == 0) || j == num_colors -1) {
                        // write out bits when we have filled accumulator
                    //}
                }
                kmer_colors[next_slot] = color_acc;
                clear_bv(color_acc);


            }
        }
        next_slot++;
    }
    return next_slot / ((kmer_num_bits/8)/sizeof(uint64_t));
}
