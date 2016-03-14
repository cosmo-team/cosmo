#include "io.hpp"
#include "kmer.hpp"
#include <vector>
#include <queue>

static inline uint64_t nibblet_reverse(const uint64_t &word)
{

    uint64_t ret_val = word;
    for (int i=0; i < 16; ++i) {
        //std::cout <<"************* "<< i << " ************" << std::endl << std::endl;
        //print_kmers(std::cout, &word , 1, 32);

        uint64_t right_mask = 3ull << (2*i);
        //print_kmers(std::cout, &right_mask , 1, 32);
        
        uint64_t left_mask = (3ull << 62ull) >> (2*i);
        //print_kmers(std::cout, &left_mask , 1, 32);
        
        uint64_t left_val = (word & left_mask) >> (62ull - (2*i));
        //print_kmers(std::cout, &left_val , 1, 32);
        
        uint64_t right_val = (word & right_mask) >> (2*i);
        //print_kmers(std::cout, &right_val , 1, 32);
        
        uint64_t new_left_bits = (right_val << 62ull) >> (2*i);
        //print_kmers(std::cout, &new_left_bits , 1, 32);
        
        uint64_t new_right_bits = left_val << (2*i);
        //print_kmers(std::cout, &new_right_bits , 1, 32);
        
        ret_val = (ret_val & ~left_mask) | new_left_bits;
        ret_val = (ret_val & ~right_mask) | new_right_bits;
        //std::cout << std::endl;
    }
    return ret_val;
}

std::vector<CKMCFile *>kmer_data_bases; //FIXME: move out of global
int kmc_read_header(std::string db_fname, uint32_t &kmer_num_bits, uint32_t &k, uint64 &_total_kmers)
{
    std::ifstream db_list(db_fname.c_str());
    std::string fname;
    while ( db_list >>  fname ) {
    
        CKMCFile * kmer_data_base = new CKMCFile;
        kmer_data_bases.push_back(kmer_data_base);

        if (!kmer_data_base->OpenForListing(fname))
        {

            return 0;
        }
        else
        {
            //TODO : check if the file is sorted FIXME
            
            //uint32 _kmer_length;
            uint32 _mode;
            uint32 _counter_size;
            uint32 _lut_prefix_length;
            uint32 _signature_len;
            uint32 _min_count;
            uint64 _max_count;


            kmer_data_base->Info(k, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);
            kmer_num_bits = 64 * ((k * 2 - 1) / 64 + 1); // FIXME: double check ceil(quotients) is what we're getting and not more, this may be too conservative
		
            //std::string str;


        
        }
    }
    assert(kmer_data_bases.size() > 0);
    return 1;
        
}

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
  std::cerr << "read k=" << k << " kmer_num_bits=" << *kmer_num_bits << " (which will be *64 momentarily) from cortex file." << std::endl;
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


// Compute the number of records by:
// 1. Compute record size
// 2. Reading the color table which follows the header
// 3. Compute EOF position - file position bytes
// 4. Divide remaining bytes by record size
int cortex_num_records(const int handle, const uint32_t kmer_num_bits, size_t &num_records, uint32_t &number_of_colours)
{
  off_t original_pos = -1;
  off_t end_pos = -1;
  int rc;
  rc  = read(handle, (char*)&number_of_colours, sizeof(uint32_t));
  size_t record_size = kmer_num_bits / 8 + number_of_colours * 5;

  // skip over per color header info
  uint32_t i;
  int  max_tot = 0;
  // read in (mean, total) read length column
  for (i=0; i<number_of_colours; i++) {
    int mean_read_len=0;
    rc = read(handle, &mean_read_len, sizeof(int));
    long long tot=0;
    rc = read(handle, &tot, sizeof(long long));
    if (tot > max_tot)
      max_tot = (int) tot;
  }
  // read in ID column
  for (i=0; i<number_of_colours; i++) {
    int sample_id_lens; 
    rc = read(handle, &sample_id_lens, sizeof(int));
    char tmp_name[100]; // hack // should null this out to avoid extra stack junk being printed
		rc = read(handle, tmp_name, sizeof(char) * sample_id_lens);
  }
  // read in seq_err column
  for (i=0; i<number_of_colours; i++) {
    long double seq_err;
    rc = read(handle, &seq_err, sizeof(long double));
  }
  // read in other stats??? and baseline graph name
  for (i=0; i<number_of_colours; i++) {
    char dummy;
    int d;

    rc = read(handle, &dummy, sizeof(char));
    rc = read(handle, &dummy, sizeof(char));
    rc = read(handle, &dummy, sizeof(char));
    rc = read(handle, &dummy, sizeof(char));
    rc = read(handle, &d, sizeof(int));
    rc = read(handle, &d, sizeof(int));
    int len_name_of_graph;
    rc = read(handle, &len_name_of_graph, sizeof(int));
    char name_of_graph_against_which_was_cleaned[100]; // hack
    rc = read(handle, name_of_graph_against_which_was_cleaned, sizeof(char) * len_name_of_graph);
 }

  char magic_number[6];
  rc = read(handle,  magic_number, sizeof(char)*6);

  if ( (original_pos = lseek(handle, 0 , SEEK_CUR)) == -1 ||
       (end_pos = lseek(handle, 0, SEEK_END)) == -1  ) {
    return -1;
  }
  if ( lseek(handle, original_pos, SEEK_SET) == -1 ) {
    return -1;
  }

  num_records = (end_pos - original_pos)/record_size;


  return 0;
}

void clear_bv(color_bv &bv)
{
    bv.reset();
}

void set_bit(color_bv &bv, uint32_t j)
{
    bv[j] = 1; // |= 1LL << j % 64;
}

// Only doing this complicated stuff to hopefully get rid of the counts in an efficient way
// (that is, read a large chunk of the file including the counts, then discard them)
size_t dsk_read_kmers(int handle, uint32_t kmer_num_bits, uint64_t * kmers_output, uint32_t kmer_size) {
  // TODO: Add a parameter to specify a limit to how many records we read (eventually multipass merge-sort?)
  // THIS IS ALSO A SECURITY CONCERN if we don't trust the DSK input (i.e. e.g. accept DSK files in a web service)

  // read the items items into the array via a buffer
  char input_buffer[BUFFER_SIZE];
  uint32_t maxcount = 0;
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
            uint32_t count = *((uint32_t*)(input_buffer + offset + 2*sizeof(uint64_t)));
            if (count >= maxcount) {
                maxcount = count;
                // std::cout << count << "\tk=" << kmer_size << "\t";
                // print_kmers(std::cout, kmers_output + next_slot, 1, kmer_size);
            }
                
        }
      }
    } while ( num_bytes_read );
  }
  else assert (kmer_num_bits <= 128);
  // Return the number of kmers read (whether 64 bit or 128 bit)
  return next_slot / ((kmer_num_bits/8)/sizeof(uint64_t));
}


// code from http://www.cplusplus.com/reference/queue/priority_queue/priority_queue/
// with the polarity reversed
// FIXME: don't compare strings!!!
typedef std::pair<unsigned, CKmerAPI > queue_entry;
class mycomparison
{
    bool reverse;
public:
    mycomparison(const bool& revparam=false)
        {reverse=revparam;}
    bool operator() (const queue_entry& lhs, const queue_entry&rhs) const
        {
            if (reverse) return (const_cast<CKmerAPI*>(&(lhs.second))->to_string() > const_cast<CKmerAPI*>(&(rhs.second))->to_string());
            else return (const_cast<CKmerAPI*>(&(lhs.second))->to_string() < const_cast<CKmerAPI*>(&(rhs.second))->to_string());
        }
};

size_t kmc_read_kmers(const int handle, const uint32_t kmer_num_bits, const uint32_t num_colors, uint32_t k, uint64_t *const &kmers_output, std::vector<color_bv>  &kmer_colors)
{

    typedef std::priority_queue<queue_entry, std::vector<queue_entry>, mycomparison> mypq_type;
    mypq_type queue;
    int numkmers = 0;
    color_bv color = 0;
    
    // initialize the queue with a file identifier (as a proxy for the input sequence itself) and the value at the head of the file for peeking
    for (unsigned i = 0; i < kmer_data_bases.size(); ++i) {
        CKmerAPI kmer_object(k);
        uint64 counter;// for coverage
        if (!kmer_data_bases[i]->ReadNextKmer(kmer_object, counter)) {
            std::cerr << "ERROR: empty file contains no kmers." << std::endl;
            assert(false);
        }
        queue.push(std::make_pair(i, kmer_object));
    }

    // pop the first element into 'current' to initialize our state (and init any other state here such as this one's color)
    queue_entry current = queue.top();
    queue.pop();
    color.set(current.first); // FIXME: make sure not using << operator elsewhere!
    //and replace that streams entry
    CKmerAPI kmer_object(k);
    uint64 counter;// for coverage    
    if (kmer_data_bases[current.first]->ReadNextKmer(kmer_object, counter)) {
        queue.push(std::make_pair(current.first, kmer_object));
    }


    // 
    while (!queue.empty()) {
        if (const_cast<CKmerAPI*>(&(queue.top().second))->to_string() == const_cast<CKmerAPI*>(&(current.second))->to_string()) { // if this is the same kmer we've seen before
            queue_entry additional_instance = queue.top();
            queue.pop();
            CKmerAPI kmer_object(k);
            uint64 counter;// for coverage
            color.set(additional_instance.first);
            if (kmer_data_bases[current.first]->ReadNextKmer(kmer_object, counter)) {
                queue.push(std::make_pair(additional_instance.first, kmer_object));
            }
        } else { // if the top of the queue contains a new instance

            // emit our current state
            std::vector<unsigned long long /*uint64*/> kmer;            
            current.second.to_long(kmer);

            kmer_colors[numkmers] = color;
            color.reset();

            for (int block=0; block < kmer.size(); ++block) {
                kmers_output[numkmers*kmer.size() + block] = kmer[block]; // FIXME: check if kmer_output is big endian or little endian

                if (block == 1) {
                    assert(k == 63); //FIXME: the following line is to fix a bug in the kmc2 API where it shifts word[0] << 64 when k=63 which results in "word[1] =  word[0] + word[1]" the next line compensates; not sure how pervasive this is in the k>32 space.
                    kmers_output[numkmers*kmer.size() + block] -=kmer[0];
                }

            }
            numkmers++;

            // now initialize our current state with the top
            current = queue.top();
            queue.pop();
            color.set(current.first); // FIXME: make sure not using << operator elsewhere!
            //and replace that streams entry
            CKmerAPI kmer_object(k);
            uint64 counter;// for coverage    
            if (kmer_data_bases[current.first]->ReadNextKmer(kmer_object, counter)) {
                queue.push(std::make_pair(current.first, kmer_object));
            }
        }
    }

    // and finally emit our current state    
    std::vector<unsigned long long /*uint64*/> kmer;            
    current.second.to_long(kmer);
        
    kmer_colors[numkmers] = color;
    color.reset();
        
    for (int block=0; block < kmer.size(); ++block) {
        kmers_output[numkmers*kmer.size() + block] = kmer[block]; // FIXME: check if kmer_output is big endian or little endian
            
        if (block == 1) {
            assert(k == 63); //FIXME: the following line is to fix a bug in the kmc2 API where it shifts word[0] << 64 when k=63 which results in "word[1] =  word[0] + word[1]" the next line compensates; not sure how pervasive this is in the k>32 space.
            kmers_output[numkmers*kmer.size() + block] -=kmer[0];
        }
            
    }
    numkmers++;
            
    

		
    // CKmerAPI kmer_object(k);



        
    // while (kmer_data_bases[0]->ReadNextKmer(kmer_object, counter)) {


    //     std::vector<unsigned long long /*uint64*/> kmer;            
    //     kmer_object.to_long(kmer);


    //     color_bv color = 1;
    //     kmer_colors[numkmers] = color;

    //     for (int block=0; block < kmer.size(); ++block) {
    //         kmers_output[numkmers*kmer.size() + block] = kmer[block]; // FIXME: check if kmer_output is big endian or little endian

    //         if (block == 1) {
    //             assert(k == 63); //FIXME: the following line is to fix a bug in the kmc2 API where it shifts word[0] << 64 when k=63 which results in "word[1] =  word[0] + word[1]" the next line compensates; not sure how pervasive this is in the k>32 space.
    //             kmers_output[numkmers*kmer.size() + block] -=kmer[0];
    //         }

    //     }
    //     numkmers++;
    // }

    for (auto kmer_data_base: kmer_data_bases) {
        kmer_data_base->Close();
        delete kmer_data_base;
    }
    return numkmers;
    
}

size_t cortex_read_kmers(const int handle, const uint32_t kmer_num_bits, const uint32_t num_colors, uint32_t k, uint64_t *const &kmers_output, std::vector<color_bv>  &kmer_colors) {
    // TODO: Add a parameter to specify a limit to how many records we read (eventually multipass merge-sort?)
    (void) k;
    size_t next_slot = 0;
    int coverage[NUM_COLS]; // hack
    unsigned long long maxcount = 0;
    int covsum = 0;
    while (1) {
        unsigned int i;
        uint64_t kmer;
        char individual_edges_reading_from_binary[NUM_COLS]; // bit field for which edges enter and leave
        int rc;

        rc = read(handle, &kmer, sizeof(uint64_t));
        if (rc <= 0)
            break;
        rc = read(handle, &coverage, sizeof(int) * num_colors);
        rc = read(handle, individual_edges_reading_from_binary, sizeof(char) * num_colors);

        // find the max coverage for this coverage
        covsum = 0;
        for (int covit=0; covit < num_colors; ++covit) {
                covsum += coverage[covit];
        }
        char edge = 0;
        for (i=0; i<num_colors; i++) {
            edge |= individual_edges_reading_from_binary[i] & 0xF;
        }
        kmers_output[next_slot] = kmer;

        if (covsum >= maxcount) {
            maxcount = covsum;
            std::cout << maxcount << "\t" << kmers_output[next_slot] << "\t";
            
            print_kmers(std::cout, kmers_output + next_slot, 1, k);
        }

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
