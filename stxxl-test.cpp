#include <fstream>      // for std::fstream
#include <stxxl.h>      // STXXL header
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>
//#include <stxxl/bits/algo/stable_ksort.h>

#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "utility.hpp"
//#include "dummies.hpp"


// TODO: Compile flag check
// typedef uint64_t kmer_t;
typedef uint128_t kmer_t;

const uint64_t MB_TO_BYTES = 1024 * 1024;

int main(int argc, char* argv[])
{
  using namespace boost::adaptors;

  // Parameter extraction
  // default 5gb for mem
  stxxl::internal_size_type M = 5 * 1024 * MB_TO_BYTES;
  std::string file_name = argv[1];
  uint32_t k = atol(argv[2]);

  // Open the file
  std::ifstream in(file_name, std::ifstream::binary);

  // Convert to our format: reverse for colex ordering, swap g/t encoding
  auto input = make_typed_input_range<kmer_t>(in)
             | transformed(swap_gt<kmer_t>())
             | transformed(reverse_nt<kmer_t>());

  // Create STXXL vector
  stxxl::vector<kmer_t> kmers;

  // Load vector with kmers and reverse complements
  auto revcomp = reverse_complement<kmer_t>(k);
  // TODO: Try out STXXL writebuf approach for async I/O
  // (Mainly useful if using multiple disks)
  for (auto x : input) {
    kmers.push_back(x);
    kmers.push_back(revcomp(x));
  }
  std::cerr << "Read " << kmers.size() << " kmers from file."<< std::endl;

  // Sort them in colex(node)-edge order
  std::cerr << "Sorting..." << std::endl;
  stxxl::ksort(kmers.begin(), kmers.end(), get_key_colex_node<kmer_t>(), M);

  // Get another copy of the table sorted by out-node (for dummy edge discovery)
  stxxl::vector<kmer_t> kmers_by_end_node(kmers);
  stxxl::ksort(kmers_by_end_node.begin(), kmers_by_end_node.end(), get_key_colex_edge<kmer_t>(), M);


  // Find incoming dummy edges
  std::cerr << "Searching for incoming dummies..." << std::endl;
  stxxl::vector<kmer_t> incoming_dummies;
  find_incoming_dummy_edges(kmers.begin(), kmers.end(),
                            kmers_by_end_node.begin(), kmers_by_end_node.end(),
                            k, std::back_inserter(incoming_dummies));
  std::cerr << "Found " << incoming_dummies.size() << " incoming dummies." << std::endl;

  //stxxl::for_each(kmers.begin(), kmers.end(),
  //    [k](const kmer_t x){ std::cout << kmer_to_string(x, k) << std::endl; }, M);
  //    [k](const kmer_t x){ std::cout << kmer_to_string(x, k, k) << std::endl; }, M);

    /*
     * if it changes the outcome, filter for aaaa and set a flag for aaaa...
     * - make small .fa file so i can control which dummy edges there are...
     *
     * TEST get_start_node, get_end_node on single kmer
     */

  return 0;
}
