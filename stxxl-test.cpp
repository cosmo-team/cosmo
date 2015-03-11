#include <typeinfo>
#include <fstream>
#include <stxxl.h>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>

#include "utility.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"

// TODO: Compile flag check
// typedef uint64_t kmer_t;
typedef uint128_t kmer_t;
typedef std::pair<kmer_t, size_t>  dummy_t;

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
  for (auto x : input) {
    kmers.push_back(x);
    kmers.push_back(revcomp(x));
  }
  std::cerr << "Read " << kmers.size()/2 << " kmers from file."<< std::endl;

  // Sort them in colex(node)-edge order
  std::cerr << "Sorting..." << std::endl;
  stxxl::ksort(kmers.begin(), kmers.end(), get_key_colex_node<kmer_t>(), M);

  // Get another copy of the table sorted by out-node (for dummy edge discovery)
  stxxl::vector<kmer_t> kmers_by_end_node(kmers);
  stxxl::ksort(kmers_by_end_node.begin(), kmers_by_end_node.end(), get_key_colex_edge<kmer_t>(), M);

  // Find nodes that require incoming dummy edges
  std::cerr << "Searching for nodes requiring incoming dummy edges..." << std::endl;
  stxxl::vector<dummy_t> incoming_dummies;
  find_incoming_dummy_edges(kmers, kmers_by_end_node, k, incoming_dummies);
  std::cerr << "Added " << incoming_dummies.size() << " incoming dummy edges." << std::endl;

  // Sort dummies
  stxxl::sort(incoming_dummies.begin(), incoming_dummies.end(), colex_dummy_less<dummy_t>(), M);

  // Merge dummies and output

  return 0;
}
