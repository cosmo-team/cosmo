#include <algorithm>    // for STL std::sort
#include <vector>       // for STL std::vector
#include <fstream>      // for std::fstream
#include <stxxl.h>      // STXXL header
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/algorithm/copy.hpp>

#include "uint128_t.hpp"
#include "kmer.hpp"
#include "utility.hpp"
#include "io.hpp"

// Compile flag check
// typedef uint64_t kmer_t;
typedef uint128_t kmer_t;
typedef stxxl::vector<kmer_t> kmer_vector_t;

const uint64_t MB_TO_BYTES = 1024 * 1024;
const uint64_t INPUT_BUFFER_SIZE = 1 * MB_TO_BYTES;

struct comp: public std::less<kmer_t> // ascending order
{
  bool operator()(const kmer_t& lhs, const kmer_t& rhs) const {
    kmer_t lhs_node = get_start_node(lhs);
    kmer_t rhs_node = get_start_node(rhs);
    if (lhs_node == rhs_node) {
      kmer_t lhs_edge = get_edge_label(lhs);
      kmer_t rhs_edge = get_edge_label(rhs);
      return (lhs_edge < rhs_edge);
    }
    else return (lhs_node < rhs_node);
  }

  static kmer_t min_value()
  { return kmer_t(); }

  static kmer_t max_value()
    // Have to change for uint64_t
  { return kmer_t(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max()); }
};

int main(int argc, char* argv[])
{
    // Parameter extraction
    // default 5gb for mem
    stxxl::internal_size_type M = 5 * 1024 * MB_TO_BYTES;
    std::string file_name = argv[1];
    uint32_t k = 56; // make parameter

    // Open the file
    std::ifstream in(file_name, std::ifstream::binary);

    // Convert to our format: reverse for colex ordering, swap g/t encoding
    auto input = get_input_range<kmer_t>(in)
               | transformed(swap_gt<kmer_t>())
               | transformed(reverse_nt<kmer_t>());

    // Create STXXL vector
    kmer_vector_t kmers;

    // Copy kmers + reverse complements into stxxl vector
    // if we don't want reverse complements: boost::copy(input, std::back_inserter(kmers));
    auto revcomp = reverse_complement<kmer_t>(k);
    for (kmer_t x : input) {
      kmers.push_back(x);
      kmers.push_back(revcomp(x));
    }

    // Sort them in colex(node)-edge order
    std::cerr << "Sorting..." << std::endl;
    // const stxxl::internal_size_type M = atol(argv[2]) * 1024 * 1024;
    stxxl::sort(kmers.begin(), kmers.end(), comp(), M);
    // stxxl::ksort(kmers.begin(), kmers.end(), M); // cant use ksort with 128bit ints

    std::cerr << kmers.size() << std::endl;

    // Output - use faster methods?
    //stxxl::for_each(kmers.begin(), kmers.end(), [k](const kmer_t x){ std::cout << kmer_to_string(x, k) << std::endl; }, M);

    return 0;
}
