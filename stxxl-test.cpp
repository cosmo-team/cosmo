#include <fstream>      // for std::fstream
#include <stxxl.h>      // STXXL header
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>
//#include <stxxl/bits/algo/stable_ksort.h>

#include "integers.hpp"
#include "kmer.hpp"
//#include "kmer.hpp"
//#include "dummies.hpp"
//#include "utility.hpp"
//#include "io.hpp"

using namespace cosmo;

// TODO: Compile flag check
// typedef uint64_t kmer_t;
typedef uint128_t kmer_t;

const uint64_t MB_TO_BYTES = 1024 * 1024;

/*
struct comp {
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

struct get_key_colex_node
{
  typedef __uint128_t key_type;
  key_type operator() (const kmer_t & obj) const {
    kmer_t temp = get_start_node(obj);
    uint64_t edge = get_edge_label(obj)>>30;
    temp._lower |= edge;
    uint64_t temp_block = temp._upper;
    temp._upper = temp._lower;
    temp._lower = temp_block;
    return *((key_type*)&temp);
  }
  key_type min_value() const
  { return std::numeric_limits<key_type>::min(); }
  key_type max_value() const
  // Have to change for uint64_t
  { return std::numeric_limits<key_type>::max(); }
};

struct get_key_colex_edge
{
  typedef __uint128_t key_type;
  key_type operator() (const kmer_t & obj) const {
    kmer_t temp;
    temp._lower = obj._upper;
    temp._upper = obj._lower;
    return *((key_type*)&temp);
  }
  key_type min_value() const
  { return std::numeric_limits<key_type>::min(); }
  key_type max_value() const
  // Have to change for uint64_t
  { return std::numeric_limits<key_type>::max(); }
};
*/

template <typename kmer_t>
std::pair<std::istream_iterator<kmer_t>, std::istream_iterator<kmer_t>>
get_input_range(std::ifstream & in) {
  auto in_begin = std::istream_iterator<kmer_t>(in);
  auto in_end   = std::istream_iterator<kmer_t>();
  auto in_range = std::make_pair(in_begin, in_end);
  return in_range;
}

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
  auto input = get_input_range<kmer_t>(in)
             | transformed(swap_gt<kmer_t>())
             | transformed(reverse_nt<kmer_t>());

  // Create STXXL vector
  stxxl::vector<kmer_t> kmers;

  //auto revcomp = reverse_complement<kmer_t>(k);
  for (auto x : input) {
    kmers.push_back(x);
    //kmers.push_back(revcomp(x));
  }
  std::cerr << "Read " << kmers.size() << " kmers from file."<< std::endl;

  // Sort them in colex(node)-edge order
  //std::cerr << "Sorting..." << std::endl;
  //stxxl::ksort(kmers.begin(), kmers.end(), get_key_colex_node(), M);

  // Get another copy of the table sorted by out-node (for dummy edge discovery)
  //kmer_vector_t kmers_by_end_node(kmers);
  //stxxl::ksort(kmers_by_end_node.begin(), kmers_by_end_node.end(), get_key_colex_edge(), M);


  // Find incoming dummy edges
  //kmer_vector_t incoming_dummies;
  //find_incoming_dummy_edges(kmers.begin(), kmers.end(),
  //                          kmers_by_end_node.begin(), kmers_by_end_node.end(),
  //                          k, std::back_inserter(incoming_dummies));

  //stxxl::for_each(kmers.begin(), kmers.end(),
  //    [k](const kmer_t x){ std::cout << kmer_to_string(x, k, 64) << std::endl; }, M);

  //std::cout << incoming_dummies.size() << std::endl;
  //stxxl::for_each(incoming_dummies.begin(), incoming_dummies.end(),
  //    [k](const kmer_t x){ std::cout << kmer_to_string(x, k, k) << std::endl; }, M);

    /*
     * aaaa is repeated... maybe the sorting is messing it up.
     * Try with normal sort?
     * if it changes the outcome, filter for aaaa and set a flag for aaaa...
     * - support 64bit... (compile flag)
     * - make small .fa file so i can control which dummy edges there are...
     * - (also have aaaa version to test the sorting thing)
     *
     * TEST get_start_node, get_end_node on single kmer
     */
  kmer_t x = kmers[0];
  //std::cerr << x.table[0] << std::endl;
  //std::cerr << x.table[1] << std::endl;
  std::cerr << std::setw(16) << std::setfill('0') << std::hex << ((uint64_t*)&x)[0] << std::endl;
  std::cerr << std::setw(16) << std::hex << ((uint64_t*)&x)[1] << std::endl;
  //kmer_t y = revcomp(kmers[0]);
  std::cerr << kmer_to_string(x, k) << std::endl;
  //std::cerr << kmer_to_string(y, k, k) << std::endl;
  // the later ones are MSB -> don't reverse them
  // a.table[0] = b.table[1] = 1 -> (a < b)

  return 0;
}
