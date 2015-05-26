#include <typeinfo>
#include <fstream>
#include <libgen.h>
#include <stxxl.h>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>

// TCLAP
#include "tclap/CmdLine.h"

#include "utility.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"

// TODO: Compile flag check
// typedef uint64_t kmer_t;
typedef uint128_t kmer_t;
typedef std::pair<kmer_t, size_t>  dummy_t;
const static string extension = ".packed";

typedef struct p
{
    //bool ascii = false;
    std::string input_filename = "";
    std::string output_prefix = "";
    size_t k = 0;
    size_t m = 0;
} parameters_t;

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  /* // Add this option after refactoring the visitors (for now just compile with DEBUG if you want printed edges)
  TCLAP::SwitchArg ascii_arg("a", "ascii",
            "Outputs *full* edges (instead of just last nucleotide) as ASCII.",
            cmd, false);
  */
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<size_t> kmer_length_arg("k", "kmer_length", "Length of edges (node is k-1).", true, 0, "length", cmd);
  size_t default_mem = 5 * 1024 * MB_TO_BYTES;
  TCLAP::ValueArg<size_t> mem_size_arg("m", "mem_size", "Internal memory to use (MB).", false, default_mem, "length", cmd);

  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Results will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );
  //params.ascii         = ascii_arg.getValue();
  params.input_filename  = input_filename_arg.getValue();
  params.k               = kmer_length_arg.getValue();
  params.m               = mem_size_arg.getValue() * MB_TO_BYTES;
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[])
{
  using namespace boost::adaptors;

  // Parameter extraction
  parameters_t params;
  parse_arguments(argc, argv, params);
  stxxl::internal_size_type M = params.m;
  std::string file_name = params.input_filename;
  char * base_name = basename(const_cast<char*>(file_name.c_str()));
  size_t k = params.k;

  // Open the file
  std::ifstream in(file_name, std::ifstream::binary);

  // Convert to our format: reverse for colex ordering, swap g/t encoding
  auto input = make_typed_input_range<kmer_t>(in)
             | transformed(swap_gt<kmer_t>())
             | transformed(reverse_nt<kmer_t>());

  // Create STXXL vector
  stxxl::vector<kmer_t> kmers;

  std::cerr << "Reading..." << std::endl;
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

  std::cerr << "Copying table for dummy edge detection..." << std::endl;
  // Get another copy of the table sorted by out-node (for dummy edge discovery)
  stxxl::vector<kmer_t> kmers_by_end_node(kmers);
  std::cerr << "Sorting second table..." << std::endl;
  stxxl::ksort(kmers_by_end_node.begin(), kmers_by_end_node.end(), get_key_colex_edge<kmer_t>(), M);

  // Find nodes that require incoming dummy edges
  std::cerr << "Searching for nodes requiring incoming dummy edges..." << std::endl;
  stxxl::vector<dummy_t> incoming_dummies;
  find_incoming_dummy_edges(kmers, kmers_by_end_node, k, incoming_dummies);
  std::cerr << "Added " << incoming_dummies.size() << " incoming dummy edges." << std::endl;

  // Sort dummies
  //std::cerr << "Sorting dummies..." << std::endl;
  //stxxl::sort(incoming_dummies.begin(), incoming_dummies.end(), colex_dummy_less<dummy_t>(), M);

  std::cerr << "Merging and writing..." << std::endl;

  // Make Outputter
  // TODO: Should probably do checking here when opening the file...
  string outfilename = (params.output_prefix == "")? base_name : params.output_prefix;
  ofstream ofs;
  ofs.open(outfilename + extension, ios::out | ios::binary);
  PackedEdgeOutputer out(ofs);
  #ifdef VAR_ORDER
  ofstream lcs;
  lcs.open(outfilename + extension + ".lcs", ios::out | ios::binary);
  #endif

  std::cerr << "Merging dummies and outputting..." << std::endl;

  // Merge dummies and output
  size_t prev_k = 0;
  merge_dummies(kmers, kmers_by_end_node, incoming_dummies, k,
    [&](edge_tag tag, const kmer_t & x, size_t this_k, size_t lcs_len, bool first_end_node) {
      #ifdef VAR_ORDER
      out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
      char l(lcs_len);
      lcs.write((char*)&l, 1);
      #else
      out.write(tag, x, this_k, lcs_len, first_end_node);
      #endif
      prev_k = this_k;
      
      #ifdef VERBOSE // print each kmer to stderr for testing
      if (tag == out_dummy) cout << kmer_to_string(get_start_node(x), k-1, k-1) << "$";
      else                  cout << kmer_to_string(x, k, this_k);
      cout << " " << lcs_len << " " << first_end_node << endl;
      #endif
    });


  out.close();
  #ifdef VAR_ORDER
  lcs.flush();
  lcs.close();
  #endif

  // Write out counts and stuff
  uint64_t t_k(k); // make uint64_t just to make parsing easier
  // (can read them all at once and take the last 6 values)
  ofs.write((char*)&t_k, sizeof(uint64_t));
  ofs.flush();
  ofs.close();

  return 0;
}
