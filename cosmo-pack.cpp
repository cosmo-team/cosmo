// TODO: fix up pointer and malloc to use smart arrays
// STL Headers
//#include <fstream>
#include <iostream>
//#include <algorithm>
#include <utility>

// TCLAP
#include "tclap/CmdLine.h"

// BOOST
#include <boost/range/adaptor/transformed.hpp>     // Map function to inputs
//#include <boost/range/algorithm/copy.hpp>
//#include <boost/tuple/tuple.hpp>
//#include <boost/tuple/tuple_comparison.hpp> // for uniqued over zipped iterators
//#include <boost/iterator/zip_iterator.hpp>
#define BOOST_RESULT_OF_USE_DECLTYPE // needed to support lambdas in transformed

// C STDLIB Headers
#include <cstdio>
#include <cstdlib>

#include <libgen.h>

// Custom Headers
#include "uint128_t.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"
#include "debug.h"


using namespace std;
using namespace boost::adaptors;

const static string extension = ".packed";

template <typename kmer_t, class Visitor>
void convert(kmer_t * kmers, size_t num_kmers, const uint32_t k, Visitor visit, bool swap, uint64_t *colors) {
  // Convert the nucleotide representation to allow tricks
  convert_representation(kmers, kmers, num_kmers, swap);

  // Append reverse complements
  #ifdef ADD_REVCOMPS
  size_t revcomp_factor = 2;
  transform(kmers, kmers + num_kmers, kmers + num_kmers, reverse_complement<kmer_t>(k));
  memcpy(colors + num_kmers, colors, sizeof(uint64_t) * num_kmers); // copy all of the color masks for the reverse kmers
  #else
  size_t revcomp_factor = 1;
  #endif

  // NOTE: There might be a way to do this recursively using counting (and not two tables)
  // After the sorting phase, Table A will in <colex(node), edge> order (as required for output)
  // and Table B will be in colex(row) order. Having both tables is helpful for detecting the missing dummy
  // edges with a simple O(N) merge-join algorithm.
  // NOTE: THESE SHOULD NOT BE FREED (the kmers array is freed by the caller)
  kmer_t * table_a = kmers;
  kmer_t * table_b = kmers + num_kmers * revcomp_factor; // x2 because of reverse complements
  uint64_t * colors_a = colors;
  uint64_t * colors_b = colors + num_kmers * revcomp_factor;
  // Sort by last column to do the edge-sorted part of our <colex(node), edge>-sorted table
  colex_partial_radix_sort<DNA_RADIX>(table_a, table_b, num_kmers * revcomp_factor, 0, 1,
                                      &table_a, &table_b, get_nt_functor<kmer_t>(),
                                      0, 0, 0, 0,
                                      colors_a, colors_b, &colors_a, &colors_b);
  // Sort from k to last column (not k to 1 - we need to sort by the edge column a second time to get colex(row) table)
  // Note: The output names are swapped (we want table a to be the primary table and b to be aux), because our desired
  // result is the second last iteration (<colex(node), edge>-sorted) but we still have use for the last iteration (colex(row)-sorted).
  // Hence, table_b is the output sorted from [hi-1 to lo], and table_a is the 2nd last iter sorted from (hi-1 to lo]
  colex_partial_radix_sort<DNA_RADIX>(table_a, table_b, num_kmers * revcomp_factor, 0, k,
                                      &table_b, &table_a, get_nt_functor<kmer_t>(),
                                      0, 0, 0, 0,
                                      colors_a, colors_b, &colors_b, &colors_a);

  // outgoing dummy edges are output in correct order while merging, whereas incoming dummy edges are not in the correct
  // position, but are sorted relatively, hence can be merged if collected in a previous pass
  // count dummies (to allocate space)
  size_t num_incoming_dummies = count_incoming_dummy_edges(table_a, table_b, num_kmers*revcomp_factor, k);
  TRACE("num_incoming_dummies: %zu\n", num_incoming_dummies);
  // allocate space for dummies -> we need to generate all the $-prefixed dummies, so can't just use an iterator for the
  // incoming dummies (the few that we get from the set_difference are the ones we apply $x[0:-1] to, so we need (k-1) more for each
  // (to make $...$x[0])... times two because we need to radix sort these bitches.
  #ifdef ALL_DUMMIES
  size_t all_dummies_factor = (k-1);
  size_t dummy_table_factor = 2;
  #else
  size_t all_dummies_factor = 1;
  size_t dummy_table_factor = 1;
  #endif
  // Don't have to alloc if we aren't preparing all dummies, but this option is only used for testing. Usually we want them
  kmer_t * incoming_dummies = (kmer_t*) malloc(num_incoming_dummies*all_dummies_factor*dummy_table_factor*sizeof(kmer_t));
  if (!incoming_dummies) {
    cerr << "Error allocating space for incoming dummies" << endl;
    exit(1);
  }
  // We store lengths because the prefix before the <length> symbols on the right will all be $ signs
  // this is a cheaper way than storing all symbols in 3 bits instead (although it means we need a varlen radix sort)
  uint8_t * incoming_dummy_lengths = (uint8_t*) malloc(num_incoming_dummies*all_dummies_factor*dummy_table_factor*sizeof(uint8_t));
  if (!incoming_dummy_lengths) {
    cerr << "Error allocating space for incoming dummy lengths" << endl;
    exit(1);
  }
  // extract dummies
  printf("num_kmers=%zu; revcomp_factor=%zu; incoming_dummmies=%zu\n", num_kmers, revcomp_factor, num_incoming_dummies);
  find_incoming_dummy_edges(table_a, table_b, num_kmers*revcomp_factor, k, incoming_dummies);
  // add extra dummies
  #ifdef ALL_DUMMIES
  prepare_incoming_dummy_edges(incoming_dummies, incoming_dummy_lengths, num_incoming_dummies, k-1);
  #else
  // Just set the lengths for merging
  memset(incoming_dummy_lengths, k-1, num_incoming_dummies);
  #endif

  kmer_t * dummies_a = incoming_dummies;
  uint8_t * lengths_a = incoming_dummy_lengths;
  // sort dummies (varlen radix)
  // dont need to sort if not adding the extras, since already sorted
  #ifdef ALL_DUMMIES
  kmer_t * dummies_b = incoming_dummies + num_incoming_dummies * (k-1);
  uint8_t * lengths_b = incoming_dummy_lengths + num_incoming_dummies * (k-1);
  colex_partial_radix_sort<DNA_RADIX>(dummies_a, dummies_b, num_incoming_dummies*(k-1), 0, 1,
                                      &dummies_a, &dummies_b, get_nt_functor<kmer_t>(),
                                      lengths_a, lengths_b, &lengths_a, &lengths_b);
  // Don't need the last iteration (i.e. dont need to go to 0) since we arent doing a set difference like above
  colex_partial_radix_sort<DNA_RADIX>(dummies_a, dummies_b, num_incoming_dummies*(k-1), 1, k-1,
                                      &dummies_a, &dummies_b, get_nt_functor<kmer_t>(),
                                      lengths_a, lengths_b, &lengths_a, &lengths_b);
  #endif
  merge_dummies(table_a, table_b, num_kmers*revcomp_factor, k,
                dummies_a, num_incoming_dummies*all_dummies_factor,
                lengths_a,
                // edge_tag needed to distinguish between dummy out edge or not...
                [=](edge_tag tag, const kmer_t & x, const uint32_t x_k, size_t first_start_node, bool first_end_node) {
                  // TODO: this should be factored into a class that prints full kmers in ascii
                  // then add a --full option
                  #ifdef VERBOSE // print each kmer to stderr for testing
                  if (tag == out_dummy)
                    cerr << kmer_to_string(get_start_node(x), k-1, k-1) << "$";
                  else
                    cerr << kmer_to_string(x, k, x_k);
                  cerr << " " << first_start_node << " " << first_end_node << endl;
                  #endif
                  visit(tag, x, x_k, first_start_node, first_end_node);
                });
  // TODO: impl external-merge (for large input. Read in chunk, sort, write temp, ext merge + add dummies to temp, 3-way-merge)
  // TODO: use SSE instructions or CUDA if available (very far horizon)

  free(incoming_dummies);
  free(incoming_dummy_lengths);
}

typedef struct p
{
    //bool ascii = false;
    std::string input_filename = "";
    std::string output_prefix = "";
    bool cortex = false;
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
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Results will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  TCLAP::SwitchArg cortex_arg("c", "cortex",
            "Input file is cortex binary",
            cmd, false);
  cmd.parse( argc, argv );
  //params.ascii         = ascii_arg.getValue();
  params.input_filename  = input_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  params.cortex          = cortex_arg.getValue();
}

int main(int argc, char * argv[]) {
  parameters_t params;
  parse_arguments(argc, argv, params);

  const char * file_name = params.input_filename.c_str();

  // Open File
  int handle = -1;
  if ( (handle = open(file_name, O_RDONLY)) == -1 ) {
    fprintf(stderr, "ERROR: Can't open file: %s\n", file_name);
    exit(EXIT_FAILURE);
  }

  // The parameter should be const... On my computer the parameter
  // isn't const though, yet it doesn't modify the string...
  // This is still done AFTER loading the file just in case
  char * base_name = basename(const_cast<char*>(file_name));

  // Read Header
  uint32_t kmer_num_bits = 0;
  uint32_t k = 0;
  if (params.cortex) {
    if ( !cortex_read_header(handle, &kmer_num_bits, &k) ) {
      fprintf(stderr, "ERROR: Error reading cortex_file %s\n", file_name);
      exit(EXIT_FAILURE);
    }
  }
  else {
    if ( !dsk_read_header(handle, &kmer_num_bits, &k) ) {
      fprintf(stderr, "ERROR: Error reading file %s\n", file_name);
      exit(EXIT_FAILURE);
    }
  }
  uint32_t kmer_num_blocks = (kmer_num_bits / 8) / sizeof(uint64_t);
  TRACE(">> READING DSK FILE\n");
  TRACE("kmer_num_bits, k = %d, %d\n", kmer_num_bits, k);
  TRACE("kmer_num_blocks = %d\n", kmer_num_blocks);

  if (kmer_num_bits > MAX_BITS_PER_KMER) {
    fprintf(stderr, "ERROR: Kmers larger than %zu bits are not currently supported."
                    " %s uses %d bits per kmer (possibly corrupt?).\n",
            MAX_BITS_PER_KMER, file_name, kmer_num_bits);
    exit(EXIT_FAILURE);
  }

  // Read how many items there are (for allocation purposes)
  size_t num_kmers = 0;
  uint32_t num_colors = 0;
  if (params.cortex) {
    printf("Num records\n");
    if ( cortex_num_records(handle, kmer_num_bits, &num_kmers, &num_colors) == -1) {
      fprintf(stderr, "Error seeking cortex file %s\n", file_name);
      exit(EXIT_FAILURE);
    }
    printf("Got num record %zu \n", num_kmers);
  }
  else {
    if ( dsk_num_records(handle, kmer_num_bits, &num_kmers) == -1) {
      fprintf(stderr, "Error seeking file %s\n", file_name);
      exit(EXIT_FAILURE);
    }
  }
  if (num_kmers == 0) {
    fprintf(stderr, "ERROR: File %s has no kmers (possibly corrupt?).\n", file_name);
    exit(EXIT_FAILURE);
  }
  TRACE("num_kmers = %zu\n", num_kmers);

  // ALLOCATE SPACE FOR KMERS (done in one malloc call)
  // x 4 because we need to add reverse complements, and then we have two copies of the table
  #ifdef ADD_REVCOMPS
  size_t revcomp_factor = 2;
  #else
  size_t revcomp_factor = 1;
  #endif
  uint64_t * kmer_blocks = (uint64_t*)malloc(num_kmers * 2 * revcomp_factor * sizeof(uint64_t) * kmer_num_blocks);
  if (!kmer_blocks) {
    cerr << "Error allocating space for kmers" << endl;
    exit(1);
  }

  uint64_t * kmer_colors = (uint64_t*)malloc(num_kmers * 2 * revcomp_factor * sizeof(uint64_t));
  if (!kmer_colors) {
    cerr << "Error allocating space for kmer colors" << endl;
    exit(1);
  }


  // READ KMERS FROM DISK INTO ARRAY
  size_t num_records_read;
  if (params.cortex) {
    printf("Reading kmers\n");
    num_records_read = cortex_read_kmers(handle, kmer_num_bits, num_colors, k, kmer_blocks, kmer_colors);
    printf("num_kmers = %zu and num_records_read=%zu\n", num_kmers, num_records_read);
    num_kmers = num_records_read;
  }
  else {
    num_records_read = dsk_read_kmers(handle, kmer_num_bits, kmer_blocks);
  }
  close(handle);
  if (num_records_read == 0) {
    fprintf(stderr, "Error reading file %s\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  TRACE("num_records_read = %zu\n", num_records_read);
  assert (num_records_read == num_kmers);

  //auto ascii_output = std::ostream_iterator<string>(std::cout, "\n");

  string outfilename = (params.output_prefix == "")? base_name : params.output_prefix;
  ofstream ofs;
  #ifdef VAR_ORDER
  ofstream lcs;
  #endif
  #if 1
  static const size_t BUFFER_LEN = 2*1024*1024;
  char buffer_a[BUFFER_LEN];
  ofs.rdbuf()->pubsetbuf(buffer_a, BUFFER_LEN);
  #ifdef VAR_ORDER
  char buffer_b[BUFFER_LEN];
  lcs.rdbuf()->pubsetbuf(buffer_b, BUFFER_LEN);
  #endif
  #endif
  // TODO: Should probably do checking here when opening the file...
  ofs.open(outfilename + extension, ios::out | ios::binary);
  #ifdef VAR_ORDER
  lcs.open(outfilename + extension + ".lcs", ios::out | ios::binary);
  #endif
  PackedEdgeOutputer out(ofs);
  ofstream cfs;
  cfs.open(outfilename + ".colors", ios::out | ios::binary);

  if (kmer_num_bits == 64) {
    typedef uint64_t kmer_t;
    size_t prev_k = 0; // for input, k is always >= 1
    size_t index = 0;
    uint64_t *colors = kmer_colors + num_kmers * revcomp_factor;
      convert(kmer_blocks, num_kmers, k,
        [&](edge_tag tag, const kmer_t & x, const uint32_t this_k, size_t lcs_len, bool first_end_node) {
          #ifdef VAR_ORDER
          out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
          char l(lcs_len);
          lcs.write((char*)&l, 1);
          #else
          out.write(tag, x, this_k, lcs_len, first_end_node);
          #endif
	  if (tag == standard) {
            // cerr << kmer_to_string(x, k, this_k) << "c" << colors[index] << "\n";
	    cfs.write((char *)&colors[index++], sizeof(uint64_t));
	  }
	  else {
	    uint64_t zero = 0;
	    cfs.write((char *)&zero, sizeof(uint64_t));
	  }
          prev_k = this_k;
	      }, !params.cortex, kmer_colors);
  }
  else if (kmer_num_bits == 128) {
    typedef uint128_t kmer_t;
    size_t prev_k = 0;
    kmer_t * kmer_blocks_128 = (kmer_t*)kmer_blocks;
    convert(kmer_blocks_128, num_kmers, k,
        [&](edge_tag tag, const kmer_t & x, const uint32_t this_k, size_t lcs_len, bool first_end_node) {
          #ifdef VAR_ORDER
          out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
          char l(lcs_len);
          lcs.write((char*)&l, 1);
          #else
          out.write(tag, x, this_k, lcs_len, first_end_node);
          #endif
          prev_k = this_k;
	    }, !params.cortex, kmer_colors);
  }

  out.close();
  #ifdef VAR_ORDER
  lcs.flush();
  lcs.close();
  #endif
  uint64_t t_k(k); // make uint64_t just to make parsing easier
  // (can read them all at once and take the last 6 values)
  ofs.write((char*)&t_k, sizeof(uint64_t));
  ofs.flush();
  ofs.close();

  cfs.close();

  free(kmer_blocks);
  free(kmer_colors);
  return 0;
}
