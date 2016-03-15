//#include <iterator>
//#include <typeinfo>
#include <fstream>
//#include <bitset>
//#include <set>

#include <stxxl.h>

/*
#include <stxxl/bits/containers/sorter.h>
#include <boost/range/adaptors.hpp>
#include <boost/range/adaptor/uniqued.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/istream_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
*/
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>

// TCLAP
#include "tclap/CmdLine.h"

#include "config.hpp"
#include "debug.hpp"
#include "utility.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"

struct parameters_t {
  std::string input_filename = "";
  std::string output_prefix = "";
  size_t k = 0;
  size_t m = 0;
  bool swap = true;
  bool variable_order = false;
  bool shift_dummies  = true;
};

parameters_t parse_arguments(int argc, char **argv) {
  parameters_t params;
  TCLAP::CmdLine cmd(banner, ' ', version);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
    "Input file (DSK or Cortex output).",
    true, "", "input_file", cmd);
  TCLAP::ValueArg<size_t> kmer_length_arg("k", "kmer_length",
    "Length of edges (node is k-1).",
    false, 0, "length", cmd);
  TCLAP::ValueArg<size_t> mem_size_arg("m", "mem_size",
    "Internal memory to use (MB).",
    false, default_mem_size, "mem_size", cmd);
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
    "Output prefix. Default prefix: basename(input_file).",
    false, "", "output_prefix", cmd);
  TCLAP::SwitchArg varord_arg("v", "variable_order", "Output .lcs file for variable order support.", cmd, false);
  TCLAP::SwitchArg swap_arg("g", "swap_gt", "Swap g and t representation in kmers (Use if DSK input).", cmd, false);
  TCLAP::SwitchArg shift_arg("d", "no_shifts", "Don't shift all incoming dummies (faster, but loses some information).", cmd, false);
  // TODO: ASCII outputter: -a for last symbol, -aa for whole kmer
  // TODO: XORed forced input format switch
  // TODO: add DSK count parser (-c)
  // TODO: add option to save full sorted kmers (e.g. for iterative dbg)
  // TODO: add option for no reverse complements (e.g. if using bcalm input)
  // TODO: add option for mutliple files being merged for colour
  cmd.parse( argc, argv );
  params.input_filename  = input_filename_arg.getValue();
  params.k               = kmer_length_arg.getValue();
  params.m               = mem_size_arg.getValue() * mb_to_bytes;
  params.variable_order  = varord_arg.getValue();
  params.swap            = swap_arg.getValue();
  params.shift_dummies   = !shift_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  return params;
}

int main(int argc, char* argv[]) {
  using namespace boost::adaptors;

  // Parameter extraction
  auto params = parse_arguments(argc, argv);
  stxxl::internal_size_type M = params.m;
  std::string file_name = params.input_filename;
  //string base_name = basename(file_name);
  size_t k = params.k;

  // Set logging level
  // TODO: make verbosity parameter
  boost::log::core::get()->set_filter (
    boost::log::trivial::severity != boost::log::trivial::debug
  );

  // Check format
  // TODO: add fastq/fasta input, and read from stdin
  bool ctx_input = (extension(file_name) == ".ctx");
  COSMO_LOG(info) << "Loading file: " << file_name;
  if (ctx_input) {
    COSMO_LOG(info) << "Input format: Cortex file";
  }
  else {
    COSMO_LOG(info) << "Input format: Raw file";
    if (k == 0) {
      COSMO_LOG(error) << "When using raw format, please provide a value for -k flag.";
      exit(1);
    }
  }

  // Check k value supported
  if (k > max_k) {
    COSMO_LOG(error) << "This version only supports k <= " << max_k << ". Try recompiling.";
    exit(1);
  }

  if (ctx_input) {
    //ifstream in_file(file_name, ios::binary | ios::in);
    COSMO_LOG(error) << "Cortex input not yet supported.";
    exit(1);
  }
  else {
    // stxxl::syscall_file out_file(file_name + ".boss", stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
    typedef kmer_sorter<kmer_t> kmer_sorter_t;
    typedef kmer_sorter_t::record_vector_t record_vector_t;
    typedef kmer_sorter_t::kmer_vector_t kmer_vector_t;
    typedef kmer_sorter_t::output_t output_t;

    stxxl::syscall_file in_file(file_name, stxxl::file::DIRECT | stxxl::file::RDONLY);
    record_vector_t in_vec(&in_file);
    record_vector_t::bufreader_type reader(in_vec);
    kmer_sorter_t sort_input;

    size_t idx = 0;
    sort_input.sort(reader, params, [&](output_t result){
      //kmer_t x = get<0>(result.record);
      // TODO: test with a colour payload
      //auto payload = result.record.get_tail();

      //char w = (result.tag == out_dummy)?'$':(DNA_ALPHA "ACGT")[get_edge_label(x) | (result.is_first_suffix<<2)];
      //cout << w << endl;

      /*
      if (result.tag == out_dummy) cout << kmer_to_string(get_start_node(x<<2), k-1, k-1) << "$";
      else                         cout << kmer_to_string(x, k, k);
      cout << endl;
      */

      /*
      cout << result.is_first_prefix
           << " "
           << result.is_first_suffix
           << endl;
      */
    });
  }

  /*
  // Make Outputter
  // TODO: Should probably do checking here when opening the file...
  string outfilename = (params.output_prefix == "")? base_name : params.output_prefix;
  ofstream ofs;
  ofs.open(outfilename + packed_ext, ios::out | ios::binary);
  PackedEdgeOutputer out(ofs);
  #ifdef VAR_ORDER
  ofstream lcs;
  lcs.open(outfilename + packed_ext + lcs_ext, ios::out | ios::binary);
  #endif
  ofstream cols;
  cols.open(outfilename + packed_ext + color_ext, ios::out | ios::binary);

  merge_dummies(xx2|transformed(capture_color), xb2, d, k,
    [&](edge_tag tag, const kmer_t & x, size_t this_k, size_t lcs_len, bool first_end_node) {
      #ifdef VAR_ORDER
      out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
      char l(lcs_len);
      lcs.write((char*)&l, 1);
      #else
      out.write(tag, x, this_k, lcs_len, first_end_node);
      #endif
      cols.write((char*)&color, sizeof(color_t));
      prev_k = this_k;

      #ifdef VERBOSE // print each kmer to stderr for testing
      if (tag == out_dummy) cout << kmer_to_string(get_start_node(x), k-1, k-1) << "$";
      else                  cout << kmer_to_string(x, k, this_k);
      cout << " " << (lcs_len != k-1) << " " << lcs_len << " " << first_end_node;
      // print color
      bitset<max_colors> bs(0); // dummy edges -> 0
      if (tag == standard) bs = color;
      else bs;
      cout << " ";
      for (int i = 0; i<number_of_colours;i++) cout << bs[i];
      cout << endl;
      #endif
      if (tag == out_dummy) num_out_dummies++;
      else if (tag == in_dummy) num_in_dummies++;
    });
  COSMO_LOG(info) << "Number of unique incoming dummies: " << num_in_dummies;
  COSMO_LOG(info) << "Number of outgoing dummies: " << num_out_dummies;

  out.close();
  #ifdef VAR_ORDER
  lcs.flush();
  lcs.close();
  #endif
  cols.flush();
  cols.close();


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
  */

  COSMO_LOG(trace) << "Done!";

  return 0;
}
