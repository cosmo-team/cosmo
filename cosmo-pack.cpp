#include <iterator>
#include <typeinfo>
#include <fstream>
#include <libgen.h>
#include <stxxl.h>
#include <bitset>
#include <set>


#include <stxxl/bits/containers/sorter.h>
#include <stxxl/bits/parallel.h>

#include <boost/range/adaptors.hpp>
#include <boost/range/adaptor/uniqued.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/istream_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

// TCLAP
#include "tclap/CmdLine.h"

#include "config.hpp"
#include "debug.hpp"
#include "utility.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"

// TODO: make generic way to do this, and make non-key fields be stored somewhere else
// (use zip iterators) since we ignore them for most of it...
// so we want to prefetch as many keys as possible
typedef std::tuple<kmer_t, color_t>  record_t;

typedef node_less<kmer_t> node_comparator_type;
typedef record_less<record_t, kmer_t> record_comparator_type;
typedef kmer_less<kmer_t> edge_comparator_type;
typedef colex_dummy_less<dummy_t, kmer_t> dummy_comparator_type;

typedef stxxl::vector<record_t, 1, stxxl::lru_pager<8>, block_size> record_vector_type;
typedef stxxl::vector<kmer_t, 1, stxxl::lru_pager<8>, block_size> kmer_vector_type;

namespace std {
  template<>
  struct iterator_traits<record_vector_type::bufreader_type::bufreader_iterator> {
    typedef input_iterator_tag iterator_category;
    typedef record_t value_type;
    typedef size_t   difference_type;
    typedef record_t& reference;
  };
  template<>
  struct iterator_traits<kmer_vector_type::bufreader_type::bufreader_iterator> {
    typedef input_iterator_tag iterator_category;
    typedef kmer_t value_type;
    typedef size_t   difference_type;
    typedef kmer_t& reference;
  };
}

typedef stxxl::sorter<record_t, record_comparator_type, block_size> record_sorter_type;
typedef stxxl::sorter<kmer_t,   edge_comparator_type, block_size> edge_sorter_type;
typedef stxxl::sorter<dummy_t,  dummy_comparator_type, block_size> dummy_sorter_type;

struct parameters_t
{
  std::string input_filename = "";
  std::string output_prefix = "";
  size_t k = 0;
  size_t m = 0;
};

void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd(banner, ' ', version);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<size_t> kmer_length_arg("k", "kmer_length", "Length of edges (node is k-1).", false, 0, "length", cmd); // NOTE: made this optional for the cortex input
  TCLAP::ValueArg<size_t> mem_size_arg("m", "mem_size", "Internal memory to use (MB).", false, default_mem_size, "mem size", cmd);

  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Default prefix: basename(input_file).",
            false, "", output_short_form, cmd);
  cmd.parse( argc, argv );
  params.input_filename  = input_filename_arg.getValue();
  params.k               = kmer_length_arg.getValue();
  params.m               = mem_size_arg.getValue() * mb_to_bytes;
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

  // Check format

  // TODO: ASCII outputter (first and last symbols, option for whole kmer)
  // TODO: stdin input
  boost::log::core::get()->set_filter (
    boost::log::trivial::severity != boost::log::trivial::debug
  );

  // Check k value supported
  if (k > K_LEN) {
    COSMO_LOG(error) << "This version only supports k <= " << K_LEN << ". Try recompiling.";
    exit(1);
  }

  // TODO: Test sort input size (with sorter ctor, overload sort func)
  // block size vs creator size?

  // Open the file
  //stxxl::syscall_file in_file(file_name, stxxl::file::DIRECT | stxxl::file::RDONLY);
  ifstream in_file(file_name, ios::binary | ios::in);
  //stxxl::syscall_file out_file(file_name + ".boss",
  //stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);

  // Consume header
  
  in_file.ignore(6+sizeof(int));
  int input_kmer_size;
  in_file.read((char*)&input_kmer_size,sizeof(int));
  int number_of_bitfields;
  int number_of_colours;
  in_file.read((char*)&number_of_bitfields,sizeof(int));
  in_file.read((char*)&number_of_colours,sizeof(int));
  k = input_kmer_size + 1;
 
  COSMO_LOG(info) << "node length (k)     : " << input_kmer_size;
  COSMO_LOG(info) << "number of bitfields : " << number_of_bitfields;
  COSMO_LOG(info) << "number of colors    : " << number_of_colours;

  in_file.ignore(number_of_colours*(sizeof(int)+sizeof(long long)));
  for (int i = 0; i < number_of_colours; i++) {
    int sample_id_lens;
    in_file.read((char*)&sample_id_lens, sizeof(int));
    in_file.ignore(sample_id_lens);
  }
  in_file.ignore(number_of_colours*sizeof(long double));
  for (int i = 0; i < number_of_colours; i++) {
    int len_name_of_graph;
    in_file.ignore(4 + 2*sizeof(int));
    in_file.read((char*)&len_name_of_graph, sizeof(int));
    in_file.ignore(len_name_of_graph);
  }
  in_file.ignore(6);
  

  // Header consumed. find how many records there are
  size_t start = in_file.tellg();
  in_file.seekg (0, in_file.end);
  size_t length = (size_t)in_file.tellg() - start;
  size_t size = length/(8 + number_of_colours * sizeof(int) + number_of_colours);
  in_file.seekg (start, in_file.beg);

  //vector_type in_vec;//(&in_file); 
  //vector_type kmers(&out_file);
  record_vector_type kmers;
  kmer_vector_type kmers_b;

  //node_sorter_type node_sorter(node_comparator_type(), M/2);
  typedef set<dummy_t> set_t;
  set_t temp_set;
  size_t max_temp_dummies = (M/2)/sizeof(dummy_t);
  record_sorter_type record_sorter(record_comparator_type(), M/2);
  edge_sorter_type edge_sorter(edge_comparator_type(), M/2);

  // Make conversion functors
  //auto swap  = swap_gt<kmer_t>();
  auto revnt = reverse_nt<kmer_t>();
  auto rc    = reverse_complement<kmer_t>(k);

  // Convert to our format: reverse for colex ordering, swap g/t encoding (DSK)
  COSMO_LOG(trace) << "Creating runs...";
  // TODO: Parallelise? #pragma omp parallel for
  // TODO: convert to asynch IO
  for (size_t i = 0; i < size; i++) {
    kmer_t in_kmer;
    int * covg = new int[number_of_colours];
    //vector<bitset<8>>  edges(number_of_colours,bitset<8>(0));
    // TODO: change these to vectors when I work out wtf is going on with the reading bug
    char * edges = new char[number_of_colours];
    vector<bitset<max_colors>> colors_per_edge(5); // 5th is node only
    
    in_file.read((char*)&in_kmer,sizeof(kmer_t));
    in_file.read((char*)covg,number_of_colours*sizeof(int));
    in_file.read((char*)edges,number_of_colours); 

    // Invert color matrix
    for (int j = 0; j < number_of_colours; j++) {
      if (covg[j] == 0) {
        continue; // skip colours with no coverage
      }
      // coverage for kmer but no edges...
      bool node_only = bitset<8>(edges[j]).none();
      // TODO: Parallelise and optimise this
      if (node_only) {
        colors_per_edge[4][j] = 1;
        //cout << "NODE ONLY (" << j << "): " 
        //     << kmer_to_string(revnt(in_kmer), input_kmer_size);
        //cout << "$" << endl;
      }
      else for (int c : {0,1,2,3}) {
        colors_per_edge[c][j] = bitset<8>(edges[j])[c];
      }
    }

    // add edges + colors
    for (int c : {0,1,2,3}) {
      if (colors_per_edge[c].any()) {
        kmer_t x = revnt((in_kmer << 2)|c);
        //cout << kmer_to_string(edge,input_kmer_size+1) << " ";
        //cout << colors_per_edge[c] << endl;
        auto y = rc(x);
        record_sorter.push(record_t(x,colors_per_edge[c].to_ulong()));
        record_sorter.push(record_t(y,colors_per_edge[c].to_ulong()));
        edge_sorter.push(x);
        edge_sorter.push(y);
      }
      // TODO: what should we do with the node_only colors? especially if no edges.
      // by BOSS dBG definition, we can ignore them, so will do that for now.
      // but maybe later can add outgoing dummies with colours using the 5th row.
    }
    delete[] covg;
    delete[] edges;
  }
  COSMO_LOG(info) << "Added " << record_sorter.size()/2 << " edges, not including revcomps.";

  /*
  for (kmer_t record : vector_type::bufreader_type(in_vec)) {
    //auto x = revnt(swap(record));
    auto x = revnt(record);
    auto y = rc(x);
    node_sorter.push(x);
    node_sorter.push(y);
    edge_sorter.push(x);
    edge_sorter.push(y);
  }
  COSMO_LOG(info) << "Added " << node_sorter.size()/2 << " kmers, and their reverse complements.";
  */

  // TODO: replace internal sort with nvidia radix sort (cub/thrust)
  // TODO: test simulated B table vs second sort
  COSMO_LOG(trace) << "Merging runs...";
  record_sorter.sort();
  edge_sorter.sort();

  // TODO: make buffered reader around sorted stream instead
  // Or keep sorting two tables and stream to range for dummy edge finding (then rewind)
  
  COSMO_LOG(trace) << "Writing to temporary storage...";
  kmers.resize(record_sorter.size());
  edge_sorter.set_merger_memory_to_use (0);
  record_sorter.set_merger_memory_to_use (M);
  stxxl::stream::materialize(record_sorter, kmers.begin(), kmers.end());
  record_sorter.finish_clear();
  edge_sorter.set_merger_memory_to_use (M);
  kmers_b.resize(edge_sorter.size());
  stxxl::stream::materialize(edge_sorter, kmers_b.begin(), kmers_b.end());
  edge_sorter.finish_clear();

  // TODO: time materialize + iterate over twice vs iterate over sorter twice
  // (One table is fine)

  // TODO: read about STXXL_PARALLEL_MULTIWAY_MERGE and other defs
  // Find nodes that require incoming dummy edges
  /*
  vector_type::bufreader_type br(kmers);
  auto f= std::function<kmer_t(void)>([&]()->kmer_t{
    kmer_t temp = *br;
    ++br;
    return temp;
  });
  typedef boost::function_input_iterator<decltype(f), size_t> fiit;
  fiit foo_0 = fiit(f, (size_t)0);
  fiit foo_1 = fiit(f, kmers.size());
  auto foo = boost::make_iterator_range(foo_0, foo_1);
  */
  // TODO: Write sort stream (or buffered reader) to iterator adapter
  // TODO: Make iterator adapter for stream/bufreader (*, ++)
  // Might not have worked before because of RAII?
  // TODO: Make range that takes two ranges and calculates set difference
  
  stxxl::vector<dummy_t> incoming_dummies;
  #define CI_RANGE(x) (boost::make_iterator_range((x).cbegin(),(x).cend()))
  auto a = CI_RANGE(kmers);
  auto b = CI_RANGE(kmers_b);
  // TODO: function output to buffered writer or sorter, for dummy edges
  // or a b-tree if generating all shifts
  COSMO_LOG(trace) << "Searching for nodes requiring incoming dummy edges...";

  //dummy_sorter_type dummy_sorter(dummy_comparator_type(), M/2);
  // TODO: redo this so we dont need the dummies (follow up with Travis)
 
  std::function<kmer_t(record_t)> record_key([](record_t x) -> kmer_t {
    return get<0>(x);
  });

  size_t num_dummies = 0; 
  typedef decltype(kmers)::bufreader_type bra_t;
  typedef decltype(kmers_b)::bufreader_type brb_t;
  bra_t bra(kmers);
  brb_t brb(kmers_b);
  auto xx = boost::make_iterator_range(bra.begin(), bra.end());
  auto xb = boost::make_iterator_range(brb.begin(), brb.end());
  find_incoming_dummy_nodes<kmer_t>(xx | transformed(record_key), xb, k,
    [&](std::pair<size_t, kmer_t> x) {
      num_dummies++;
      //cout << get<0>(x) << ": " << kmer_to_string(get<1>(x),k) << endl;
      //dummy_sorter.push(dummy_t(x, k-1));
    }
  );
  COSMO_LOG(info)  << "Added " << num_dummies << " incoming dummy edges.";
  //COSMO_LOG(trace) << "Sorting dummies...";
  //dummy_sorter.sort();

  //incoming_dummies.resize(dummy_sorter.size());
  //stxxl::stream::materialize(dummy_sorter, incoming_dummies.begin(), incoming_dummies.end());
  //dummy_sorter.finish_clear();

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

  COSMO_LOG(trace) << "Merging dummies and outputting...";
  // Merge dummies and output
  size_t prev_k = 0;
  color_t color;
  std::function<kmer_t(record_t)> capture_color([&](record_t x){
    color = get<1>(x);
    return get<0>(x);
  });
  size_t num_in_dummies = 0;
  size_t num_out_dummies = 0;
  auto d = CI_RANGE(incoming_dummies);
  bra.rewind();
  brb.rewind();
  bra_t bra2(kmers);
  brb_t brb2(kmers_b);
  auto xx2 = boost::make_iterator_range(bra2.begin(), bra2.end());
  auto xb2 = boost::make_iterator_range(brb2.begin(), brb2.end());
  // TODO: refactor
  // TODO: Use STXXL ASIO for input/output
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
  COSMO_LOG(trace) << "Done!";

  return 0;
}
