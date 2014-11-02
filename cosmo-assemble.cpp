#include <iostream>
#include <fstream>

#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
//#include <boost_iterator/zip_iterator.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"
#include "debruijn_hypergraph.hpp"
#include "algorithm.hpp"
#include "wt_algorithm.hpp"

using namespace std;
using namespace sdsl;

string graph_extension = ".dbg";
string contig_extension = ".fasta";

struct parameters_t {
  std::string input_filename = "";
  std::string output_prefix = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".dbg file (output from cosmo-build).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Contigs will be written to [" + output_short_form + "]" + contig_extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  // -d flag for decompression to original kmer biz
  params.input_filename  = input_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  // The parameter should be const... On my computer the parameter
  // isn't const though, yet it doesn't modify the string...
  // This is still done AFTER loading the file just in case
  char * base_name = basename(const_cast<char*>(p.input_filename.c_str()));
  string outfilename = ((p.output_prefix == "")? base_name : p.output_prefix);

  // TO LOAD:
  debruijn_graph<> g;
  load_from_file(g, p.input_filename);

  cerr << "k             : " << g.k << endl;
  cerr << "num_nodes()   : " << g.num_nodes() << endl;
  cerr << "num_edges()   : " << g.num_edges() << endl;
  cerr << "W size        : " << size_in_mega_bytes(g.m_edges) << " MB" << endl;
  cerr << "L size        : " << size_in_mega_bytes(g.m_node_flags) << " MB" << endl;
  cerr << "Total size    : " << size_in_mega_bytes(g) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(g) << " Bits" << endl;

  #ifdef VAR_ORDER
  wt_int<rrr_vector<63>> lcs;
  load_from_file(lcs, p.input_filename + ".lcs.wt");

  cerr << "LCS size      : " << size_in_mega_bytes(lcs) << " MB" << endl;
  cerr << "LCS bits/edge : " << bits_per_element(lcs) << " Bits" << endl;

  typedef debruijn_hypergraph<> dbh;
  dbh h(g, lcs);

  typedef dbh::node_type node_type;
  //node_type v(12,215855); // ...aaa*
  //node_type v(3422777,3422778); // ...ttattccgtagc->[t,g]. Other than that, totally different nodes.
  node_type v(3422774,3422778, 11); // ....tattccgtagc->[t3,g2]
  //node_type v(3422771,3422796); // ......ttccgtagc->[acgt]
  node_type u = h.maxlen(v);
  cout << "maxlen: " << get<0>(u) << ", " << get<1>(u) << endl;
  auto y = h.maxlen(v,2);
  if (!y) cout << "NONE" << endl;
  else cout << "maxlen: " << get<0>(*y) << ", " << get<1>(*y) << endl;
  //size_t plte = prev_lte(lcs, 1, 2);
  //cout << plte << endl;
  //cout << "bit: " << get_bit_at_level(lcs, 1, 15) << endl;


  //construct_im(wt, int_vector<>({2, 4, 1, 1, 0, 1, 1, 20}));
  //wt_int<rrr_vector<63>> wt; //0 1  2  3  4  5  6  7  8
  //construct_im(wt, int_vector<>({0, 6, 6, 2, 5, 0}));
  //auto plte = prev_lte(wt, 6, 100);
  //auto nlte = next_lte(wt, 5, 4); // 5
  //cout << plte << endl;
  //cout << nlte << endl;
  //auto v_s = h.shorter(v, 3);
  //cout << "shorter: " get<0>(v_s) << ", " << get<1>(v_s) << endl;
  //auto v_l = h.longer(v, 3);
  auto r = h.longer(v,12);

  for (auto x : r) {
    cout << get<0>(x) << ", " << get<1>(x) << endl;
  }

  cout << "lastchar: " << "$acgt"[(int)h.lastchar(v)] << endl;

  // node_type v(3422774,3422778, g.k-1); // ....tattccgtagc->[t3,g2]
  // 01234
  // $acgt
  auto q = h.outgoing(v, 4);
  if (!q) cout << "NONE" << endl;
  else cout << get<0>(*q) << ", " << get<1>(*q) << endl;

  cout << "backward:" <<endl;
  auto prevs = h.backward(v);
  for (auto prev: prevs) {
    cout << get<0>(prev) << ", " << get<1>(prev) << endl;
  }

  #endif

  /*
  int num_queries = 10;
  typedef boost::mt19937 rng_type;
  rng_type rng(time(0));
  boost::uniform_int<size_t> node_distribution(0,g.num_nodes()); // make go up to size of graph
  boost::uniform_int<size_t> k_distribution(0, g.k); // make go up to size of graph
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_node(rng, node_distribution);
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_k(rng, k_distribution);
  cout << random_node() << ", " << random_k() << endl;
  vector<size_t> q_nodes(boost::make_function_input_iterator(random_node,0),
                         boost::make_function_input_iterator(random_node,num_queries));
  //vector<size_t>    q_ks(boost::make_function_input_iterator(random_node,0),
  //                       boost::make_function_input_iterator(random_node,10));
  for (auto x: q_nodes) {
    cout << x << endl;
  }
  */
}

