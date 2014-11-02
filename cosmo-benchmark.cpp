#include <iostream>
#include <fstream>
#include <algorithm>
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
  typedef dbh::node_type node_type;
  dbh h(g, lcs);
  #endif


  int num_queries = 1e5;
  size_t min_k = 0;
  typedef boost::mt19937 rng_type;
  rng_type rng(time(0));
  boost::uniform_int<size_t> node_distribution(0,g.num_nodes()-1); // make go up to size of graph
  boost::uniform_int<size_t> k_distribution(min_k, g.k-1); // make go up to size of graph
  boost::uniform_int<size_t> symbol_distribution(1, 4); // make go up to size of graph
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_node(rng, node_distribution);
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_k(rng, k_distribution);
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_symbol(rng, symbol_distribution);
  //auto random_higher_k = [&](size_t low)  { return boost::uniform_int<size_t>(low+1,g.k-1)(rng); }; // make go up to size of graph
  //auto random_lower_k  = [&](size_t high) { return boost::uniform_int<size_t>(min_k,high-1)(rng); }; // make go up to size of graph

  vector<size_t> query_nodes(boost::make_function_input_iterator(random_node,0),
                             boost::make_function_input_iterator(random_node,num_queries));
  vector<size_t> query_syms(boost::make_function_input_iterator(random_symbol,0),
                            boost::make_function_input_iterator(random_symbol,num_queries));
  #ifdef VAR_ORDER
  // Convert to variable order nodes
  vector<size_t> query_ks(boost::make_function_input_iterator(random_k,0),
                          boost::make_function_input_iterator(random_k,num_queries));

  vector<node_type> query_varnodes;
  //transform(query_nodes.begin(), query_nodes.end(), query_varnodes.begin(),[&](size_t v){ return h.get_node(v); });
  for (auto u:query_nodes) {
    query_varnodes.push_back(h.get_node(u));
  }

  // set shorter for each one that has a shorter k
  for (int i=0;i<num_queries;i++) {
    if (query_ks[i] == g.k-1) continue;
    query_varnodes[i] = h.shorter(query_varnodes[i], query_ks[i]);
  }
  #endif

  #ifndef VAR_ORDER // standard dbg
  // backward
  for (auto v : query_nodes) { g.all_preds(v); }

  // forward
  for (size_t i=0;i<(size_t)num_queries;i++) {
    size_t v = query_nodes[i];
    auto   x = query_syms[i];
    g.outgoing(v, x);
  }

  // last_char
  for (auto v : query_nodes) { g.lastchar(v); }
  #else

  // backward
  for (auto v : query_varnodes) {
    h.backward(v);
  }

  // forward
  for (size_t i=0;i<(size_t)num_queries;i++) {
    auto v = query_varnodes[i];
    auto x = query_syms[i];
    auto r = h.outgoing(v, x);
  }

  // last_char
  for (auto v : query_varnodes) { h.lastchar(v); }
  #endif
}

