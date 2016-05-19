#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <libgen.h> // basename
#include <sys/mman.h> // mlockall

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
//#include <boost_iterator/zip_iterator.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"
#include "debruijn_graph_shifted.hpp"
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

  auto base_path = boost::filesystem::path(p.input_filename).parent_path().string();
  auto base_name = base_path + "/" + boost::filesystem::path(p.input_filename).stem().string();
  COSMO_LOG(debug) << base_name;

  // TO LOAD:
  debruijn_graph_shifted<> g;
  load_from_file(g, p.input_filename);

  auto dbg_size = size_in_mega_bytes(g);
  cerr << "k             : " << g.k << endl;
  cerr << "num_nodes()   : " << g.num_nodes() << endl;
  cerr << "num_edges()   : " << g.num_edges() << endl;
  cerr << "W size        : " << size_in_mega_bytes(g.m_edges) << " MB" << endl;
  cerr << "L size        : " << size_in_mega_bytes(g.m_node_flags) << " MB" << endl;
  cerr << "DBG size      : " << dbg_size << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(g) << " Bits" << endl;

  #ifdef VAR_ORDER
  wt_int<rrr_vector<63>> lcs;
  load_from_file(lcs, base_name + ".lcs");

  auto lcs_size = size_in_mega_bytes(lcs);
  auto total_size = dbg_size + lcs_size;
  cerr << "LCS size      : " << lcs_size << " MB" << endl;
  cerr << "LCS bits/edge : " << bits_per_element(lcs) << " Bits" << endl;
  cerr << "Total size    : " << total_size << " MB" << endl;

  typedef debruijn_hypergraph<> dbh;
  typedef dbh::node_type node_type;
  dbh h(g, lcs);
  #endif

  /*
  cerr << "symbol ends:" << endl;
  for (size_t i = 0; i<5; ++i) {
    cerr << i << ": " << g.m_symbol_ends[i] << endl;
  }
  */

  int num_queries = 2e4;
  size_t min_k = 8; // this could be 0, but it affects longer too much
  size_t max_k = g.k-2; // K is the edge length, and K-1 node lenght.
  // Time will be identical if we measure K-1

  cerr << "Generating " << num_queries << " random queries" << endl;
  // set up RNGs
  typedef boost::mt19937 rng_type;
  rng_type rng(time(0));
  boost::uniform_int<size_t> node_distribution(0, g.num_nodes()-1); // Randomly choose nodes from the graph
  boost::uniform_int<size_t> k_distribution(min_k, max_k);
  boost::uniform_int<size_t> symbol_distribution(1, 4); // A C G T in our WT (ignoring $=0)
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_node(rng, node_distribution);
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_k(rng, k_distribution);
  boost::variate_generator<rng_type, boost::uniform_int<size_t>> random_symbol(rng, symbol_distribution);

  vector<size_t> query_nodes(boost::make_function_input_iterator(random_node,0),
                             boost::make_function_input_iterator(random_node,num_queries));
  vector<size_t> query_syms(boost::make_function_input_iterator(random_symbol,0),
                            boost::make_function_input_iterator(random_symbol,num_queries));

  #ifdef VAR_ORDER
  vector<node_type> query_varnodes;
  for (auto u:query_nodes) {
    query_varnodes.push_back(h.get_node(u));
  }

  // Randomly change the order of random nodes
  // TODO: WHY is this in this range?
  // -2 because we need to remove the last character
  auto random_low_bound_k = [&](size_t low)  {
    assert(low >= min_k);
    return boost::uniform_int<size_t>(low, g.k-1)(rng);
  };
  auto random_high_bound_k  = [&](size_t high) {
    assert(high <= g.k-1);
    return boost::uniform_int<size_t>(min_k, high)(rng);
  };

  // Make random nodes that have an order that we can increase a certain amount
  const vector<size_t> query_order_deltas{1,2,4,8}; // Try all ks?
  size_t num_deltas = query_order_deltas.size();
  vector<vector<node_type>> shorter_query_varnodes(query_order_deltas.size());
  vector<vector<node_type>> longer_query_varnodes(query_order_deltas.size());

  for (size_t i=0; i<query_order_deltas.size();++i) {
    shorter_query_varnodes[i].reserve(num_queries);//= vector<node_type>(num_queries);
    longer_query_varnodes[i].reserve(num_queries);//= vector<node_type>(num_queries);
    for (size_t node_idx=0; node_idx < (size_t)num_queries; ++node_idx) {
      auto delta = query_order_deltas[i];
      auto node = query_varnodes[node_idx];
      // Nodes that are of at least a certain k so that we can shorten them later
      // Note that we dont adhere to the max_k since these operations
      // arent possible on the standard dbg
      auto shorter_k = random_low_bound_k(min_k+delta);
      assert(min_k+delta <= shorter_k && shorter_k < g.k);
      auto longer_k = random_high_bound_k(g.k-1-delta);
      node_type new_shorter_node;
      if (shorter_k == g.k-1) new_shorter_node = node;
      else new_shorter_node = h.shorter(node, shorter_k);
      shorter_query_varnodes[i].push_back(new_shorter_node);
      longer_query_varnodes[i].push_back(h.shorter(node, longer_k));
    }
  }

  // Convert to variable order nodes
  // Need: random lows for longer of each delta
  // random highs for shorter of each delta
  // random nodes for maxlen, maxlen*, backward, forward, etc
  vector<size_t> query_ks(boost::make_function_input_iterator(random_k,0),
                          boost::make_function_input_iterator(random_k,num_queries));
  vector<size_t> maxlen_syms(boost::make_function_input_iterator(random_symbol,0),
                            boost::make_function_input_iterator(random_symbol,num_queries));

    // set shorter for each one that has a shorter k
  for (int i=0;i<num_queries;i++) {
    if (query_ks[i] == g.k-1) continue;
    query_varnodes[i] = h.shorter(query_varnodes[i], query_ks[i]);
  }
  #else

  vector<debruijn_graph_shifted<>::node_type> query_rangenodes;
  //transform(query_nodes.begin(), query_nodes.end(), query_varnodes.begin(),[&](size_t v){ return h.get_node(v); });
  for (auto u:query_nodes) {
    query_rangenodes.push_back(g.get_node(u));
  }
  #endif

  mlockall(MCL_CURRENT);
  typedef chrono::nanoseconds unit;
  string unit_s = " ns";

  #ifndef VAR_ORDER // standard dbg
  // backward
  auto t1 = chrono::high_resolution_clock::now();
  for (auto v : query_rangenodes) {
    g.all_preds(v);
  }
  auto t2 = chrono::high_resolution_clock::now();
  auto dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "backward total : " << dur << " ns" <<endl;
  cerr << "backward : " << (double)dur/num_queries << unit_s <<endl;

  // forward
  t1 = chrono::high_resolution_clock::now();
  for (size_t i=0;i<(size_t)num_queries;i++) {
    auto v = query_rangenodes[i];
    auto x = query_syms[i];
    g.interval_node_outgoing(v, x);
  }
  t2 = chrono::high_resolution_clock::now();
  dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "forward total : " << dur << " ns" <<endl;
  cerr << "forward : " << (double)dur/num_queries << unit_s <<endl;

  // last_char
  t1 = chrono::high_resolution_clock::now();
  for (auto v : query_rangenodes) { g.lastchar(v); }
  t2 = chrono::high_resolution_clock::now();
  dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "lastchar total : " << dur << " ns" <<endl;
  cerr << "lastchar : " << (double)dur/num_queries << unit_s <<endl;
  #else
  auto t1 = chrono::high_resolution_clock::now();
  // backward
  for (auto v : query_varnodes) {
    h.backward(v);
  }
  auto t2 = chrono::high_resolution_clock::now();
  auto dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "backward total : " << dur << " ns" <<endl;
  cerr << "backward  : " << (double)dur/num_queries << unit_s <<endl;

  // forward
  t1 = chrono::high_resolution_clock::now();
  for (size_t i=0;i<(size_t)num_queries;i++) {
    auto v = query_varnodes[i];
    auto x = query_syms[i];
    //cerr << "("<< get<0>(v) << ", " << get<1>(v) << ", "<< get<2>(v) << ") " << (int)x << endl;
    h.outgoing(v, x);
  }
  t2 = chrono::high_resolution_clock::now();
  dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "forward total : " << dur << " ns" <<endl;
  cerr << "forward   : " << (double)dur/num_queries << unit_s <<endl;

  // last_char
  t1 = chrono::high_resolution_clock::now();
  for (auto v : query_varnodes) { h.lastchar(v); }
  t2 = chrono::high_resolution_clock::now();
  dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "lastchar total : " << dur << " ns" <<endl;
  cerr << "lastchar  : " << (double)dur/num_queries << unit_s <<endl;

  // shorter
  for (size_t i=0; i<num_deltas; ++i) {
    size_t delta = query_order_deltas[i];
    t1 = chrono::high_resolution_clock::now();
    for (size_t idx=0;idx<(size_t)num_queries;idx++) {
      auto v = shorter_query_varnodes[i][idx];
      size_t k = get<2>(v);
      size_t new_k = k-delta;
      assert(new_k >= min_k);
      h.shorter(v, new_k);
    }
    t2 = chrono::high_resolution_clock::now();
    dur = chrono::duration_cast<unit>(t2-t1).count();
    cerr << "shorter_" << delta << " : " << (double)dur/(num_queries) << unit_s <<endl;
  }

  // longer
  for (size_t i=0; i<num_deltas; ++i) {
    size_t delta = query_order_deltas[i];
    t1 = chrono::high_resolution_clock::now();
    size_t total = 0;
    for (size_t idx=0;idx<(size_t)num_queries;idx++) {
      auto v = longer_query_varnodes[i][idx];
      size_t k = get<2>(v);
      size_t new_k = k+delta;
      assert(new_k >= min_k);
      assert(new_k < g.k);
      auto result = h.longer(v, new_k);
      total += result.size();
    }
    t2 = chrono::high_resolution_clock::now();
    dur = chrono::duration_cast<unit>(t2-t1).count();
    cerr << "longer_" << delta << "  : " << (double)dur/(num_queries) << unit_s <<endl;
    cerr << "per node  : " << (double)dur/(total) << unit_s << " (" << total << ")" << endl;
  }


  /*
  size_t skipped = 0;
  for (size_t k : {1,2,4,8}) {
    t1 = chrono::high_resolution_clock::now();
    for (size_t i=0;i<(size_t)num_queries;i++) {
      auto v = query_varnodes[i];
      auto k = higher_ks[i];
      if (get<2>(v) > k) {
        skipped++;
        continue;
      }
      //cout << get<0>(v) << ", " << get<1>(v) << ", " << get<2>(v) << endl;
      h.longer(v, k);
    }
    t2 = chrono::high_resolution_clock::now();
    dur = chrono::duration_cast<unit>(t2-t1).count();
    //cerr << "longer(v,+"<<k<<") total : " << dur << " ns" <<endl;
    cerr << "longer"<<" mean  : " << (double)dur/(num_queries-skipped) << unit_s <<endl;
  //}
  */

  // maxlen with symbol
  t1 = chrono::high_resolution_clock::now();
  for (size_t i=0;i<(size_t)num_queries;i++) {
    auto v = query_varnodes[i];
    auto x = maxlen_syms[i];
    h.maxlen(v, x);
  }
  t2 = chrono::high_resolution_clock::now();
  dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "maxlen(v,c) total : " << dur << " ns" <<endl;
  cerr << "maxlen    : " << (double)dur/num_queries << unit_s <<endl;

  // Regular maxlen
  t1 = chrono::high_resolution_clock::now();
  for (auto v : query_varnodes) { h.maxlen(v); }
  t2 = chrono::high_resolution_clock::now();
  dur = chrono::duration_cast<unit>(t2-t1).count();
  //cerr << "maxlen(v,*) total : " << dur << " ns" <<endl;
  cerr << "maxlen*   : " << (double)dur/num_queries << unit_s <<endl;
  #endif
  munlockall();
}

