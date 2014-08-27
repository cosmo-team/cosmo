#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <array>
#include <type_traits>

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"

using namespace std;
using namespace sdsl;

typedef uint64_t kmer_t;

struct parameters_t {
  std::string input_filename = "";
};

#define VERSION "0.0"
void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("SDBG-Lite Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from kramer).", true, "", "input_file", cmd);
  cmd.parse( argc, argv );
  // -d flag for decompression to original kmer biz
  params.input_filename  = input_filename_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  cout << "Opening " << p.input_filename << "..." << endl;
  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);

  /*
  // check length
  streampos filesize = input.tellg();
  if (filesize/sizeof(uint64_t) < 6) {
    cerr << "Error: file doesn't have necessary footer information (cumulative symbol counts + k)." << endl;
    exit(1);
  }

  input.seekg(-6 * sizeof(uint64_t), ios::end);
  uint64_t counts[5] = {0, 0, 0, 0, 0};
  uint64_t k = 0;
  input.read((char*)&counts[0], 5 * sizeof(uint64_t));
  input.read((char*)&k, sizeof(uint64_t));
  size_t num_edges = counts[4];

  cout << "counts = [ ";
  cout << counts[0] << ", ";
  cout << counts[1] << ", ";
  cout << counts[2] << ", ";
  cout << counts[3] << ", ";
  cout << counts[4] << " ]" << endl;
  cout << "k = " << k << endl;

  size_t num_blocks = size_t(filesize)/sizeof(uint64_t) - 6;
  input.seekg(0, ios::beg); // rewind

  // Loop through this while reading from file instead
  vector<uint64_t> blocks(num_blocks,0);
  input.read((char*)&blocks[0], sizeof(uint64_t) * num_blocks);
  input.close();

  bit_vector first(num_edges,0);
  // would be nice to fix wavelet trees so the constructor
  // can accept a int_vector<4> instead (which is all we need)
  int_vector<8> edges(num_edges);

  for (size_t i = 0; i < num_edges; i++) {
    auto x = get_edge(blocks.begin(), i);
    //cout << get<1>(x) << endl;// << " " << "$acgt"[get<0>(x)] << " " << get<2>(x) << endl;
    first[i] = !get<1>(x); // convert 0s to 1s so we can have a sparse bit vector
    edges[i] = (get<0>(x) << 1) | !get<2>(x);
    //cout << "$acgt"[edges[i] >> 1] <<endl;
  }

  //rrr_vector<67> rrr(first); // pass to debruijn_graph as param
  sd_vector<> sbv(first);
  // These are swapped because we invert the bits to make it a sparse bit vector instead of a dense one
  // (which sdsl doesnt seem to have an implementation for, fair enough though...)
  sd_vector<>::rank_0_type   node_rank(&sbv);
  sd_vector<>::select_0_type node_select(&sbv);

  wt_huff<rrr_vector<67>> wt;
  construct_im(wt, edges);
  cout << int(wt[1]) << endl;
  // the ORing here is to add a minus flag
  cout << wt.select(3,(3<<1)|1)+1 << endl;
  cout << wt.rank(9,(1<<1)) << endl;

  //cout << "Size of bit_vector in MB: " << size_in_mega_bytes(first) << endl;
  //cout << "Size of rrr_vector in MB: " << size_in_mega_bytes(rrr) << endl;
  //cout << "Size of int_vector in MB: " << size_in_mega_bytes(edges) << endl;
  cout << "Size of wavelet_t in MB : " << size_in_mega_bytes(wt) << endl;
  cout << "Size of sd_vector in  B : " << size_in_bytes(sbv) << endl;
  cout << "Size of   sd_rank in  B : " << size_in_bytes(node_rank) << endl;
  cout << "Size of sd_select in  B : " << size_in_bytes(node_select) << endl;
  size_t total_bits = (size_in_bytes(sbv) + size_in_bytes(wt) + 4*sizeof(uint64_t))*8;
  cout << "          Bits per edge : " << double(total_bits)/num_edges << endl;
  */
  /*
  cout << "k              = " << dbg.k << endl;
  cout << "num_nodes()    = " << dbg.num_nodes() << endl;
  cout << "num_edges()    = " << dbg.num_edges() << endl;
  cout << "node_label(4)  = " << dbg.node_label(4) << endl;
  cout << "node_label(-1) = " << dbg.node_label(dbg.num_nodes()-1) << endl;
  cout << "edge_label(0)  = " << dbg.edge_label(0) << endl;
  cout << "edge_label(1)  = " << dbg.edge_label(1) << endl;
  cout << "edge_label(-1) = " << dbg.edge_label(dbg.num_edges()-1) << endl;
  */

  /*
  for (size_t node = 0; node < dbg.num_nodes()-1; node++) {
    size_t edge = dbg._node_to_edge(node);
    size_t indegree = dbg.indegree(node);
    cout << "indegree(" << node << "<" << edge << ">) = " << indegree << endl;
  }
  */

  /*
  size_t edge = 1262;
  size_t node = dbg._edge_to_node(edge);
  for (uint8_t x = 0; x < 5; x++) {
    ssize_t result = dbg.incoming(node, x);
    cout << "dbg.incoming(" << node << "<" << edge << ">, " << "$acgt"[x] << ") = " << result << "<";
    if (result==-1) cout << "_";
    else cout << dbg._node_to_edge(result);
    cout << ">" << endl;
  }
  */
  /*
  array<size_t, 5> a{};
  a[0] = 0; a[1] = 1; a[2] = 2; a[3] = 2; a[4] = 4;
  store_to_file(a, "array.sdsl");
  array<size_t, 5> b{};
  load_from_file(b, "array.sdsl");
  for (auto &x : b) {
    cout << x << endl;
  }
  */
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$acgt");
  debruijn_graph<> dbg2;
  store_to_file(dbg, "dbg.dbg");
  load_from_file(dbg2, "dbg.dbg");
  cout << "Total size:    " << size_in_mega_bytes(dbg2) << " MB" << endl;
  cout << "Bits per edge: " << bits_per_element(dbg2) << " Bits" << endl;
  //cout << is_pod<array<size_t, 5>>::value << endl;
  //write_structure<JSON_FORMAT>(dbg2, cout);
}
