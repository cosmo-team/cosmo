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
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from kramer).", true, "", "input_file", cmd);
  cmd.parse( argc, argv );
  // -d flag for decompression to original kmer biz
  params.input_filename  = input_filename_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$acgt");
  input.close();

  cout << "k              = " << dbg.k << endl;
  cout << "num_nodes()    = " << dbg.num_nodes() << endl;
  cout << "num_edges()    = " << dbg.num_edges() << endl;

  for (size_t node = 0; node < dbg.num_nodes()-1; node++) {
    size_t edge = dbg._node_to_edge(node);
    size_t indegree = dbg.indegree(node);
    cout << "indegree(" << node << "<" << edge << ">) = " << indegree << endl;
  }

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
  debruijn_graph<> dbg2;
  store_to_file(dbg, "dbg.dbg");
  load_from_file(dbg2, "dbg.dbg");
  cout << "Total size:    " << size_in_mega_bytes(dbg2) << " MB" << endl;
  cout << "Bits per edge: " << bits_per_element(dbg2) << " Bits" << endl;
  //cout << is_pod<array<size_t, 5>>::value << endl;
  //write_structure<JSON_FORMAT>(dbg2, cout);
}
