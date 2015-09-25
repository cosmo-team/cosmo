#include <iostream>
#include <fstream>

#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"
#include "algorithm.hpp"

using namespace std;
using namespace sdsl;

string extension = ".dbg";

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string output_prefix = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color",
            ".color file (output from pack-edges).", true, "", "color_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}

void find_bubbles(debruijn_graph<> dbg, uint64_t * colors);
void find_bubbles(debruijn_graph<> dbg, uint64_t * colors) {
  //char visited[dbg.num_nodes()];
  //bzero(visited, dbg.num_nodes());
  printf("Got colors %llx \n", *colors);
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
    cout << dbg.node_label(i) << "\n";
      /*
    if (dbg.outdegree(i) > 1) {
      // start of a bubble
      cout << dbg.node_label(i) << "\n";
      for (unsigned long x = 0; x<dbg.sigma+1;x++) {
	// follow each strand
	ssize_t pos = dbg.outgoing(i, x);
	if (pos == -1)
	  continue;
	int len = 1;
	while (dbg.indegree(pos) == 1) {
	  printf("pos is %zu\n", pos);
	  pos = dbg._edge_to_node(dbg._forward(pos));
	  len++;
	}
	printf("Found bubble from %zu to %zu with %d length\n", i, pos, len); 
      }
    }
    */
  }
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  // Can add this to save a couple seconds off traversal - not really worth it.
  //vector<size_t> minus_positions;
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
  input.close();

  ifstream colorfile(p.color_filename, ios::in|ios::binary|ios::ate);
  uint64_t * colors = (uint64_t *) malloc(dbg.num_edges() * sizeof(uint64_t));
  colorfile.read((char *)colors, dbg.num_edges() * sizeof(uint64_t));
  colorfile.close();

  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;

  find_bubbles(dbg, colors);

  // The parameter should be const... On my computer the parameter
  // isn't const though, yet it doesn't modify the string...
  // This is still done AFTER loading the file just in case
  char * base_name = basename(const_cast<char*>(p.input_filename.c_str()));
  string outfilename = ((p.output_prefix == "")? base_name : p.output_prefix) + extension;
  store_to_file(dbg, outfilename);

  #ifdef VAR_ORDER
  wt_int<rrr_vector<63>> lcs;
  construct(lcs, base_name + string(".lcs"), 1);
  cerr << "LCS size      : " << size_in_mega_bytes(lcs) << " MB" << endl;
  cerr << "LCS bits/edge : " << bits_per_element(lcs) << " Bits" << endl;
  store_to_file(lcs, outfilename + ".lcs.wt");
  // TODO: Write compressed LCS
  #endif
}
