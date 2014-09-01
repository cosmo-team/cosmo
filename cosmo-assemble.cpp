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
  debruijn_graph<> dbg;
  load_from_file(dbg, p.input_filename);
  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;

  // This actually does take a little while (10 sec) to build. Should maybe add flags to
  // pre-build it during graph construction and pass it in (faster)
  // TODO: TIME it on ch14
  // e.g. cosmo-build -u (unipaths)
  // cosmo-assemble -u unipaths-file
  sd_vector<> b = make_branch_vector(dbg);
  cerr << "Branch size   : " << size_in_mega_bytes(b) << " MB" << endl;

  // Traversal
  ofstream out;
  out.open(outfilename + ".fasta", ios::out);
  debruijn_graph<>::label_type s{};
  size_t threshold = 100;
  size_t id = 1;
  visit_unipaths(dbg, b, [&](char x) {
    if (x == '$') {
      if(s.length() >= threshold) {
        out << ">cosmo_" << id++ << endl;
        out << s << endl;
      }
      s = debruijn_graph<>::label_type{};
    }
    else s.push_back(x);
  });
  out.flush();
  out.close();
}

