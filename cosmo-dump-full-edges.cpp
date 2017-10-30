#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph_shifted.hpp"
#include "algorithm.hpp"
#include "cosmo-dump-full-edges.hpp"
#include "sort.hpp"
//using namespace std;
//using namespace sdsl;

#include <sys/timeb.h>


int getMilliCount(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}


int getMilliSpan(int nTimeStart){
  int nSpan = getMilliCount() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}

void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);

  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();


}


void dump_edges(const debruijn_graph_shifted<> &dbg) {
  for (size_t i = 0; i < dbg.size(); i++) {
      cout <<  /* << "e:" <<*/ dbg.edge_label(i) << std::endl;
  }
}




int main(int argc, char* argv[]) {
    parameters_t p;
    parse_arguments(argc, argv, p);
    cerr << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
    //ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
    // Can add this to save a couple seconds off traversal - not really worth it.
    cerr << "loading dbg "<< p.input_filename  << std::endl;
    debruijn_graph_shifted<> dbg;
    load_from_file(dbg, p.input_filename);
    //input.close();

 
    cerr << "k             : " << dbg.k << endl;
    cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << dbg.num_edges() << endl;
    cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
    // dumpcolumns(dbg, p);
    dump_edges(dbg);


      

 }
