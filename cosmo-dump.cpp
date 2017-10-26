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
#include "cosmo-dump.hpp"

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

static char base[] = {'?','A','C','G','T'};


void test_symmetry(debruijn_graph_shifted<> dbg) {
  for (unsigned long x = 0; x<dbg.sigma+1;x++) {
    ssize_t in = dbg.incoming(43, x);
    if (in == -1)
      continue;
    for (unsigned long y = 0; y<dbg.sigma+1;y++) {
      ssize_t out = dbg.outgoing(in, y);
      if (out == -1)
	continue;
      cout << "Incoming " << in <<  ":" << out <<"\n";
    }
  }
}



void dump_nodes(debruijn_graph_shifted<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
    cout << i << ":" << dbg.node_label(i) << colors[dbg._node_to_edge(i)] << "\n";
  }
}


void dump_edges(debruijn_graph_shifted<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.size(); i++) {
    cout << i << "e:" << dbg.edge_label(i) << colors[i] << "\n";
  }
}





int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  cerr << "loading dbg" << std::endl;
  debruijn_graph_shifted<> dbg;
  load_from_file(dbg, p.input_filename);





  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;

  std::cout << "edges:" <<std::endl;
  for (uint64 i = 0; i < dbg.num_edges(); ++i) {
      std::cout << dbg._map_symbol(dbg._strip_edge_flag(dbg.m_edges[i])) << std::endl;
  }
  std::cout << std::endl;
  
  std::cout << "flags:" << std::endl;
  for (uint64 i = 0; i < dbg.num_edges(); ++i) {
      std::cout << (int)(dbg.m_edges[i] & 1) << std::endl;
  }
  std::cout << std::endl;

  std::cout << "node flags:" << std::endl;
  for (uint64 i = 0; i < dbg.num_edges(); ++i) {
      std::cout << (int)(dbg.m_node_flags[i]) << std::endl;
  }
  std::cout << std::endl;

  std::cout << "m_symbol_ends:" << std::endl;
  for (uint64 i = 0; i < dbg.sigma + 1; ++i) {
      std::cout << (dbg.m_symbol_ends[i]) << std::endl;
  }
 

}
