#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"
#include "algorithm.hpp"

using namespace std;
using namespace sdsl;

#include <sys/timeb.h>

int getMilliCount();
int getMilliCount(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}

int getMilliSpan(int nTimeStart);
int getMilliSpan(int nTimeStart){
  int nSpan = getMilliCount() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}

string extension = ".dbg";

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string output_prefix = "";
  std::string color_mask1 = "";
  std::string color_mask2 = "";
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
  string color_mask1 = "color_mask1";
  TCLAP::ValueArg<std::string> color_mask1_arg("a", "color_mask1",
	    "Color mask 1, color1 [" + color_mask1 + "]", false, "", color_mask1, cmd);
  string color_mask2 = "color_mask2";
  TCLAP::ValueArg<std::string> color_mask2_arg("b", "color_mask2",
	    "Color mask 2, color2 [" + color_mask2 + "]", false, "", color_mask2, cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  params.color_mask1     = color_mask1_arg.getValue();
  params.color_mask2     = color_mask2_arg.getValue();
}

static char base[] = {'?','A','C','G','T'};

void test_symmetry(debruijn_graph<> dbg);
void test_symmetry(debruijn_graph<> dbg) {
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


void dump_nodes(debruijn_graph<> dbg, uint64_t * colors);
void dump_nodes(debruijn_graph<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
    cout << i << ":" << dbg.node_label(i) << colors[dbg._node_to_edge(i)] << "\n";
  }
}

void dump_edges(debruijn_graph<> dbg, uint64_t * colors);
void dump_edges(debruijn_graph<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.num_edges(); i++) {
    cout << i << "e:" << dbg.edge_label(i) << colors[i] << "\n";
  }
}

void find_bubbles(debruijn_graph<> dbg, rrr_vector<63> &colors, uint64_t color_mask1, uint64_t color_mask2);
void find_bubbles(debruijn_graph<> dbg, rrr_vector<63> &colors, uint64_t color_mask1, uint64_t color_mask2) {
    int t = getMilliCount();
    int num_colors = colors.size() / dbg.num_edges();
    //uint64_t combined_mask = color_mask1 | color_mask2;
    bit_vector visited = bit_vector(dbg.num_nodes(), 0);
    cout << "Starting to look for bubbles\n";
    for (size_t i = 0; i < dbg.num_nodes(); i++) {
        // cout << "Node " << i << ":" << dbg.node_label(i) << " color: " << color_mask << "\n";
        if (!visited[i] && dbg.outdegree(i) == 2) { //FIXME: why do we only care about outdegree == 2?
            // initialize bubble tracking variables
            int branch_num = 0;
            ssize_t end[2]; // place to store end of branch kmer
//#define MAX_BRANCH 6000
            std::vector<std::string> branch(2);
            //char branch[2][MAX_BRANCH+1] = {{0}, {0}}; // the branch sequence for printing out later
            int branch_offset = 0;
            uint64_t branch_color[2];
            ssize_t start = i; // place to store start of branch kmer

            // start of a bubble handling
            for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
                // follow each strand or supernode
                ssize_t edge = dbg.outgoing_edge(i, x);
                if (edge == -1)
                    continue;
                branch[branch_num] += base[x];
                // build color mask
                uint64_t color_mask = 0;
                for (int c = 0; c < num_colors; c++)
                    color_mask |= colors[edge * num_colors + c] << c;
                branch_color[branch_num] = color_mask;

                // walk along edges until we encounter 
                ssize_t pos = dbg._edge_to_node(edge);
                while (dbg.indegree(pos) == 1 && dbg.outdegree(pos) == 1 /*&& branch_offset < MAX_BRANCH*/) {
                    visited[pos] = 1;
                    ssize_t next_edge = 0;
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(pos, x2);
                        if (next_edge != -1) {
                            branch[branch_num] += base[x2];
                            break;
                        }
                    }
                    pos = dbg._edge_to_node(next_edge);
                    //cout << pos << ":" << dbg.node_label(pos) << "\n";
                }
                // cout << "Stopped due to : " << dbg.indegree(pos) << ":" << dbg.outdegree(pos) << ":" << branch_offset << "\n";
                // if (branch_offset > MAX_BRANCH)
                //     break;

                end[branch_num++] =  (dbg.indegree(pos) > 1) ? pos : 0;
                branch_offset = 0;
            }
            // check if both branches ended on the same kmer and they pass the requested color masks
            //cout << "Trying " << branch_color[0] << ":" << branch_color[1] << " " << end[0] << ":" << end[1] <<"\n";
            //cout << color_mask1 << ":" << color_mask2 << "\n";
            //cout << "PutativeStart flank: " << dbg.node_label(start) << " c: " << branch_color[0] << ":" << branch_color[1] << "\n";

            if ((end[0] && end[0] == end[1]) &&
                ((color_mask1 & branch_color[0] && !(~color_mask1 & branch_color[0]) &&
                  color_mask2 & branch_color[1] && !(~color_mask2 & branch_color[1])) || 
                 (color_mask1 & branch_color[1] && !(~color_mask1 & branch_color[1]) &&
                  color_mask2 & branch_color[0] && !(~color_mask2 & branch_color[0])))) {
                cout << "\nStart flank: " << dbg.node_label(start) << " c: " << branch_color[0] << ":" << branch_color[1] << "\n";
                cout << "Branch: " << branch[0] << "\n";
                cout << "Branch: " << branch[1] << "\n";
                cout << "End flank: " << dbg.node_label(end[0]) << "\n";
            }
        }
    }
    cerr << "Find bubbles time: " << getMilliSpan(t) << "\n";
}


int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  // Can add this to save a couple seconds off traversal - not really worth it.
  //vector<size_t> minus_positions;
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
  input.close();

  rrr_vector<63> colors;
  load_from_file(colors, p.color_filename);

  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "colors        : " << colors.size() / dbg.num_edges() << endl; 
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
  cerr << "Color size    : " << size_in_mega_bytes(colors) << " MB" << endl;

  //dump_nodes(dbg, colors);
  //dump_edges(dbg, colors);
  uint64_t mask1 = (p.color_mask1.length() > 0) ? atoi(p.color_mask1.c_str()) : -1;
  uint64_t mask2 = (p.color_mask2.length() > 0) ? atoi(p.color_mask2.c_str()) : -1;
  find_bubbles(dbg, colors, mask1, mask2);
}
