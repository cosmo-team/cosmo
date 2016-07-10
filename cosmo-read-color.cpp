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
#include "cosmo-color-pd.hpp"
#include <future>
using namespace std;
using namespace sdsl;

#include <sys/timeb.h>



//////////// pack-color headers


#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>

// TCLAP
#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <cstdio>

#include <cstdlib>

#include <libgen.h>

// Custom Headers
//#include "uint128_t.hpp"
//#include "debug.h"
#include "kmer.hpp"

//using namespace std;
//using namespace sdsl;

#include <cstdlib>
#include <sys/timeb.h>
//#include "pack-color.hpp"








bool trace = false;
string file_extension = ".dbg";

unsigned long long global_perseq_t;
unsigned long long global_t;
static char base[] = {'?','A','C','G','T'};
char dna_bases[] = "$ACGT";
//unsigned colorgroups[] = {1, 25, 49, 57, 64, 80};
unsigned colorgroups[] = {1, 25, 49, 57, 72, 88};
unsigned num_colorgroups = 5;
debruijn_graph_shifted<>* gdbg;
//sd_vector<>* gcolors;
//rank_support_sd<>* gcolor_ranks;
const char *const starts[] = {"GCCATACTGCGTCATGTCGCCCTGACGCGC","GCAGGTTCGAATCCTGCACGACCCACCAAT","GCTTAACCTCACAACCCGAAGATGTTTCTT","AAAACCCGCCGAAGCGGGTTTTTACGTAAA","AATCCTGCACGACCCACCAGTTTTAACATC","AGAGTTCCCCGCGCCAGCGGGGATAAACCG","GAATACGTGCGCAACAACCGTCTTCCGGAG"};

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
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color",
            ".color file (output from pack-edges).", true, "", "color_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + file_extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  string ref_color = "ref_color";
  TCLAP::ValueArg<std::string> ref_color_arg("a", "ref_color",
	    "Ref color, ref_color [" + ref_color + "]", false, "", ref_color, cmd);
  string sample_mask = "sample_mask";
  TCLAP::ValueArg<std::string> sample_mask_arg("b", "sample_mask",
	    "Sample mask, sample_mask [" + sample_mask + "]", false, "", sample_mask, cmd);
  string ref_fasta = "ref_fasta";
  TCLAP::ValueArg<std::string> ref_fasta_arg("r", "ref_fasta",
	    "Reference FASTA filename, ref_fasta [" + ref_fasta + "]", false, "", ref_fasta, cmd);

  string output_matrix = "output_matrix";
  TCLAP::ValueArg<std::string> output_matrix_arg("m", "output_matrix",
	    "Reference FASTA filename, output_matrix [" + output_matrix + "]", false, "", output_matrix, cmd);
  
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  params.ref_color     = ref_color_arg.getValue();
  params.sample_mask     = sample_mask_arg.getValue();
    params.ref_fasta = ref_fasta_arg.getValue();
  params.output_matrix = output_matrix_arg.getValue();
}




void test_symmetry(debruijn_graph_shifted<>& dbg)
{
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



void dump_nodes(debruijn_graph_shifted<>& dbg, uint64_t * colors)
{
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
    cout << i << ":" << dbg.node_label(i) << colors[dbg._node_to_edge(i)] << "\n";
  }
}


void dump_edges(debruijn_graph_shifted<>& dbg, uint64_t * colors)
{
  for (size_t i = 0; i < dbg.num_edges(); i++) {
    cout << i << "e:" << dbg.edge_label(i) << colors[i] << "\n";
  }
}


ssize_t get_first_node(debruijn_graph_shifted<>& dbg, sd_vector<> &colors, uint64_t ref_color, std::string& ref_fasta_content)
{

    
    // find the first edge which is colored 'ref_color' and whose label we can find in the ref genome
    // (since we generate revcomps and color them the same as given, the first edge's label may not
    //  be found in the ref genome; we'll just keep trying till we find one)
    int num_colors = colors.size() / dbg.num_edges();
    ssize_t zeroth_rank_edge = (unsigned long long)-1;
    // size_t pos =  std::string::npos;
    ssize_t node_num = 0;
    std::string node_label;
    std::string query = ref_fasta_content.substr(0, dbg.node_label(0).size()/*hack to get the edge kmer size*/);
    for (; node_num < (size_t)dbg.num_edges(); ++ node_num) {
        if (colors[node_num * num_colors + ref_color]) {
            zeroth_rank_edge = node_num;
            node_label = dbg.node_label(zeroth_rank_edge);
            if (node_label == query) {
                //if ((pos = ref_fasta_content.find(edge_label)) == 0 /*!= std::string::npos*/){
                std::cerr << "Found fasta start at edge num " << node_num << std::endl;
                return node_num;
                break;
            }
        }
    }

    // if (pos != std::string::npos) {
    //     std::cout << "Found edge " << edge_label << " number " << edge_num << " in reference genome at position " << pos << "." << std::endl;
    //     std::cout << "ref_genome[" << pos - 1 << ":" << pos + edge_label.size() << " = "
    //               << ref_fasta_content.substr(pos - 1, 1 + edge_label.size()) << std::endl;
    // }else{
    //     std::cerr << "ERROR: Can't find any ref color edges in the ref FASTA" << std::endl;
    //     exit(EXIT_FAILURE);
    // }

    // while (pos != 0) {
        
        

    // now traverse the graph backward till we get back to the start of the ref genome
    //while (pos != 0) {
        

    return -1;

}

unsigned dna_ord(char c)
{
    switch(c) {
    case 'A' : return 1;
    case 'C': return 2;
    case 'G': return 3;
    case 'T': return 4;
    default: std::cerr << "ERROR: Unknown base '" << c << "'." << std::endl;
        exit(EXIT_FAILURE);
        
    };
}

void advance(debruijn_graph_shifted<>& dbg, const std::string& ref_fasta_content, const unsigned amount, ssize_t& node_k, ssize_t& node_k_pos)
{
    unsigned node_label_size = dbg.k - 1;  
    for (unsigned i = 0; i < amount; ++i) {
        
        ssize_t edge = dbg.outgoing_edge(node_k, dna_ord(ref_fasta_content[node_k_pos + node_label_size]));
        if (edge == -1) {
            std::cerr << "Reference fasta guides graph traversal through invalid path at position " << node_k_pos + node_label_size << "." << std::endl;
            exit(EXIT_FAILURE);
        }
        node_k = dbg._edge_to_node(edge);
        node_k_pos++;
    }
    
    
}

    
void dumping_advance(debruijn_graph_shifted<>& dbg, const std::string& ref_fasta_content, const unsigned amount, ssize_t& node_k, ssize_t& node_k_pos, rank_support_sd<1>& color_ranks, std::vector<unsigned>& group_counts, unsigned num_colors, std::stringstream& ret)
{
    //std::cout << dbg.node_label(node_k) << " ";
    unsigned node_label_size = dbg.k - 1;  
    for (unsigned i = 0; i < amount; ++i) {
        
        ssize_t edge = dbg.outgoing_edge(node_k, dna_ord(ref_fasta_content[node_k_pos + node_label_size]));
        if (edge == -1) {
            ret <<  ref_fasta_content.substr(node_k_pos, node_label_size) << " Reference fasta guides graph traversal through invalid path at position " << node_k_pos + node_label_size ;
        }

        for (unsigned i=0; i< num_colorgroups; ++i) {
            group_counts[i] += color_ranks(edge * num_colors + colorgroups[i+1]-1) - color_ranks(edge * num_colors + colorgroups[i]-1);
        }

        //std::cout << "color group counts: ";
        // for (unsigned i=0; i< num_colorgroups; ++i) {
        //     std::cout << group_counts[i] / (float)(colorgroups[i+1] - colorgroups[i]) / (float)ref_fasta_content.size() << " ";
        // }
        
        //std::cout << " time: " << (getMilliCount() - global_perseq_t) / 1000.0 << " s"  << std::endl;
        

        node_k = dbg._edge_to_node(edge);
        node_k_pos++;
    }
    
    
}

inline void dumping_advance2(debruijn_graph_shifted<>& dbg, const std::string& ref_fasta_content, const unsigned amount, ssize_t& node_k, ssize_t& node_k_pos,  std::vector<unsigned>& group_counts, std::vector<std::set<unsigned long long> >& found_kmers)
{
    //std::cout << dbg.node_label(node_k) << " ";
    unsigned node_label_size = dbg.k - 1;  

        
    ssize_t edge = dbg.outgoing_edge(node_k, dna_ord(ref_fasta_content[node_k_pos + node_label_size]));
    if (edge == -1) {
        std::cerr <<  ref_fasta_content.substr(node_k_pos, node_label_size) << " Reference fasta guides graph traversal through invalid path at position "
                  << node_k_pos + node_label_size // we want the right end of the current position k-mer
                  << std::endl ;
    }
    if ((found_kmers.end()-1)->find(edge) != (found_kmers.end()-1)->end()) {
        std::set<unsigned long long> newset;
        found_kmers.push_back(newset);
    }
    (found_kmers.end()-1)->insert(edge);        
    
    node_k = dbg._edge_to_node(edge);
    node_k_pos++;

    
    
}

std::vector<std::set<unsigned long long> > walk_refs(std::string ref_fasta_content )
{
    std::vector<std::set<unsigned long long> > found_kmers;
    std::set<unsigned long long> initial;
    found_kmers.push_back(initial);
        
    unsigned long long t = getMilliCount();

    unsigned node_label_size = gdbg->k - 1;

    std::string first_edge_label(ref_fasta_content.substr(0, node_label_size + 1));
    
    auto edge = gdbg->index(first_edge_label.begin());

    if (!edge) {
        assert(false);
        //return  "ERROR: could not locate the first node specified by the reference sequence in the graph!";
    }
//    found_kmers.insert(get<0>(*node));
    ssize_t first_node =  gdbg->_edge_to_node(get<0>(*edge));

    
    ssize_t node_i = first_node; // cdbg node labeled with a k-mer existing in the reference sequence
    ssize_t node_i_pos = 0;  // starting position in the reference sequence for the above k-mer

    std::vector<unsigned> group_counts(num_colorgroups, 0);
    
    while(node_i_pos < (ssize_t)ref_fasta_content.size() - gdbg->k + 1 ) {
        dumping_advance2(*gdbg, ref_fasta_content, 1, node_i, node_i_pos, group_counts,  found_kmers);
        
    }

    return found_kmers;

}






void dump_node(debruijn_graph_shifted<>& dbg, sd_vector<> &colors, ssize_t v)
{
    std::cout << dbg.node_label(v) ;
    int num_colors = colors.size() / dbg.num_edges();

    
    //int c = 1;
    
    // out edges
    std::cout << " orientation 0 { ";
    // for each outgoing edge
    std::set<char> outedges;
    unsigned long long last_outgoing_symbol = 0;

    for (unsigned long symbol_iter = 1; symbol_iter < dbg.sigma + 1; symbol_iter++) {

            // if there exists an outgoing edge for that symbol
            ssize_t outgoing_edge = dbg.outgoing_edge(v, symbol_iter);
            // for each color
//            for (int c = 0; c < 6; ++c) {
                if (outgoing_edge != -1) {
                    last_outgoing_symbol = symbol_iter;
                    uint64_t node_colors = 0;
                    for (int c = 0; c < num_colors; c++)
                        node_colors |= colors[outgoing_edge * num_colors + c] << c;

            // and if any colors of that edge match the sample set of colors, increment the out degree counter
//                    if (node_colors & sample_mask) {
                    std::cout << dna_bases[symbol_iter] << node_colors << "(" << outgoing_edge <<")";
                    
                    //                  if (colors[outgoing_edge * num_colors + c]) {
                        //outedges.insert(dna_bases[symbol_iter]);
//                    }
//                }
            }
        }

    // for(std::set<char>::iterator it = outedges.begin(); it != outedges.end(); ++it)
    //     std::cout << *it;

    

    // in edges
    std::cout << "}";
    std::cout << " orientation 1 { ";
    std::set<char> inedges;
    for (unsigned long symbol_iter = 1; symbol_iter < dbg.sigma + 1; symbol_iter++) {
            ssize_t incoming_edge = dbg.incoming(v, symbol_iter);
            if (incoming_edge != -1) {
            // if there exists an outgoing edge for that symbol
                ssize_t last_outgoing_edge = dbg.outgoing_edge(v, last_outgoing_symbol);
                std::string incoming_label = dbg.edge_label(incoming_edge);
                std::string outgoing_label = dbg.edge_label(last_outgoing_symbol);
                
                if (last_outgoing_edge != -1 && (outgoing_label.substr(1).compare(incoming_label.substr(0, incoming_label.size() - 1)))) {
                    std::cout << std::endl << "ERROR: incompatible edges incident to the same node" << std::endl;
                    std::cout << " outgoing edge:" << outgoing_label << std::endl;
                    std::cout << " incoming edge:" << incoming_label << std::endl;
                    assert(!"invalid incoming/outgoing edge combination\n");
                }
                uint64_t node_colors = 0;
                for (int c = 0; c < num_colors; c++)
                    node_colors |= colors[incoming_edge * num_colors + c] << c;

                std::cout << dna_bases[symbol_iter] << node_colors << "(" << incoming_edge <<")";
                    
                //              for (int c = 0; c < 6; ++c) {
//                    if( colors[incoming_edge * num_colors + c]) {
//                        inedges.insert(dna_bases[symbol_iter]);
//                    }
//                }
            }
    }

    // for(std::set<char>::iterator it = inedges.begin(); it != inedges.end(); ++it)
    //     std::cout << *it;


    std::cout << "}" << std::endl;

}

void dump_graph_shifted(debruijn_graph_shifted<>& dbg, sd_vector<> &colors)
{

    
    for (size_t node_i = 0; node_i < dbg.num_nodes(); node_i++) {
        dump_node(dbg, colors, node_i);
    }
}


    
void find_bubbles(debruijn_graph_shifted<>& dbg, sd_vector<> &colors, uint64_t ref_color, uint64_t sample_mask)
{
    int t = getMilliCount();
    int num_colors = colors.size() / dbg.num_edges();
    //uint64_t combined_mask = ref_color | sample_mask;
    bit_vector visited = bit_vector(dbg.num_nodes(), 0);
    cout << "Starting to look for bubbles\n";
    std::vector<std::string> branch(2);
    bool found_miss = false;
    for (size_t node_i = 0; node_i < dbg.num_nodes(); node_i++) {
        ssize_t start_node = node_i; // place to store start of branch kmer
        std::string start_label(dbg.node_label(start_node));
        found_miss = false;
        // for (int si = 0; si < 7; ++si) {
        //     if (!start_label.compare(starts[si])) {
        //         std::cerr << "Found missing start node " << starts[si] << " outdegree: " << dbg.outdegree(i) << std::endl;
        //         found_miss = true;
        //     }
        // }

        // cout << "Node " << i << ":" << dbg.node_label(i) << " color: " << color_mask << "\n";
        if (!visited[node_i] && dbg.outdegree(node_i) == 2) { //FIXME: why do we only care about outdegree == 2?
            // initialize bubble tracking variables
            int branch_num = 0;
            ssize_t end[2]; // place to store end of branch kmer

            branch[0].clear();
            branch[1].clear();

            int branch_offset = 0;
            uint64_t branch_color[2];


            // start of a bubble handling
            for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
                // follow each strand or supernode
                ssize_t edge = dbg.outgoing_edge(node_i, x);
                if (edge == -1)
                    continue;
                branch[branch_num] += base[x];
                // build color mask
                uint64_t color_mask = 0;
                for (int c = 0; c < num_colors; c++)
                    color_mask |= colors[edge * num_colors + c] << c;
                branch_color[branch_num] = color_mask;

                // walk along edges until we encounter 
                ssize_t node_pos = dbg._edge_to_node(edge);
                while (dbg.indegree(node_pos) == 1 && dbg.outdegree(node_pos) == 1) {
                    visited[node_pos] = 1;
                    ssize_t next_edge = 0;
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(node_pos, x2);
                        if (next_edge != -1) {
                            branch[branch_num] += base[x2];
                            break;
                        }
                    }
                    node_pos = dbg._edge_to_node(next_edge);
                    //cout << node_pos << ":" << dbg.node_label(node_pos) << "\n";
                }
                if (found_miss) {
                    std::cerr << "dbg.indegree(node_pos) = " << dbg.indegree(node_pos)  << " dbg.outdegree(node_pos) = " << dbg.outdegree(node_pos)  << std::endl;
                    ssize_t next_edge = 0;
                    std::cerr << "outgoing bases: ";
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(node_pos, x2);
                        uint64_t color_mask = 0;

                        if (next_edge != -1) {
                            for (int c = 0; c < num_colors; c++)
                                color_mask |= colors[next_edge * num_colors + c] << c;


                            std::cerr << base[x2] << " (c " << color_mask << ") (p " << dbg._edge_to_node(next_edge) << ")" << std::endl;
;
                        }
                    }
                    std::cerr << std::endl;
                    
                }
                // cout << "Stopped due to : " << dbg.indegree(node_pos) << ":" << dbg.outdegree(node_pos) << ":" << branch_offset << "\n";

                end[branch_num++] =  (dbg.indegree(node_pos) > 1) ? node_pos : 0;
                branch_offset = 0;
            }
            // check if both branches ended on the same kmer and they pass the requested color masks
            //cout << "Trying " << branch_color[0] << ":" << branch_color[1] << " " << end[0] << ":" << end[1] <<"\n";
            //cout << ref_color << ":" << sample_mask << "\n";
            //cout << "PutativeStart flank: " << dbg.node_label(start) << " c: " << branch_color[0] << ":" << branch_color[1] << "\n";
            if (found_miss) {
                std::cerr << "arm sizes: " << branch[0].size() << "  " << branch[1].size() << std::endl;
                std::cerr << branch[0] << std::endl << branch[1] << std::endl;
            }
                                    
            // check same end node
            if ((end[0] && end[0] == end[1]) ) {
                if (found_miss)  std::cerr << "Missing bubble passed end check" << std::endl;
                // check color:
                if (true || ((ref_color & branch_color[0] && !(~ref_color & branch_color[0]) &&
                  sample_mask & branch_color[1] && !(~sample_mask & branch_color[1])) || 
                 (ref_color & branch_color[1] && !(~ref_color & branch_color[1]) &&
                  sample_mask & branch_color[0] && !(~sample_mask & branch_color[0])))) {
                    cout << "\nStart flank: " << dbg.node_label(start_node) << " c: " << branch_color[0] << ":" << branch_color[1] << "\n";
                    cout << "Branch: " << branch[0] << "\n";
                    cout << "Branch: " << branch[1] << "\n";
                    cout << "End flank: " << dbg.node_label(end[0]) << "\n";
                    if (found_miss) std::cerr << "Reported 'missing' bubble" << std::endl;
                }
            }
        }
    }
    cerr << "Find bubbles time: " << getMilliSpan(t) << std::endl;
}

//FIXME : deal with uppercase/lowercase/unspecified NT in FASTA
unsigned parse_fasta(const std::string& ref_fasta_fname, std::vector<std::string>& ref_fastas, std::vector<std::string>& ids)
{

    std::ifstream f(ref_fasta_fname);
    std::string s;
    //process first line
    if (!std::getline(f, s) || !(s[0] == '>'||s[0] == '@')) {
        std::cerr << "ERROR: File " << ref_fasta_fname << " doesn't seem like a FASTA file." << std::endl;
        exit(EXIT_FAILURE);
    }

    bool capture_mode = true; // capture mode
    ids.push_back(s);
    s = "";
    ref_fastas.push_back(s);

    //process rest of lines
    while(std::getline(f, s) && s.size() > 0){


        if (s[0] == '>' || s[0] == '@') {
            ids.push_back(s);
            capture_mode = true;
            s = "";
            ref_fastas.push_back(s);
        } else if (s[0] == '+') {
            capture_mode = false;
        } else if (capture_mode == true) {
            *(ref_fastas.end()-1) += s;
        }
    }
    return ref_fastas.size();
        

}

//KMC discards reads with non DNA letters, so we do the same
bool read_okay(const std::string& readstr)
{
    for (const char& c : readstr) {
        if (c == 'A' || c == 'G' || c == 'C' || c == 'T') {
            continue;
        } else {
            return false;
        }
    }

    return true;
}
    
int main(int argc, char* argv[])
{
    global_t = getMilliCount();
    parameters_t p;
    parse_arguments(argc, argv, p);



    std::cerr << "=== Loading data structures from disk ==" << std::endl;
    debruijn_graph_shifted<> dbg;
    unsigned long long t = getMilliCount();
    std::cerr << "Loading succinct de Bruijn graph...";
    load_from_file(dbg, p.input_filename);
    std::cerr << "done" << std::endl << "graph load wall time: " << (getMilliCount() - t) / 1000.0 << " s" << std::endl;


    std::cerr << "=== Calculating/reporting stats ===" << std::endl;
    cerr << "k             : " << dbg.k << endl;
    cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << dbg.num_edges() << endl;
//  cerr << "colors        : " << colors.size() / dbg.num_edges() << endl; 
    cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
//  cerr << "Color size    : " << size_in_mega_bytes(colors) << " MB" << endl;
    std::cerr << "Total wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;




    uint64_t mask2 = (p.sample_mask.length() > 0) ? atoi(p.sample_mask.c_str()) : -1;
    std::vector<std::string> ref_fastas;
    std::vector<std::string> ids;

    std::cerr << "Loading reference FASTA file " << p.ref_fasta  << "...";
    std::cerr << "Loaded " << parse_fasta(p.ref_fasta, ref_fastas, ids) << " sequences" << std::endl;
    std::cerr << "Total wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;
    int i = 0;
    std::cerr << "Creating Rank support..." ;


    std::cerr << "done" << std::endl << "Total wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;

    std::cerr << "=== making sd_vector_builder ==" << std::endl;

  
    std::cerr << "=== determining size of color matrix ==" << std::endl;
    std::vector<std::future<std::set<unsigned long long> > > retvals;
    // make global version since async can't pass reference args
  

    gdbg = &dbg;

    size_t sd_base_pos = 0;
    unsigned count = 0;
    unsigned discarded = 0; // for having non DNA symbols
    unsigned long long totkmers = 0;
    for (std::vector<std::string >::iterator rf = ref_fastas.begin(); rf != ref_fastas.end(); ++i, ++rf) {
        if (count % 100000 == 0) std::cerr << "on read " << count << " discarded: " << discarded
                                           << " read stats progress: " << count << "/" << ref_fastas.size()
                                           << "(" << count / (float)ref_fastas.size() << "%)"<< std::endl;

        //KMC discards reads with non DNA letters, so we do the same
        if (read_okay(*rf)) {
            totkmers += rf->size() - dbg.k + 1;
            count += 1;
        } else {
            discarded += 1;
        }
    }
    std::cerr << "total kmers: " << totkmers << std::endl;
    size_t n = count /*reads*/ * dbg.num_edges()/*edges*/;  //params.n;

    size_t m = totkmers; /*total k-mers*/ //1611287698ull; //params.m;
    std::cerr << "stack allocing sdsl::sd_vector_builder base object with n=" << n
              << " m=" << m << std::endl;
    sdsl::sd_vector_builder b_builder(n, m);// = bit_vector(num_edges*num_color, 0);

    std::cerr << "===  traversing graph ==" << std::endl;
  
    discarded = 0;
    count = 0;
    std::vector<bool> okay_reads;
    unsigned long long num_okay_reads = 0;
    std::vector<bool> read_starts;
    unsigned long long num_read_starts = 0;
    for (std::vector<std::string >::iterator rf = ref_fastas.begin(); rf != ref_fastas.end(); ++i, ++rf) {
        if (count % 100000 == 0) std::cerr << "on read "
                                           << count << "/" << ref_fastas.size()
                                           << "(" << count / (float)ref_fastas.size() << "%)"
                                           << " discarded: " << discarded
                                           << " sd_vector progress: " << b_builder.items() << "/" << b_builder.capacity()
                                           << "(" << (float)b_builder.items() / (float)b_builder.capacity() << "%)"<< std::endl;
        //KMC discards reads with non DNA letters, so we do the same
        if (read_okay(*rf)) {
            okay_reads.push_back(true);
            num_okay_reads++;
            std::vector<std::set<unsigned long long> > found_kmers = walk_refs(*rf);
            unsigned subread=0;
            for(std::vector<std::set<unsigned long long> >::iterator subcolorit = found_kmers.begin(); subcolorit != found_kmers.end(); ++subcolorit) {
                if (subread == found_kmers.size() - 1) {
                    read_starts.push_back(true);
                    num_read_starts++;
                    
                } else {
                    read_starts.push_back(false);
                }
                for (std::set<unsigned long long>::iterator it = subcolorit->begin(); it != subcolorit->end(); ++it) {
                    b_builder.set(sd_base_pos + *it);
                }
                sd_base_pos += dbg.num_edges();
                subread +=1;
            }
            count++;

        } else {
            okay_reads.push_back(false);
            discarded += 1;
        }
    }

    sdsl::sd_vector_builder valid_read_builder(okay_reads.size(), num_okay_reads);// = bit_vector(num_edges*num_color, 0);


    for (unsigned long long i = 0; i < okay_reads.size(); ++i) {
        if (okay_reads[i]) {

            valid_read_builder.set(i);
        }
    }
    sdsl::sd_vector<> b_valid(valid_read_builder);
    std::string validoutfilename = p.output_matrix + ".valid"; //"97_GTGGCC_L008_reverse_paired.colors.sd_vector";
    sdsl::store_to_file(b_valid, validoutfilename);

    
    sdsl::sd_vector_builder read_start_builder(read_starts.size(), num_read_starts);// = bit_vector(num_edges*num_color, 0);
    for (unsigned long long i = 0; i < read_starts.size(); ++i) {
        if (read_starts[i]) {

            read_start_builder.set(i);
        }
    }
    sdsl::sd_vector<> b_start(read_start_builder);
    std::string startoutfilename = p.output_matrix + ".starts"; //"97_GTGGCC_L008_reverse_paired.colors.sd_vector";
    sdsl::store_to_file(b_start, startoutfilename);


    
    std::cerr << "sd_vector items: " << b_builder.items()  << " sd_vector capacity " << b_builder.capacity() << std::endl;
    std::cerr << "Total wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;
    sdsl::sd_vector<> b(b_builder);
 
    std::string outfilename = p.output_matrix; //"97_GTGGCC_L008_reverse_paired.colors.sd_vector";
    sdsl::store_to_file(b, outfilename);
    std::cerr << "done!" << std::endl;

}
