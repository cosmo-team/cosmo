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

debruijn_graph_shifted<>* gdbg;
//sd_vector<>* gcolors;
//rank_support_sd<>* gcolor_ranks;


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


    
// Follow 
inline bool advance(debruijn_graph_shifted<>& dbg, const std::string& read_seq, ssize_t& node_i, ssize_t& node_i_pos,   std::vector<std::set<unsigned long long> >& found_kmers)
{


    unsigned node_label_size = dbg.k - 1;  

    // Follow the read_seq by 1bp to get the next edge
    ssize_t edge = dbg.outgoing_edge(node_i, dna_ord(read_seq[node_i_pos + node_label_size]));
    if (edge == -1) {
        std::cerr <<  read_seq.substr(node_i_pos, node_label_size) << " Reference fasta guides graph traversal through invalid path at position "
                  << node_i_pos + node_label_size // we want the right end of the current position k-mer
                  << std::endl ;

        // FIXME: In one dataset, this error condition occurred.  Because we were hoping to make RECOMB at the time, I just had it ignore the problem read using the 'return false' below
        // I've reinstated the bail-out here because it may cause much greater headaches later if not accounted for correctly
        // It might still be okay to just ignore the problem read as long as it doesn't throw off counting read pairs, etc.
        assert(false);
        exit(1);
        //return false;
    }
    
    // if this k-mer has already been seen in the current set
    if ((found_kmers.end()-1)->find(edge) != (found_kmers.end()-1)->end()) {
        // create a new set and add it to the vector
        std::set<unsigned long long> newset;
        found_kmers.push_back(newset);
    }

    // Add the kmer to the current set
    (found_kmers.end()-1)->insert(edge);        

    // now look up the next node
    node_i = dbg._edge_to_node(edge);
    node_i_pos++;
    return true;
    
    
}

// This function takes a string (the read) and returns a vector of sets.
// If the read has repeated k-mers within it, we want to be able to
// recover the full walk of the read in the de Bruijn graph later

// To do this, we generate a /vector/ of k-mer sets.
// Each k-mer set comprises all the k-mers in a non-repeating section of the read
// The first repeated k-mer will be found in the second set of the vector
std::vector<std::set<unsigned long long> > walk_refs(std::string read_seq )
{
    std::vector<std::set<unsigned long long> > found_kmers;
    std::set<unsigned long long> initial;
    found_kmers.push_back(initial);
        

    // First, find the edge corresponding to the first k nucleotides of the read
    unsigned node_label_size = gdbg->k - 1;
    std::string first_edge_label(read_seq.substr(0, node_label_size + 1));
    auto edge = gdbg->index(first_edge_label.begin());
    assert(edge);

    ssize_t first_node =  gdbg->_edge_to_node(get<0>(*edge));

    // In order to follow a read sequence through the graph, we need to keep track of positions
    // of both where we are in the read and where we are in the graph:
    ssize_t node_i = first_node; // node_i is a cdbg node labeled with a k-mer existing in the read_seq
    ssize_t node_i_pos = 0;  // starting position in the reference sequence for the above k-mer


    
    while(node_i_pos < (ssize_t)read_seq.size() - gdbg->k + 1 ) {
        if (!advance(*gdbg, read_seq, node_i, node_i_pos, found_kmers)) break;
    }

    return found_kmers;

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
// this function checks to see if a read is okay (i.e. only contains uppercase nucleotides)
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
    cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
    std::cerr << "Total wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;





    std::vector<std::string> ref_fastas;
    std::vector<std::string> ids;

    std::cerr << "Loading reference FASTA file " << p.ref_fasta  << "...";
    std::cerr << "Loaded " << parse_fasta(p.ref_fasta, ref_fastas, ids) << " sequences" << std::endl;
    std::cerr << "Total wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;
    int i = 0;

    std::cerr << "=== determining size of color matrix ==" << std::endl;

    // make the global pointer point to the automatic dbg instance
    // FIXME: just pass dbg by reference where needed, the global pointer was for legacy reasons
    gdbg = &dbg;

    // We fill the color matrix in column major order
    // sd_base_pos marks the index in the serialized matrix that starts a column
    // so each '1' in the serialized matrix will be filled in increasing order
    // of the base (for that column) plus an offset (for that row)
    
    size_t sd_base_pos = 0; 
    unsigned count = 0; // number of reads processes
    unsigned discarded = 0; // number of reads having non DNA symbols
    unsigned long long tot_kmers = 0; // total k-mers, not total unique k-mers

    // sdsl::sd_vector needs to know parameters n and m before any bits are added
    // we know how many distinct k-mer values there are (rows) but not how many subreads (columns) there are
    // So we scan through all the reads first to determine how big and how dense the matrix needs to be
    for (std::vector<std::string >::iterator rf = ref_fastas.begin(); rf != ref_fastas.end(); ++i, ++rf) {
        if (count % 100000 == 0) std::cerr << "on read " << count << " discarded: " << discarded
                                           << " read stats progress: " << count << "/" << ref_fastas.size()
                                           << "(" << count / (float)ref_fastas.size() << "%)"<< std::endl;

        //KMC discards reads with non DNA letters, so we do the same
        if (read_okay(*rf)) {
            tot_kmers += rf->size() - dbg.k + 1;
            count += 1;
        } else {
            discarded += 1;
        }
    }
    std::cerr << "total kmers: " << tot_kmers << std::endl;

    // now compute the actual size of the sd_vector
    size_t n = count /*reads*/ * dbg.num_edges()/*edges*/;  // FIXME: this really should be the number of subreads, not reads
    size_t m = tot_kmers; 
    std::cerr << "stack alloc'ing sdsl::sd_vector_builder base object with n=" << n
              << " m=" << m << std::endl;
    sdsl::sd_vector_builder b_builder(n, m);

    std::cerr << "===  traversing graph ==" << std::endl;

    
    // reset the counters from the initial scan
    discarded = 0;
    count = 0;

    // we keep track of which reads are valid and which subread (column) is the start of a read
    // you can think of these as a header above the matrix
    // these will eventually be turned into sd_vectors too, but they are relatively small, and we don't know how big they will be
    // FIXME: it would be smarter to just get rid of all the invalid reads (and their mates) before running this.
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
        // KMC2 discards reads with non DNA letters, so we do the same
        if (read_okay(*rf)) {
            okay_reads.push_back(true);
            num_okay_reads++;

            // walk the graph using a read
            std::vector<std::set<unsigned long long> > found_kmers = walk_refs(*rf);

            // for each subread returned from the walk, populate the header and the column
            unsigned subread = 0;
            for(std::vector<std::set<unsigned long long> >::iterator subreadit = found_kmers.begin(); subreadit != found_kmers.end(); ++subreadit) {
                if (subread == found_kmers.size() - 1) {
                    read_starts.push_back(true);
                    num_read_starts++;
                    
                } else {
                    read_starts.push_back(false);
                }
                for (std::set<unsigned long long>::iterator it = subreadit->begin(); it != subreadit->end(); ++it) {
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

    // now that we've completed the traversal for all reads, we can generate the succinct headers
    // first do the 'valid' vector
    sdsl::sd_vector_builder valid_read_builder(okay_reads.size(), num_okay_reads);// = bit_vector(num_edges*num_color, 0);


    for (unsigned long long i = 0; i < okay_reads.size(); ++i) {
        if (okay_reads[i]) {

            valid_read_builder.set(i);
        }
    }
    sdsl::sd_vector<> b_valid(valid_read_builder);
    std::string validoutfilename = p.output_matrix + ".valid"; //"97_GTGGCC_L008_reverse_paired.colors.sd_vector";
    sdsl::store_to_file(b_valid, validoutfilename);

    // then do the 'start' vector
    sdsl::sd_vector_builder read_start_builder(read_starts.size(), num_read_starts);// = bit_vector(num_edges*num_color, 0);
    for (unsigned long long i = 0; i < read_starts.size(); ++i) {
        if (read_starts[i]) {

            read_start_builder.set(i);
        }
    }
    sdsl::sd_vector<> b_start(read_start_builder);
    std::string startoutfilename = p.output_matrix + ".starts"; //"97_GTGGCC_L008_reverse_paired.colors.sd_vector";
    sdsl::store_to_file(b_start, startoutfilename);


    // finally, finish off the main color matrix
    std::cerr << "sd_vector items: " << b_builder.items()  << " sd_vector capacity " << b_builder.capacity() << std::endl;
    std::cerr << "Total wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;
    sdsl::sd_vector<> b(b_builder);
 
    std::string outfilename = p.output_matrix; //"97_GTGGCC_L008_reverse_paired.colors.sd_vector";
    sdsl::store_to_file(b, outfilename);
    std::cerr << "done!" << std::endl;

}
