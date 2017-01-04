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
  
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();

  params.ref_color       = ref_color_arg.getValue();
  params.sample_mask     = sample_mask_arg.getValue();
  // params.ref_fasta       = ref_fasta_arg.getValue();
  //params.output_matrix   = output_matrix_arg.getValue();
}


    
int main(int argc, char* argv[])
{
    global_t = getMilliCount();
    parameters_t p;
    parse_arguments(argc, argv, p);
    

    // load the color matrix
    sd_vector<> colors;
    std::cerr << "=== Loading data structures from disk ==" << std::endl;
    load_from_file(colors, p.input_filename);
    std::cerr << "Elapsed wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;

    unsigned long long rows = atoll(p.ref_color.c_str());
    unsigned long long cols = atoll(p.sample_mask.c_str());
    std::cerr << "rows*cols = " << rows*cols << " sd_vector.size(): " << colors.size() << std::endl;

    rank_support_sd<1> color_ranks(&colors);
    select_support_sd<1> color_selects(&colors);

    unsigned long long num_ones = color_ranks(colors.size()-1);
    std::cerr << "num 1's: " << num_ones << std::endl;
    std::cerr << "first 1: " << color_selects(1) << std::endl;
    std::cerr << "last 1: " << color_selects(num_ones) << std::endl;

    // for each 1 in the matrix, compute the position in the new serialized matrix 'elems'
    std::vector<unsigned long long> elems;

    std::cerr << "=== gen int vector ==" << std::endl;
    for (unsigned long long i = 0; i < num_ones; ++i) {

        unsigned long long elemval = color_selects(i+1);

        unsigned long long old_row = elemval % rows;
        unsigned long long old_col = (elemval - old_row) / rows;
        unsigned long long newelem = old_row * cols + old_col;
        if (i < 10) std::cerr << i <<"th set bit is at " << elemval << "changing to " << newelem << std::endl;        
        elems.push_back(newelem) ;

    }
    std::cerr << "Elapsed wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;



    // sort the 1 positions in elems
    std::cerr << "=== sort int vector ==" << std::endl;
    std::sort(elems.begin(), elems.end());
    std::cerr << "Elapsed wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;

    // then build the new succinct matrix
    std::cerr << "=== alloc builder ==" << std::endl;
    sdsl::sd_vector_builder b_builder(colors.size(), num_ones);
    std::cerr << "Elapsed wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;


    std::cerr << "=== elias-fano encode to sd_vector from int vector ==" << std::endl;
    for (unsigned long long i = 0; i < num_ones; ++i) {
        if (i < 10) std::cerr << i <<"th set bit is at " << elems[i] << std::endl;
        b_builder.set(elems[i]);
    }

    std::cerr << "sd_vector items: " << b_builder.items()  << " sd_vector capacity " << b_builder.capacity() << std::endl;
    std::cerr << "Elapsed wall time: " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;    
    sdsl::sd_vector<> b(b_builder);
    select_support_sd<1> b_color_selects(&b);
    rank_support_sd<1> b_color_ranks(&b);

    
    unsigned long long b_num_ones = b_color_ranks(b.size()-1);
    std::cerr << "num 1's: " << b_num_ones << std::endl;
    std::cerr << "first 1: " << b_color_selects(1) << std::endl;
    std::cerr << "last 1: " << b_color_selects(b_num_ones) << std::endl;


    
    for (unsigned long long i = 0; i < num_ones; ++i) {

        unsigned long long elemval = b_color_selects(i+1);
        if (i < 10) std::cerr << i <<"th set bit is at " << elemval << std::endl;
    }
    sdsl::store_to_file(b, p.color_filename );

  


}
