//#include <iterator>
//#include <typeinfo>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
//#include <bitset>
//#include <set>

#include <stxxl.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

// #include <boost/log/trivial.hpp>
// #include <boost/log/core.hpp>
// #include <boost/log/expressions.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>

// TCLAP
#include "tclap/CmdLine.h"

#include "config.hpp"
#include "debug.hpp"
#include "utility.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"
#include "debruijn_graph.hpp"
#include "debruijn_graph_shifted.hpp"
#include "kmc_api/kmc_file.h"



struct parameters_t {
    std::string input_filename = "";
    std::string output_prefix  = "";
    std::string output_base    = "";
    size_t k = 0;
    size_t m = 0;
    bool swap = false;
    bool variable_order = false;
    bool shift_dummies  = false;
};

parameters_t parse_arguments(int argc, char **argv) {
    parameters_t params;
    TCLAP::CmdLine cmd(banner, ' ', version);
    TCLAP::ValueArg<size_t> kmer_length_arg("k", "kmer_length", "Length of edges (node is k-1). Needed for raw/DSK input.", false, 0, "length", cmd);
    TCLAP::ValueArg<size_t> mem_size_arg("m", "mem_size", "Internal memory to use (MB).", false, default_mem_size, "mem_size", cmd);
    TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix", "Output prefix.", false, "", "output_prefix", cmd);
    TCLAP::SwitchArg varord_arg("v", "variable_order", "Output .lcs file for variable order support.", cmd, false);
    TCLAP::SwitchArg shift_arg("d", "shift_dummies", "Shift all incoming dummies (slower but compresses better, and necessary for variable order without losing information).", cmd, false);
    // TODO: make this detect if directory by reading it
    //TCLAP::SwitchArg mono_arg("", "monochrome", "If multiple files are specified, don't add color data.");
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", "Input file.", true, "", "input_file", cmd);
    //TCLAP::UnlabeledMultiArg<string> input_filenames_arg("input", "file names", true, "", "input_file", cmd);
    //cmd.add( input_filename_arg );
    // TODO: add option for no reverse complements (e.g. if using bcalm input)
    // TODO: ASCII outputter: -a for last symbol, -aa for whole kmer
    // TODO: add DSK count parser (-c)
    // TODO: add option to keep temporaries (e.g. for iterative dbg)
    // TODO: add option for mutliple files being merged for colour
    cmd.parse( argc, argv );
    params.input_filename  = input_filename_arg.getValue();
    //params.temp_dir        = "";
    params.k               = kmer_length_arg.getValue();
    params.m               = mem_size_arg.getValue() * mb_to_bytes;
    params.variable_order  = varord_arg.getValue();
    params.shift_dummies   = shift_arg.getValue();
    params.output_prefix   = output_prefix_arg.getValue();
    return params;
}


int main(int argc, char* argv[])
{
    auto params = parse_arguments(argc, argv);
    std::string file_name = params.input_filename;
    params.output_base = boost::filesystem::path(file_name).stem().string();


    std::vector<CKMCFile*> kmer_data_bases;
    COSMO_LOG(trace) << "Reading KMC2 database list file..." << std::endl;
    size_t num_colors;

    size_t min_union;
    size_t max_union;
    uint32_t kmc_k;

    if ( !kmc_read_header(file_name, kmc_k, min_union, max_union, num_colors, kmer_data_bases) ) {
        COSMO_LOG(error) << "Error reading databases listed in KMC2 list file" << file_name;
        exit(EXIT_FAILURE);
    }
    size_t k = kmc_k;

    ofstream cfs;
    cfs.open(file_name + ".kmers", ios::out | ios::binary);

    /*! \param n Vector size.
     *  \param m The number of 1-bits. */
    size_t n = max_union * num_colors /* maybe find a better predictor of how many k-mers there will be */;
    size_t m = max_union;
    sdsl::sd_vector_builder b_builder(n, m);

    std::cerr << "Using Elias-Fano encoding" << std::endl;
    std::cerr << "stack allocing sdsl::sd_vector_builder base object with n=" << n
              << " m=" << m << std::endl;
    std::string file_extension = ".sd_vector";


    size_t builder_base = 0;
    int num_set = 0;
    size_t num_kmers_read = kmc_read_kmers(kmer_data_bases, k, [&](auto x, auto c) {
            cfs.write((char*)&x, sizeof(x)); //builder.push(x,c);
            for (unsigned int i = 0; i < num_colors; ++i) {
                if (c[i]) {
                    num_set++;
                    //if (builder_base + i >= m) std::cerr << "uh oh, trying to set too many 1's in sd_vector" << std::endl;
                    b_builder.set(builder_base + i);
                    
                }
            }
            builder_base += num_colors;
        });

    std::cerr << "m: " << m << " n: " << n << " num_set: " << num_set << std::endl;
    std::cerr << " sd_vector progress: " << b_builder.items() << "/" << b_builder.capacity()
              << "(" << (float)b_builder.items() / (float)b_builder.capacity() << "%)" << std::endl << std::flush;

    sdsl::sd_vector<> b(b_builder);

    
    sdsl::store_to_file(b, file_name + ".mergecolors");


    COSMO_LOG(info) << "Percentage of min union : " << num_kmers_read/(double)min_union * 100 << "%";
    COSMO_LOG(info) << "Percentage of max union : " << num_kmers_read/(double)max_union * 100 << "%";

    cfs.close();
    
}
