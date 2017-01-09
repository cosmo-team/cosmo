
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
#include "pack-color.hpp"

int getMilliCount()
{
    timeb tb;
    ftime(&tb);
    int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
    return nCount;
}


int getMilliSpan(int nTimeStart)
{
    int nSpan = getMilliCount() - nTimeStart;
    if(nSpan < 0)
        nSpan += 0x100000 * 1000;
    return nSpan;
}

std::string file_extension = ".<extension>";


void parse_arguments(int argc, char **argv, parameters_t & params)
{
    TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
                                                             "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> num_colors_arg("num_colors",
                                                         "Number of colors", true, "", "num colors", cmd);
    TCLAP::UnlabeledValueArg<std::string> n_arg("n",
                                                         "sd_vector n=Vector size", true, "", "sd_vector n=Vector size", cmd);
    TCLAP::UnlabeledValueArg<std::string> m_arg("m",
                                                         "sd_vector m=The number of 1-bits", true, "", "sd_vector m=The number of 1-bits", cmd);
    std::string output_short_form = "output_prefix";
    TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
                                                   "Output prefix. Graph will be written to [" + output_short_form + "]" + file_extension + ". " +
                                                   "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
    cmd.parse( argc, argv );
    params.input_filename  = input_filename_arg.getValue();
    params.num_colors  = atoi(num_colors_arg.getValue().c_str());
    params.n  = atoll(n_arg.getValue().c_str());
    params.m  = atoll(m_arg.getValue().c_str());
    params.output_prefix   = output_prefix_arg.getValue();

}

void deserialize_color_bv(std::ifstream &colorfile, color_bv &value)
{
    colorfile.read((char *)&value, sizeof(color_bv));
}

int main(int argc, char * argv[])
{
    const bool rrr = false;
    std::cerr << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
    std::cerr <<"Starting" << std::endl;
    parameters_t params;
    parse_arguments(argc, argv, params);

    const char * file_name = params.input_filename.c_str();
    std::cerr << "file name: " << file_name << std::endl;

// Open File
    std::ifstream colorfile(file_name, std::ios::in|std::ios::binary);

    colorfile.seekg(0, colorfile.end);
    size_t end = colorfile.tellg();
    colorfile.seekg(0, colorfile.beg);
    std::cerr << "file size: " << end << std::endl;
    std::cerr << "sizeof(color_bv): " << sizeof(color_bv) << std::endl;
    size_t num_color = params.num_colors;
    size_t num_edges = end / sizeof(color_bv);
    //size_t n = 4000000000000;
    //size_t n = 201136208922 + 7967979875878; // 80 color set -d
    // size_t n = 1049958 + 54489174; // ecoli6 -d
    size_t n = params.n;
    //size_t m = 120543830541;
    // size_t m = 201136208922; // 80 color set -d
    //size_t m = 54489174; // ecoli6 -d
    size_t m = params.m;
    sdsl::bit_vector *b = NULL;
    sdsl::sd_vector_builder *b_builder = NULL;
    if (rrr) {
        std::cerr << "Using RRR encoding.  Allocating a vector for " << num_edges*num_color << " bits." << std::endl;
        file_extension = ".rrr";
        b = new sdsl::bit_vector(num_edges*num_color, 0);
    } else{
        std::cerr << "Using Elias-Fano encoding" << std::endl;
        std::cerr << "stack allocing sdsl::sd_vector_builder base object with n=" << n
                  << " m=" << m << std::endl;
        file_extension = ".sd_vector";
        b_builder = new sdsl::sd_vector_builder(n, m);// = bit_vector(num_edges*num_color, 0);
        std::cerr << "builder size: " << b_builder->size() << " capacity: " << b_builder->capacity() << std::endl;
    }
    std::cerr << "Succinct builder object allocated." << std::endl;

    size_t cnt0 = 0;
    size_t cnt1 = 0;
    for (size_t i=0; i < num_edges; i++) {
        if (i % 100000000 == 0) {
            std::cerr << "deserializing edge " << i
                      << " cnt0: " << cnt0
                      << " cnt1: " << cnt1 << std::endl;
        }
        color_bv value;
//        std::cerr << "about to deserialize" << std::endl;
        deserialize_color_bv(colorfile, value);
//        std::cerr << "done deserialize" << std::endl;
        for (size_t j=0; j < num_color; j++) {
            if (value[j]){
//b[i*num_color + j] = value[j];
//                std::cerr << "setting a bit at pos" <<  i * num_color + j<< std::endl;
                if (rrr) {
                    (*b)[i*num_color + j] = value[j];
                } else {
                    b_builder->set(i*num_color + j);
                }
                cnt1++;
            } else {
                if (rrr) {
                    (*b)[i*num_color + j] = value[j];
                }
                cnt0++;

            }
        }
    }
    std::cerr << "edges: " << num_edges << " colors: " << num_color << " Total: " << num_edges * num_color << std::endl;
    std::cerr << cnt0  << ":" << cnt1 << std::endl;

    int sysTime = getMilliCount();
/*
  bit_vector bv(b);
  sysTime = getMilliCount();
  std::cerr << "BV Creation Time: " << getMilliSpan(sysTime) << endl;
  for (size_t i=0; i < num_edges*num_color; i++) {
  bv[i];
  }
  std::cerr << "BV Access Time: " << getMilliSpan(sysTime) << endl;
  std::cerr << "BV Size (MB): " << size_in_mega_bytes(b) << endl;
*/
    if (rrr) {
        sdsl::rrr_vector<63> rrrb(*b);
        std::cout << "RRR Creation Time: " << getMilliSpan(sysTime) << std::endl;
        sysTime = getMilliCount();
        // for (size_t i=0; i < num_edges*num_color; i++) {
        //   rrrb[i];
        // }
        std::cout << "RRR AccessTime: " << getMilliSpan(sysTime) << std::endl;
        std::cout << "RRR Size (MB): " << size_in_mega_bytes(rrrb) << std::endl;
        char * base_name = basename(const_cast<char*>(params.input_filename.c_str()));
        std::string outfilename = ((params.output_prefix == "")? base_name : params.output_prefix) + file_extension;
        store_to_file(rrrb, outfilename);
        delete b;
    } else {
    
        sysTime = getMilliCount();
        sdsl::sd_vector<> b(*b_builder);

// rrr_vector<63> rrrb(b);
// std::cerr << "RRR Creation Time: " << getMilliSpan(sysTime) << endl;
// sysTime = getMilliCount();
// for (size_t i=0; i < num_edges*num_color; i++) {
// rrrb[i];
// }
// std::cerr << "RRR AccessTime: " << getMilliSpan(sysTime) << endl;
// std::cerr << "RRR Size (MB): " << size_in_mega_bytes(rrrb) << endl;
        char * base_name = basename(const_cast<char*>(params.input_filename.c_str()));
        std::string outfilename = ((params.output_prefix == "")? base_name : params.output_prefix) + file_extension;
        sdsl::store_to_file(b, outfilename);
        delete b_builder;
    }
/*
  sysTime = getMilliCount();
  sd_vector<> sdb(b);
  std::cerr << "SD Creation Time: " << getMilliSpan(sysTime) << endl;
  sysTime = getMilliCount();
  for (size_t i=0; i < num_edges*num_color; i++) {
  sdb[i];
  }
  std::cerr << "SD Access Time: " << getMilliSpan(sysTime) << endl;
  std::cerr << "SD Size (MB): " << size_in_mega_bytes(sdb) << endl;

  sysTime = getMilliCount();
  hyb_vector<> hyb(b);
  std::cerr << "Hyb Creation Time: " << getMilliSpan(sysTime) << endl;
  sysTime = getMilliCount();
  for (size_t i=0; i < num_edges*num_color; i++) {
  hyb[i];
  }
  std::cerr << "Hyb Access Time: " << getMilliSpan(sysTime) << endl;
  std::cerr << "Hyb Size (MB): " << size_in_mega_bytes(hyb) << endl;
*/
}
