// this file is for merging two uncompressed color matrices

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
#include "color-merge.hpp"

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
    TCLAP::UnlabeledValueArg<std::string> plan_filename_arg("plan",
                                                             "Plan file. Currently only supports DSK's binary format (for k<=64).", true, "", "plan_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> matrix1_filename_arg("matrix1",
                                                             "Matrix1 file. Currently only supports DSK's binary format (for k<=64).", true, "", "matrix1_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> matrix2_filename_arg("matrix2",
                                                             "Matrix2 file. Currently only supports DSK's binary format (for k<=64).", true, "", "matrix2_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> num_colors1_arg("num_colors1",
                                                         "Number of colors 1", true, "", "num colors 1", cmd);
    TCLAP::UnlabeledValueArg<std::string> num_colors2_arg("num_colors2",
                                                         "Number of colors 2", true, "", "num colors 2", cmd);
    std::string output_short_form = "output_prefix";
    TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
                                                   "Output prefix. Graph will be written to [" + output_short_form + "]" + file_extension + ". " +
                                                   "Default prefix: basename(plan_file).", false, "", output_short_form, cmd);
    cmd.parse( argc, argv );
    params.plan_filename  = plan_filename_arg.getValue();
    params.matrix1_filename  = matrix1_filename_arg.getValue();
    params.matrix2_filename  = matrix2_filename_arg.getValue();    
    params.num_colors1  = atoi(num_colors1_arg.getValue().c_str());
    params.num_colors2  = atoi(num_colors2_arg.getValue().c_str());    
    params.output_prefix   = output_prefix_arg.getValue();

}

void deserialize_color_bv(std::ifstream &colorfile, color_bv &value)
{
    colorfile.read((char *)&value, sizeof(color_bv));
}

int main(int argc, char * argv[])
{
    const bool rrr = false;
    parameters_t params;
    parse_arguments(argc, argv, params);

    const char * plan_file = params.plan_filename.c_str();
    std::ifstream planfile(plan_file, std::ios::in|std::ios::binary);

    planfile.seekg(0, planfile.end);
    size_t end = planfile.tellg();
    planfile.seekg(0, planfile.beg);
    
    const char * matrix1_file = params.matrix1_filename.c_str();
    sdsl::sd_vector<> colors1;
    load_from_file(colors1, params.matrix1_filename);




    const char * matrix2_file = params.matrix2_filename.c_str();
    sdsl::sd_vector<> colors2;
    load_from_file(colors2, params.matrix2_filename);

    size_t out_colors = params.num_colors1 + params.num_colors2;  
    size_t n = out_colors * end;// (colors1.size() + colors2.size()) * 1.5; // total bits FIXME: we can figure this out
    size_t m = colors1.low.size() + colors2.low.size();
  


    sdsl::sd_vector_builder *b_builder = NULL;
    

    file_extension = ".sd_vector";
    b_builder = new sdsl::sd_vector_builder(n, m);// = bit_vector(num_edges*num_color, 0);
    std::cerr << "builder size: " << b_builder->size() << " capacity: " << b_builder->capacity() << std::endl;

    SDIter color1iter(&colors1);
    SDIter color2iter(&colors2);

    size_t color1_row = 0;
    size_t color2_row = 0;
    size_t output_row = 0;
    
    char planstep = 0;
    while (planfile >> planstep) {
        size_t next_pos = 0;        
        if (planstep == 3) {
            next_pos = color1iter.peek();
            while (next_pos < (color1_row + 1) * params.num_colors1 && next_pos != -1) {
                color1iter.advance();
                size_t orig_col = next_pos - color1_row * params.num_colors1;
                b_builder->set(output_row * out_colors + orig_col);
                next_pos = color1iter.peek();
            }
            color1_row++;

            
            next_pos = color2iter.peek();
            while (next_pos < (color2_row + 1) * params.num_colors2 && next_pos != -1) {
                color2iter.advance();
                size_t orig_col = next_pos - color2_row * params.num_colors2;
                b_builder->set(output_row * out_colors + params.num_colors1 + orig_col);
                next_pos = color2iter.peek();
            }
            color2_row++;

            output_row++;
        } else if (planstep == 1) {
            next_pos = color1iter.peek();
            while (next_pos < (color1_row + 1) * params.num_colors1 && next_pos != -1) {
                color1iter.advance();
                size_t orig_col = next_pos - color1_row * params.num_colors1;
                b_builder->set(output_row * out_colors + orig_col);
                next_pos = color1iter.peek();
            }
            color1_row++;

            output_row++;
        } else if (planstep == 2) {
            next_pos = color2iter.peek();
            while (next_pos < (color2_row + 1) * params.num_colors2 && next_pos != -1) {
                color2iter.advance();
                size_t orig_col = next_pos - color2_row * params.num_colors2;
                b_builder->set(output_row * out_colors + params.num_colors1 + orig_col);
                next_pos = color2iter.peek();
            }
            color2_row++;

            output_row++;
        }
    }
        
    int sysTime = getMilliCount();
    

        sdsl::sd_vector<> b(*b_builder);

// rrr_vector<63> rrrb(b);
// std::cerr << "RRR Creation Time: " << getMilliSpan(sysTime) << endl;
// sysTime = getMilliCount();
// for (size_t i=0; i < num_edges*num_color; i++) {
// rrrb[i];
// }
// std::cerr << "RRR AccessTime: " << getMilliSpan(sysTime) << endl;
// std::cerr << "RRR Size (MB): " << size_in_mega_bytes(rrrb) << endl;
        const char * base_name = "merged_colors"; //basename(const_cast<char*>(params.input_filename.c_str()));
        std::string outfilename = ((params.output_prefix == "")? base_name : params.output_prefix) + file_extension;
        sdsl::store_to_file(b, outfilename);
        delete b_builder;
    
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
