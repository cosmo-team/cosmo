#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>
// TCLAP
#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>

// C STDLIB Headers
#include <cstdio>
#include <cstdlib>

#include <libgen.h>

// Custom Headers
#include "uint128_t.hpp"
#include "debug.h"

using namespace std;
using namespace sdsl;

#include <cstdlib>
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

typedef struct p
{
  std::string input_filename = "";
  int num_colors;
} parameters_t;

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> num_colors_arg("num_colors",
            "Number of colors", true, "", "num colors", cmd);
  string output_short_form = "output_prefix";
  cmd.parse( argc, argv );
  params.input_filename  = input_filename_arg.getValue();
  params.num_colors  = atoi(num_colors_arg.getValue().c_str());
}

int main(int argc, char * argv[]) {
  cout <<"Starting\n";
  parameters_t params;
  parse_arguments(argc, argv, params);

  const char * file_name = params.input_filename.c_str();
  cout << file_name << "\n";

  // Open File
  ifstream colorfile(file_name, ios::in|ios::binary);

  colorfile.seekg(0, colorfile.end);
  size_t end = colorfile.tellg();
  colorfile.seekg(0, colorfile.beg);

  size_t num_color = params.num_colors;
  size_t num_edges = end / 8;
  bit_vector b = bit_vector(num_edges*num_color, 0);
  size_t cnt0 = 0;
  size_t cnt1 = 0;
  for (size_t i=0; i < num_edges; i++) {
    uint64_t value;
    colorfile.read((char *)&value, sizeof(uint64_t));
    for (size_t j=0; j < num_color; j++) {
      b[i*num_color + j] = (value & 1 << j) ? 0 : 1;
      if (b[i*num_color + j] == 0)
	cnt0++;
      else
	cnt1++;
    }
  }
  cout << "edges: " << num_edges << " colors: " << num_color << " Total: " << num_edges * num_color << endl;
  cout << cnt0  << ":" << cnt1 << endl;

  int sysTime = getMilliCount();
  bit_vector bv(b);
  sysTime = getMilliCount();
  cout << "BV Creation Time: " << getMilliSpan(sysTime) << endl;
  for (size_t i=0; i < num_edges; i++) {
    for (size_t j=0; j < num_color; j++) {
      bv[i + j];
    }
  }
  cout << "BV Access Time: " << getMilliSpan(sysTime) << endl;
  cout << "BV Size (MB): " << size_in_mega_bytes(b) << endl;

  sysTime = getMilliCount();
  rrr_vector<63> rrrb(b);
  cout << "RRR Creation Time: " << getMilliSpan(sysTime) << endl;
  sysTime = getMilliCount();
  for (size_t i=0; i < num_edges; i++) {
    for (size_t j=0; j < num_color; j++) {
      rrrb[i + j];
    }
  }
  cout << "RRR AccessTime: " << getMilliSpan(sysTime) << endl;
  cout << "RRR Size (MB): " << size_in_mega_bytes(rrrb) << endl;

  sysTime = getMilliCount();
  sd_vector<> sdb(b);
  cout << "SD Creation Time: " << getMilliSpan(sysTime) << endl;
  sysTime = getMilliCount();
  for (size_t i=0; i < num_edges; i++) {
    for (size_t j=0; j < num_color; j++) {
      sdb[i + j];
    }
  }
  cout << "SD Access Time: " << getMilliSpan(sysTime) << endl;
  cout << "SD Size (MB): " << size_in_mega_bytes(sdb) << endl;
}
