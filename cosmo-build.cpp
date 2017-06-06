//#include <iterator>
//#include <typeinfo>
#include <sys/timeb.h>
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
void serialize_color_bv(std::ofstream &cfs, const color_bv &color)//std::vector<color_bv>::iterator &colors, uint64_t index)
{
    cfs.write((char *)&color, sizeof(color_bv)); //FIXME: Is this the right way to serailize std::bitset?
}

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

int main(int argc, char* argv[]) {
    int merge_done_time = 0;
    int start_time = getMilliCount();
  using namespace boost::adaptors;
  COSMO_LOG(info) << "cosmo-build compiled with supported colors=" << NUM_COLS << std::endl;
  // Parameter extraction
  auto params = parse_arguments(argc, argv);
  std::string file_name = params.input_filename;
  params.output_base = boost::filesystem::path(file_name).stem().string();
  size_t k = params.k;
  
  stxxl::stats * stats = stxxl::stats::get_instance();
  stxxl::stats_data stats_begin(*stats);
  stxxl::block_manager * bm = stxxl::block_manager::get_instance();

  // Set logging level
  // TODO: make verbosity parameter
  // boost::log::core::get()->set_filter (
  //   boost::log::trivial::severity != boost::log::trivial::debug
  //);

  // TEST FILE
  // TODO: add fastq/fasta input support, and read from stdin if no file
  if (!boost::filesystem::exists(file_name)) {
    COSMO_LOG(error) << "Trouble opening " << file_name;
    exit(1);
  }
  input_format fmt;
  if (probably_list_of_files(file_name, &fmt)) {
    auto fmt_s = input_format_strings[(int)fmt];
    COSMO_LOG(info) << file_name << " looks like a list of " << fmt_s << " files.";
    if (fmt == input_format::dsk) {
      COSMO_LOG(error) << "Lists of DSK files aren't yet supported (we rely on KMC2's sorted output).";
      exit(1);
    }
  } else {
    fmt = get_format(file_name);
    auto fmt_s = input_format_strings[(int)fmt];
    COSMO_LOG(info) << file_name << " looks like a " << fmt_s << " file.";
    // TODO: support multiple input files, add to vector
    if (fmt == input_format::kmc) {
      COSMO_LOG(error) << "Single KMC files aren't yet supported (we're lazy). But you can put it in a list file.";
      exit(1);
    }
  }
  params.swap = (fmt == input_format::dsk);
  if ((fmt == input_format::dsk || fmt == input_format::raw) && k == 0) {
    COSMO_LOG(error) << "Need to specify k value when dealing with raw or DSK input files.";
    exit(EXIT_FAILURE);
  }

  // If its a named pipe, we need to make a new output name
  bool is_regular_file = boost::filesystem::is_regular_file(file_name);
  if (!is_regular_file) {
    namespace pt = boost::posix_time;
    auto now = pt::second_clock::local_time();
    params.output_base = pt::to_iso_string(now);
    COSMO_LOG(info) << file_name << " is not a regular file. Using \"" << params.output_base << "\" as the base name.";
  }

  // Check k value supported
  if (k > max_k) {
    COSMO_LOG(error) << "This version only supports k <= " << max_k << ". Try recompiling.";
    exit(1);
  }

  if (fmt == input_format::dsk || fmt == input_format::raw) {
    // TODO: refactor this pleasssse
    if (params.shift_dummies) {
      typedef debruijn_graph_shifted<> dbg_t;
      //typedef debruijn_graph<> dbg_t;
      // stxxl::syscall_file out_file(file_name + ".boss", stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
      typedef dbg_builder<dbg_t, kmer_t> builder_t;
      typedef typename builder_t::record_vector_t record_vector_t;

      stxxl::syscall_file in_file(file_name, stxxl::file::DIRECT | stxxl::file::RDONLY);
      record_vector_t in_vec(&in_file);
      typename record_vector_t::bufreader_type reader(in_vec);
      builder_t builder(params);

      COSMO_LOG(trace) << "Reading input and creating runs...";
      for (auto & x : reader) {
        builder.push(x);
      }
      // TODO: move dbg decl out of block and create copy constructor/operator etc
      //auto dbg = builder.build();
      vector<string> flags({"-", " "});
      string temp_lcs_file = params.output_prefix + params.output_base + ".lcs.temp";
      stxxl::syscall_file lcs_file(temp_lcs_file, stxxl::file::DIRECT | stxxl::file::RDWR /*WRONLY*/ | stxxl::file::CREAT | stxxl::file::TRUNC);
      stxxl::vector<uint8_t> * lcs_v;
      typename stxxl::vector<uint8_t>::bufwriter_type * lcs_writer;
      size_t idx = 0;
      {
        auto dbg = builder.build([&](size_t num_rows) {
          //lcs_v = sdsl::int_vector<8>((params.variable_order)*num_rows);
          lcs_file.set_size((params.variable_order)*num_rows);
          // a bit hacky. Should seperate finding dummies from merging them in builder,
          // so I can get the number of rows as a function call
          lcs_v = new stxxl::vector<uint8_t>(&lcs_file);
          lcs_writer = new stxxl::vector<uint8_t>::bufwriter_type(*lcs_v);
          //lcs_v.resize((params.variable_order)*num_rows);
        }, [&](auto x){
          size_t l;
          if (params.variable_order) {
            l = x.lcs;
            //lcs_v[idx++] = l;
            idx++;
            *lcs_writer << (uint8_t)l;
          }
          // Comment left in for potential verbose mode
          /*
          auto kmer = x.edge;
          string flag = flags[x.is_first_suffix];
          cerr << (int) x.is_first_prefix << " ";
          if (x.tag == in_dummy) {
            cerr << kmer_to_string(kmer, k, x.k) << flag << " ";
          } else if (x.tag == out_dummy) {
            assert(x.is_first_prefix);
            cerr << kmer_to_string(kmer<<2, k-1) << "$ ";
          } else {
            cerr << kmer_to_string(kmer, k) << flag << " ";
          }
          cerr << l << endl;
          */
        });
        //COSMO_LOG(info) << "size of DBG: " << size_in_mega_bytes(dbg) << " MB";
        sdsl::store_to_file(dbg, params.output_prefix + params.output_base + ".dbg");
      }

      delete lcs_writer;
      delete lcs_v;

      if (params.variable_order) {
        sdsl::wt_int<sdsl::rrr_vector<63>> lcs_wt;
        // TODO: make SDSL accept a templated type for forward iteration
        sdsl::construct(lcs_wt, temp_lcs_file, 1);
        //COSMO_LOG(info) << "size of LCS WT: " << size_in_mega_bytes(lcs_wt) << " MB";
        sdsl::store_to_file(lcs_wt, params.output_prefix + params.output_base + ".lcs");
        //lcs_v->clear();
        lcs_file.close_remove();
      }
    } else {
      // stxxl::syscall_file out_file(file_name + ".boss", stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
      typedef debruijn_graph<> dbg_t;
      typedef dbg_builder<dbg_t, kmer_t> builder_t;
      typedef typename builder_t::record_vector_t record_vector_t;

      stxxl::syscall_file in_file(file_name, stxxl::file::DIRECT | stxxl::file::RDONLY);
      record_vector_t in_vec(&in_file);
      typename record_vector_t::bufreader_type reader(in_vec);
      builder_t builder(params);

      COSMO_LOG(trace) << "Reading input and creating runs...";
      for (auto & x : reader) {
        builder.push(x);
      }

      auto dbg = builder.build();
      //COSMO_LOG(info) << "size of DBG: " << size_in_mega_bytes(dbg) << " MB";
      sdsl::store_to_file(dbg, params.output_prefix + params.output_base + ".dbg");
    }
  }
  else if (fmt == input_format::kmc) {
    typedef debruijn_graph_shifted<> dbg_t;
    typedef dbg_builder<dbg_t, kmer_t, color_bv> builder_t;
    if (!params.shift_dummies) {
      COSMO_LOG(error) << "Must use -d with KMC2 flow.";
      exit(EXIT_FAILURE);
    }
        
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
    k = params.k = kmc_k;

    // Check k value supported
    if (k > max_k) {
      COSMO_LOG(error) << "This version only supports k <= " << max_k << ". Try recompiling.";
      exit(EXIT_FAILURE);
    }

    if (num_colors > NUM_COLS) {
      COSMO_LOG(error) << "KMC file " << file_name << " contains " << num_colors << " colors which exceeds the compile time limit of "
                       << NUM_COLS << ". Please recompile with colors=" << NUM_COLS << " (or larger).";
      exit(EXIT_FAILURE);
    }

    builder_t builder(params);

    size_t num_kmers_read = kmc_read_kmers(kmer_data_bases, k, [&](auto x, auto c) {
      builder.push(x,c);
    });

    COSMO_LOG(info) << "Percentage of min union : " << num_kmers_read/(double)min_union * 100 << "%";
    COSMO_LOG(info) << "Percentage of max union : " << num_kmers_read/(double)max_union * 100 << "%";
    merge_done_time = getMilliSpan(start_time);
    COSMO_LOG(info) << "Color merge time:" << merge_done_time << std::endl;


    ofstream cfs;
    cfs.open(file_name + ".colors", ios::out | ios::binary);
    
    // TODO: improve the memory use for colors here (can SDSL use external vectors?)
//    bit_vector color_bv;
    size_t num_set = 0;
    size_t edge_idx = 0;
    auto dbg = builder.build([&](auto x) { // Pre merge
            // color_bv = bit_vector(x * num_colors);
    },[&](auto x) { // Merge visitor
      //auto kmer = x.edge;
      //cerr << kmer_to_string(kmer,k) << endl;
      auto color = get<0>(x.payload);
      serialize_color_bv(cfs, color);
      for (size_t color_idx = 0; color_idx < num_colors; color_idx++) {
        // TODO: try complemented bits in color-major form in a *sd_vector* for large data set.
        //color_bv[color_idx * num_colors + edge_idx] = !color[color_idx];
//        color_bv[edge_idx * num_colors + color_idx] = color[color_idx];
        num_set += color[color_idx];
      }
      edge_idx++;
    });

    COSMO_LOG(info) << "size of DBG: " << size_in_mega_bytes(dbg) << " MB";
    sdsl::store_to_file(dbg, params.output_prefix + params.output_base + ".dbg");

//    rrr_vector<63> color_rrr(color_bv);
    //sd_vector<> color_sd(color_bv);
    size_t total_colors = edge_idx * num_colors;
    COSMO_LOG(info) << "Color density : " << num_set/(double)total_colors * 100 << "%";
    COSMO_LOG(info) << "Set bits : " << num_set << std::endl;
    COSMO_LOG(info) << "total bits : " << total_colors << std::endl;
//    COSMO_LOG(info) << "size of color_bv  : " << size_in_mega_bytes(color_bv) << " MB";
//    COSMO_LOG(info) << "size of color_rrr : " << size_in_mega_bytes(color_rrr) << " MB";
    //COSMO_LOG(info) << "size of color_sd  : " << size_in_mega_bytes(color_sd) << " MB";
//    sdsl::store_to_file(color_rrr, params.output_prefix + params.output_base + ".rrr");
    cfs.close();
  }
    
  else {
    COSMO_LOG(error) << "Unsupported operation.";
    exit(EXIT_FAILURE);
  }

  COSMO_LOG(trace) << "Done!";
  COSMO_LOG(info) << endl
    << (stxxl::stats_data(*stats) - stats_begin)
    << " Peak disk allocs                           : " << bm->get_maximum_allocation()/1048576.0
    << " MB";
  int total_time = getMilliSpan(start_time);
      COSMO_LOG(info) << "total time since start:" << total_time << std::endl;
      COSMO_LOG(info) << "merge fraction of total time" << (float)merge_done_time / (float)total_time << std::endl;
  return 0;
}
